import pandas as pd
samples = pd.read_csv('resources/2024-06-18-SAP-metadata.csv', index_col=0, comment='#')


def get_fwd_reads(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.short_fwd}'

def get_rev_reads(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.short_rev}'

def get_host(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.host}'

def get_title(wildcards):
	row = samples.loc[wildcards.sample]
	return f'{row.host_genus} phage {row.given_name if row.given_name else ""} {wildcards.sample}'

rule all:
	input:
		expand("output/ncbi/reads/{sample}-fwd.qc.fq.gz", sample=samples.index),
		expand("output/{sample}/host-mapped-reads.bam", sample=samples.index),
		expand("output/{sample}/host_mapping.json", sample=samples.index),
		expand("output/{sample}/assembly/{sample}-shovill-contigs.gfa.gz", sample=samples.index),
		expand("output/{sample}/assembly/{sample}-unicycler-subsample.gfa.gz", sample=samples.index),
		#expand("output/{sample}/assembly/{sample}-unicycler-unmapped.gfa.gz", sample=samples.index),
		#expand("output/{sample}/checkv/unicycler-quality_summary.tsv", sample=samples.index),
		#expand("output/{sample}/checkv/shovill-quality_summary.tsv", sample=samples.index),
		#expand("output/{sample}/pharokka/mapped-reads.bam", sample=samples.index),
		#expand("output/{sample}/pharokka/unicycler/{sample}.gbk", sample=samples.index),
		#expand("output/{sample}/pharokka/coverage.json", sample=samples.index),

		#"output/viridic"

rule download_reads:
	output:
		fwd="scratch/{sample}/fwd.fq.gz",
		rev="scratch/{sample}/rev.fq.gz"
	params:
		fwd=get_fwd_reads,
		rev=get_rev_reads
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.fwd} {output.fwd}
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.rev} {output.rev}
		"""

rule download_host:
	output:
		"output/{sample}/01_reads/host-mapping/host.fa"
	params:
		host=get_host
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.host} {output}
		"""

rule fastqc_reads:
	input:
		fwd=rules.download_reads.output.fwd,
		rev=rules.download_reads.output.rev
	output:
		fwd="scratch/{sample}/01_reads/orig_reads/fastqc-fwd.html",
		rev="scratch/{sample}/01_reads/orig_reads/fastqc-rev.html"
	conda:
		"../envs/qc.yml"
	log:
		"../logs/{sample}/fastqc.log"
	threads: 16
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/fastqc-out
			fastqc \
			--threads {threads} \
			--outdir scratch/{wildcards.sample}/fastqc-out \
			{input.fwd} \
			{input.rev} 2>&1 | tee {log}

			mv scratch/{wildcards.sample}/fastqc-out/fwd_fastqc.html {output.fwd}
			mv scratch/{wildcards.sample}/fastqc-out/rev_fastqc.html {output.rev}

			rm -rf scratch/{wildcards.sample}/fastqc-out
		"""


rule qc_reads:
	input:
		fwd=rules.download_reads.output.fwd,
		rev=rules.download_reads.output.rev
	output:
		fwd="output/{sample}/01_reads/orig_reads/{sample}-fwd.qc.fq.gz",
		rev="output/{sample}/01_reads/orig_reads/{sample}-rev.qc.fq.gz",
		report="scratch/{sample}/01_reads/orig_reads/qc-report.json"
	conda:
		"../envs/qc.yml"
	threads: 16
	log: "logs/{sample}/qc-reads.log"
	shell:
		"""
			fastp \
			--in1 {input.fwd} \
			--in2 {input.rev} \
			--out1 {output.fwd} \
			--out2 {output.rev} \
			--dedup \
			--dup_calc_accuracy 6 \
			--length_required 30 \
			--correction \
			--json {output.report} \
			--thread {threads} 2>&1 | tee {log}


		"""

rule archive_qc_data:
	input:
		fwd_qc=rules.fastqc_reads.output.fwd,
		rev_qc=rules.fastqc_reads.output.rev,
		report=rules.qc_reads.output.report
	output:
		"output/{sample}/01_reads/orig_reads/qc-report.tgz"
	shell:
		"""
		tar czvf {output} {input.fwd_qc} {input.rev_qc} {input.report}
		"""



rule map_reads_against_host:
	input:
		fwd=rules.qc_reads.output.fwd,
		rev=rules.qc_reads.output.rev,
		host=rules.download_host.output
	output:
		mapped="output/{sample}/01_reads/host-mapping/host-mapped-reads.bam",
		unmapped="scratch/{sample}/unmapped-reads.bam",
		mapped_fwd="output/{sample}/01_reads/host-mapping/host-mapped-reads-fwd.fq.gz",
		mapped_rev="output/{sample}/01_reads/host-mapping/host-mapped-reads-rev.fq.gz",
		mapped_single="scratch/{sample}/mapped-reads-single.fq.gz",
		unmapped_fwd="output/{sample}/01_reads/host-mapping/unmapped-reads-fwd.fq.gz",
		unmapped_rev="output/{sample}/01_reads/host-mapping/unmapped-reads-rev.fq.gz",
		unmapped_single="scratch/{sample}/unmapped-reads-single.fq.gz"

	conda:
		"../envs/map-reads.yml"
	threads: 16
	log: "logs/{sample}/map-reads.log"
	shell:
		"""
		minimap2 \
		-a -x sr -t {threads} \
		-o scratch/{wildcards.sample}/mapping.sam \
		{input.host} \
		{input.fwd} \
		{input.rev} 2>&1 | tee {log}

		samtools view -bS -F 4 \
		-@ {threads} scratch/{wildcards.sample}/mapping.sam | 
		samtools sort -@ {threads} -o {output.mapped}

		samtools view -bS -f 4 \
		-@ {threads} scratch/{wildcards.sample}/mapping.sam | 
		samtools sort -@ {threads} -o {output.unmapped}

		samtools fastq \
		-0 /dev/null \
		-1 {output.mapped_fwd} \
		-2 {output.mapped_rev} \
		-s {output.mapped_single} \
		-@ {threads} \
		-n {output.mapped}

		samtools fastq \
		-0 /dev/null \
		-1 {output.unmapped_fwd} \
		-2 {output.unmapped_rev} \
		-s {output.unmapped_single} \
		-@ {threads} \
		-n {output.unmapped}

		rm scratch/{wildcards.sample}/mapping.sam
		"""

rule parse_mapped_reads:
	input:
		mapped_fwd=rules.map_reads_against_host.output.mapped_fwd,
		mapped_rev=rules.map_reads_against_host.output.mapped_rev,
		unmapped_fwd=rules.map_reads_against_host.output.unmapped_fwd,
		unmapped_rev=rules.map_reads_against_host.output.unmapped_rev
	output:
		"output/{sample}/01_reads/host-mapping/host_mapping.json"
	shell:
		"""
			python scripts/parse_mapped_reads.py \
			--mapped_fwd {input.mapped_fwd} \
			--mapped_rev {input.mapped_rev} \
			--unmapped_fwd {input.unmapped_fwd} \
			--unmapped_rev {input.unmapped_rev} \
			--output {output}
		"""


rule assemble_reads_with_shovill:
	input:
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev
	output:
		contigs="output/{sample}/02_assembly/shovill/{sample}-shovill-contigs.fa.gz",
		graph="output/{sample}/02_assembly/shovill/{sample}-shovill-contigs.gfa.gz",
		sub_fwd="output/{sample}/01_reads/subsample/shovill-subsample-fwd.fq.gz",
		sub_rev="output/{sample}/01_reads/subsample/shovill-subsample-rev.fq.gz",
		sub_merged="output/{sample}/01_reads/subsample/shovill-subsample-merged.fq.gz",
		report ="output/{sample}/02_assembly/shovill/shovill-subsampled-contigs-report.json",
	conda:
		"../envs/assembly.yml"
	log: "logs/{sample}/shovill-assembly.log"
	threads: 16
	shell:
		"""
		shovill \
		--R1 {input.fwd} \
		--R2 {input.rev} \
		--outdir scratch/{wildcards.sample}/shovill \
		--minlen 10000 \
		--depth 500 \
		--mincov 20 \
		--keepfiles \
		--cpus {threads} \
		--force \
		--noreadcorr 2>&1 | tee {log}

		python scripts/rename_contigs.py \
		--input scratch/{wildcards.sample}/shovill/contigs.fa \
		--outfile scratch/{wildcards.sample}/shovill/renamed.fa \
		--report {output.report} \
		--prefix {wildcards.sample} \
		--assembly_method shovill

		pigz -c scratch/{wildcards.sample}/shovill/renamed.fa > {output.contigs}
		pigz -c scratch/{wildcards.sample}/shovill/contigs.gfa > {output.graph}



		mv scratch/{wildcards.sample}/shovill/flash.notCombined_1.fastq.gz {output.sub_fwd}
		mv scratch/{wildcards.sample}/shovill/flash.notCombined_2.fastq.gz {output.sub_rev}
		mv scratch/{wildcards.sample}/shovill/flash.extendedFrags.fastq.gz {output.sub_merged}

		rm -rf scratch/{wildcards.sample}/shovill
		"""

rule reassemble_shovill_reads_with_unicycler:
	input:
		fwd=rules.assemble_reads_with_shovill.output.sub_fwd,
		rev=rules.assemble_reads_with_shovill.output.sub_rev,
		merged=rules.assemble_reads_with_shovill.output.sub_merged
	output:
		contigs="output/{sample}/02_assembly/unicycler/unicycler-subsampled-contigs.fa",
		report ="output/{sample}/02_assembly/unicycler/unicycler-subsampled-contigs-report.json",
		graph="output/{sample}/02_assembly/unicycler/{sample}-unicycler-subsample.gfa.gz"
	conda:
		"../envs/assembly.yml"
	log: "logs/{sample}/unicycler-subsample-assembly.log"
	threads: 16
	shell:
		"""
		unicycler \
		-1 {input.fwd} \
		-2 {input.rev} \
		-s {input.merged} \
		--min_fasta_length 10000 \
		--threads {threads} \
		--start_genes resources/terL.fa \
		--out scratch/{wildcards.sample}/unicycler-subsample 2>&1 | tee {log}

		python scripts/rename_contigs.py \
		--input scratch/{wildcards.sample}/unicycler-subsample/assembly.fasta \
		--outfile scratch/{wildcards.sample}/unicycler-subsample/renamed.fa \
		--report {output.report} \
		--prefix {wildcards.sample} \
		--assembly_method unicycler-subsample

		cp scratch/{wildcards.sample}/unicycler-subsample/renamed.fa {output.contigs}
		pigz -c scratch/{wildcards.sample}/unicycler-subsample/assembly.gfa > {output.graph}


		"""


rule assemble_mapped_reads_with_unicycler:
	input:
		fwd=rules.map_reads_against_host.output.mapped_fwd,
		rev=rules.map_reads_against_host.output.mapped_rev,
	output:
		contigs="output/{sample}/02_assembly/host-mapped-reads/unicycler-host-mapped-contigs.fa",
		graph="output/{sample}/02_assembly/host-mapped-reads/unicycler-host-mapped.gfa.gz",
		report ="output/{sample}/02_assembly/host-mapped-reads/unicycler-host-mapped-contigs-report.json"
	conda:
		"../envs/assembly.yml"
	log: "logs/{sample}/unicycler-host-mapped-assembly.log"
	threads: 16
	shell:
		"""
		unicycler \
		-1 {input.fwd} \
		-2 {input.rev} \
		--min_fasta_length 1000 \
		--threads {threads} \
		--out scratch/{wildcards.sample}/unicycler-host-mapped 2>&1 | tee {log}

		python scripts/rename_contigs.py \
		--input scratch/{wildcards.sample}/unicycler-host-mapped/assembly.fasta \
		--outfile scratch/{wildcards.sample}/unicycler-host-mapped/renamed.fa \
		--report {output.report} \
		--prefix {wildcards.sample} \
		--assembly_method unicycler-all

		cp scratch/{wildcards.sample}/unicycler-host-mapped/assembly.fasta {output.contigs}
		pigz -c scratch/{wildcards.sample}/unicycler-host-mapped/assembly.gfa > {output.graph}

		rm -rf scratch/{wildcards.sample}/unicycler-host-mapped

		"""




rule setup_checkv:
	output:
		directory("scratch/checkv-db")
	conda: "../envs/checkv.yml"
	log: "logs/setup-checkv.log"
	shell:
		"""
		checkv download_database {output} 2>&1 | tee {log}
		"""

rule run_checkv_on_unicycler:
	input:
		assembly=rules.reassemble_shovill_reads_with_unicycler.output.contigs,
		db=rules.setup_checkv.output
	output:
		viruses="scratch/{sample}/checkv/unicycler-viruses.fna",
		quality_summary="output/{sample}/02_assembly/unicycler/checkv-quality_summary.tsv"
	conda:
		"../envs/checkv.yml"
	threads: 16
	log: "logs/{sample}/checkv.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/unicycler-checkv

		checkv end_to_end {input.assembly} \
		scratch/{wildcards.sample}/unicycler-checkv \
		-t {threads} \
		-d {input.db}/checkv-db-v* \
		--remove_tmp 2>&1 | tee {log}

		mv scratch/{wildcards.sample}/unicycler-checkv/quality_summary.tsv {output.quality_summary}
		mv scratch/{wildcards.sample}/unicycler-checkv/viruses.fna {output.viruses}
		"""

rule setup_pharokka:
	output:
		directory("scratch/pharokka-db")
	conda: "../envs/pharokka.yml"
	log: "logs/setup-pharokka.log"
	shell:
		"""
		install_databases.py -o {output} 2>&1 | tee {log}
		"""


rule annotate_unicycler_output:
	input:
		contigs=rules.reassemble_shovill_reads_with_unicycler.output.contigs,
		db=rules.setup_pharokka.output
	output:
		gbk="output/{sample}/03_annotation/{sample}-pharokka.gbk",
		inphared="output/{sample}/03_annotation/{sample}_top_hits_mash_inphared.tsv"
	conda:
		"../envs/pharokka.yml"
	params:
		title=get_title
	threads: 16
	log: "logs/{sample}/pharokka-unicycler.log"
	shell:
		"""
			pharokka.py \
			-i {input.contigs} \
			-o scratch/{wildcards.sample}/pharokka-unicycler \
			-d {input.db} \
			-t {threads} \
			-l {wildcards.sample} \
			-g prodigal \
			--force \
			--dnaapler \
			--prefix {wildcards.sample} 2>&1 | tee {log}

			cp scratch/{wildcards.sample}/pharokka-unicycler/{wildcards.sample}.gbk {output.gbk}
			cp scratch/{wildcards.sample}/pharokka-unicycler/{wildcards.sample}_top_hits_mash_inphared.tsv {output.inphared}

	
		"""

rule setup_phold:
	output:
		directory("scratch/phold-db")
	conda: "../envs/phold.yml"
	log: "logs/setup-phold.log"
	shell:
		"""
		phold install -d {output} 2>&1 | tee {log}
		"""

rule run_phold:
	input:
		gbk=rules.annotate_genome.output.gbk,
		db=rules.setup_phold.output
	output:
		plot_png="output/{sample}/03_annotation/{sample}-phold.png",
		plot_svg="output/{sample}/03_annotation/{sample}-phold.svg",
		gbk="output/{sample}/03_annotation/{sample}-phold.gbk"
	conda:
		"../envs/phold.yml"
	threads: 16
	log: "logs/{sample}/pharokka-phold.log"
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/pharokka-unicycler/phold
			phold run -i {input.gbk} \
			-p {wildcards.sample} \
			-o scratch/{wildcards.sample}/pharokka-unicycler/phold/{wildcards.sample} \
			-t {threads} \
			--database {input.db} \
			--force

			cp scratch/{wildcards.sample}/pharokka-unicycler/phold/{wildcards.sample}/{wildcards.sample}.gbk {output.gbk}

			phold plot -i {output.gbk} \
			-p {wildcards.sample} \
			-o scratch/{wildcards.sample}/pharokka-unicycler/phold/plots \
			-t {wildcards.sample} \
			--force

			cp scratch/{wildcards.sample}/pharokka-unicycler/phold/plots/{wildcards.sample}.svg {output.plot_svg}
			cp scratch/{wildcards.sample}/pharokka-unicycler/phold/plots/{wildcards.sample}.png {output.plot_png}


		"""


rule get_nearest_ncbi_hit:
	input:
		rules.annotate_genome.output.inphared
	output:
		genome="output/{sample}/03_annotation/nearest-neighbour.fa",
		report="output/{sample}/03_annotation/nearest-neighbour-report.json"
	shell:
		"""
			python scripts/get_closest_hit_from_inphared.py \
			--inphared_file {input} \
			--report {output.report} \
			--outfile {output.genome}
		"""



rule extract_fasta_from_genbank:
	input:
		assembly=rules.run_phold.output.gbk
	output:
		"output/{sample}/03_annotation/{sample}.fa"
	shell:
		"""
		mkdir -p scratch/{wildcards.sample}
		python ../Assembly-and-annotation/scripts/extract_fasta_from_genbank.py \
		--gbk TempertonLab-PAO1-2024-04-30.gbk \
		--outfile TempertonLab-PAO1-2024-04-30.fa \
		--name TempertonLab-PAO1-2024-04-30

		"""

rule map_reads_against_assembly:
	input:
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev,
		assembly=rules.extract_fasta_from_genbank.output
	output:
		mapped="output/{sample}/03_annotation/assembly-mapped-reads.bam"
	conda:
		"../envs/map-reads.yml"
	threads: 16
	shell:
		"""
		
		minimap2 \
		-a -x sr -t {threads} \
		-o scratch/{wildcards.sample}/assembly-mapping.sam \
		{input.assembly} \
		{input.fwd} \
		{input.rev} 2>&1 | tee {log}

		samtools view -bS -F 4 \
		-@ {threads} scratch/{wildcards.sample}/assembly-mapping.sam | 
		samtools sort -@ {threads} -o {output.mapped}

		samtools index {output.mapped}


		rm scratch/{wildcards.sample}/assembly-mapping.sam
		"""

rule plot_coverage:
	input:
		bam=rules.map_reads_against_assembly.output.mapped,
		genome=rules.extract_fasta_from_genbank.output
	output:
		json="output/{sample}/03_annotation/coverage.json"
	shell:
		"""
		python scripts/plot_genome_coverage.py \
		--bam {input.bam} \
		--prefix output/{wildcards.sample}/03_annotation/coverage \
		--genome {input.genome} \
		--json {output.json}
		"""

rule check_for_variants:
	input:
		gbk=rules.run_phold.output.gbk,
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev,
	output:
		report="output/{sample}/04_variants/{sample}-snps.html",
		bam="output/{sample}/04_variants/{sample}-snps.bam"
	conda:
		"../envs/snippy.yml"
	threads: 16
	shell:
		"""
			snippy --cpus {threads} --outdir scratch/{wildcards.sample}/snippy/ \
			--ref {input.gbk} \
			--R1 {input.fwd} \
			--R2 {input.rev} \
			--force 

			cp scratch/{wildcards.sample}/snippy/snps.html {output.report}
			cp scratch/{wildcards.sample}/snippy/snps.bam {output.bam}
		"""

rule stage_data:
	input:
		fwd=rules.qc_reads.output.fwd,
		rev=rules.qc_reads.output.rev,
		gbk=rules.annotate_genome.output.gbk,
	output:
		gbk="output/ncbi/GOSH-paper/genomes/{sample}/{sample}.gbk",
		fwd="output/ncbi/GOSH-paper/reads/{sample}/{sample}.fwd.qc.fq.gz",
		rev="output/ncbi/GOSH-paper/reads/{sample}/{sample}.rev.qc.fq.gz"
	shell:
		"""
			cp {input.fwd} {output.fwd}
			cp {input.rev} {output.rev}
			cp {input.gbk} {output.gbk}
		"""

rule generate_report:
	input:
		unicycler_log="logs/{sample}/unicycler-subsample-assembly.log",
		shovill_contigs="output/{sample}/assembly/{sample}-shovill-contigs.fa.gz",
		checkv_log="output/{sample}/checkv/shovill-quality_summary.tsv",
		inphared_log="output/{sample}/pharokka/unicycler/closest-inphared-hit.tsv",
		coverage_log="output/{sample}/pharokka/coverage.json",
		pharokka_log="logs/{sample}/pharokka-unicycler.log"

	output:
		"output/{sample}/{sample}-report.json"
	shell:
		"""
		python scripts/generate_report.py \
		--output {output} \
		--unicycler_log {input.unicycler_log} \
		--shovill_contigs {input.shovill_contigs} \
		--checkv_log {input.checkv_log} \
		--inphared_log {input.inphared_log} \
		--coverage_log {input.coverage_log} \
		--pharokka_log {input.pharokka_log}
		"""

