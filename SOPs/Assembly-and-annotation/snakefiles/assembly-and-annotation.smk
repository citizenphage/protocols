import pandas as pd
samples = pd.read_csv('resources/metadata.tbl', index_col=0, comment='#')

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
		expand("output/{sample}/pharokka/mapped-reads.bam", sample=samples.index),
		expand("output/{sample}/pharokka/{sample}.gbk", sample=samples.index)

rule download_reads:
	output:
		fwd="scratch/{sample}/fwd.fq.gz",
		rev="scratch/{sample}/rev.fq.gz"
	params:
		fwd=get_fwd_reads,
		rev=get_rev_reads
	shell:
		"""
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.fwd} {output.fwd}
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.rev} {output.rev}
		"""

rule download_host:
	output:
		temp("scratch/{sample}/host.fa")
	params:
		host=get_host
	shell:
		"""
			rsync -avh --progress  bt273@login.isca.ex.ac.uk:{params.host} {output}
		"""

rule fastqc_reads:
	input:
		fwd=rules.download_reads.output.fwd,
		rev=rules.download_reads.output.rev
	output:
		fwd="output/{sample}/fastqc-fwd.html",
		rev="output/{sample}/fastqc-rev.html"
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
		fwd=temp("scratch/{sample}/fwd.qc.fq.gz"),
		rev=temp("scratch/{sample}/rev.qc.fq.gz"),
		report="output/{sample}/qc-report.json"
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

rule map_reads_against_host:
	input:
		fwd=rules.qc_reads.output.fwd,
		rev=rules.qc_reads.output.rev,
		host=rules.download_host.output
	output:
		mapped="output/{sample}/host-mapped-reads.bam",
		unmapped=temp("scratch/{sample}/unmapped-reads.bam"),
		mapped_fwd=temp("scratch/{sample}/host-mapped-reads-fwd.fq.gz"),
		mapped_rev=temp("scratch/{sample}/host-mapped-reads-rev.fq.gz"),
		mapped_single=temp("scratch/{sample}/mapped-reads-single.fq.gz"),
		unmapped_fwd=temp("scratch/{sample}/unmapped-reads-fwd.fq.gz"),
		unmapped_rev=temp("scratch/{sample}/unmapped-reads-rev.fq.gz"),
		unmapped_single=temp("scratch/{sample}/unmapped-reads-single.fq.gz")

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


rule assemble_reads:
	input:
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev
	output:
		contigs="output/{sample}/assembly/contigs.fa.gz",
		graph="output/{sample}/assembly/{sample}-contigs.gfa.gz"
	conda:
		"../envs/assembly.yml"
	log: "logs/{sample}/assembly.log"
	threads: 16
	shell:
		"""
		shovill \
		--R1 {input.fwd} \
		--R2 {input.rev} \
		--outdir scratch/{wildcards.sample}/shovill \
		--minlen 10000 \
		--mincov 20 \
		--cpus {threads} \
		--noreadcorr 2>&1 | tee {log}

		pigz -c scratch/{wildcards.sample}/shovill/contigs.fa > {output.contigs}
		pigz -c scratch/{wildcards.sample}/shovill/contigs.gfa > {output.graph}

		rm -rf scratch/{wildcards.sample}/shovill
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

rule run_checkv:
	input:
		assembly=rules.assemble_reads.output.contigs,
		db=rules.setup_checkv.output
	output:
		viruses=temp("scratch/{sample}/checkv/viruses.fna"),
		quality_summary="output/{sample}/checkv/quality_summary.tsv"
	conda:
		"../envs/checkv.yml"
	threads: 16
	log: "logs/{sample}/checkv.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/checkv

		checkv end_to_end {input.assembly} \
		scratch/{wildcards.sample}/checkv \
		-t {threads} \
		-d {input.db}/checkv-db-v* \
		--remove_tmp 2>&1 | tee {log}

		mv scratch/{wildcards.sample}/checkv/quality_summary.tsv output/{wildcards.sample}/checkv/quality_summary.tsv
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


rule annotate_checkv_output:
	input:
		contigs=rules.run_checkv.output.viruses,
		db=rules.setup_pharokka.output
	output:
		temp("scratch/{sample}/pharokka/{sample}.gff")
	conda:
		"../envs/pharokka.yml"
	threads: 16
	log: "logs/{sample}/pharokka.log"
	shell:
		"""
		pharokka.py \
		-i {input.contigs} \
		-o scratch/{wildcards.sample}/pharokka \
		-d {input.db} \
		-t {threads} \
		--force \
		--prefix {wildcards.sample} 2>&1 | tee {log}

		"""

rule reorientate_pharokka_output:
	input:
		contigs=rules.run_checkv.output.viruses,
		db=rules.setup_pharokka.output,
		gff=rules.annotate_checkv_output.output
	output:
		gbk="output/{sample}/pharokka/{sample}.gbk",
		genome_map="output/{sample}/pharokka/{sample}.png",
		genome="output/{sample}/pharokka/{sample}.fa",
		inphared="output/{sample}/pharokka/closest-inphared-hit.tsv"
	conda:
		"../envs/pharokka.yml"
	params:
		title=get_title
	threads: 16
	log: "logs/{sample}/pharokka-rd2.log"
	shell:
		"""

		TERM_START=$(cat {input.gff} | grep "terminase large subunit" | cut -f 4)
		STRAND=$(cat {input.gff} | grep "terminase large subunit" | cut -f 7)

		if [ "$STRAND" == "+" ]; then
		STRAND='pos'
		else
		STRAND='neg'
		fi

		pharokka.py \
		-i {input.contigs} \
		-o scratch/{wildcards.sample}/pharokka-rd2 \
		-d {input.db} \
		-t {threads} \
		-l {wildcards.sample} \
		--force \
		--terminase \
		--terminase_strand $STRAND \
		--terminase_start $TERM_START \
		--prefix {wildcards.sample} 2>&1 | tee {log}

		cp scratch/{wildcards.sample}/pharokka-rd2/{wildcards.sample}.gbk {output.gbk}

		pharokka_plotter.py \
		-i {input.contigs} \
		--plot_name scratch/{wildcards.sample}/{wildcards.sample} \
		-p {wildcards.sample} \
		-o scratch/{wildcards.sample}/pharokka-rd2 \
		-t "{params.title}" 2>&1 | tee -a {log}

		mv scratch/{wildcards.sample}/{wildcards.sample}.png {output.genome_map}
		mv scratch/{wildcards.sample}/pharokka-rd2/{wildcards.sample}_genome_terminase_reoriented.fasta {output.genome}
		mv scratch/{wildcards.sample}/pharokka-rd2/{wildcards.sample}_top_hits_mash_inphared.tsv {output.inphared}




		"""

rule extract_fasta_from_genbank:
	input:
		assembly=rules.reorientate_pharokka_output.output.gbk
	output:
		"scratch/{sample}/{sample}.fa"
	shell:
		"""
		mkdir -p scratch/{wildcards.sample}
		python scripts/extract_fasta_from_genbank.py \
		--gbk {input.assembly} \
		--outfile scratch/{wildcards.sample}/{wildcards.sample}.fa

		"""

rule map_reads_against_assembly:
	input:
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev,
		assembly=rules.extract_fasta_from_genbank.output
	output:
		mapped="output/{sample}/pharokka/mapped-reads.bam"
	conda:
		"../envs/map-reads.yml"
	threads: 16
	shell:
		"""
		
		minimap2 \
		-a -x sr -t {threads} \
		-o scratch/{wildcards.sample}/assembly-mapping.sam \
		scratch/{wildcards.sample}/{wildcards.sample}.fa \
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
		"output/{sample}/pharokka/coverage.png"
	shell:
		"""
		python scripts/plot_genome_coverage.py \
		--bam {input.bam} \
		--prefix output/{wildcards.sample}/pharokka/coverage \
		--genome {input.genome} \
		"""