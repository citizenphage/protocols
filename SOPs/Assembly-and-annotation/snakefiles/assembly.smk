import pandas as pd
import json

with open('resources/config.json', 'r') as handle:
	config = json.load(handle)

samples = pd.read_csv(config['metadata_file'], index_col=0, comment='#')


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


def get_read_server(wildcards):
	return config['read_store_server']


def get_sequence_date(wildcards):
	row = samples.loc[wildcards.sample]
	return f'{row.date_sequenced}'

def get_sequencing_centre(wildcards):
	row = samples.loc[wildcards.sample]
	return f'{row.sequencing_centre}'

def get_sequencing_type(wildcards):
	row = samples.loc[wildcards.sample]
	return f'{row.sequencing_type}'


rule all:
	input:
		expand("output/{sample}/01_reads/orig_reads/stats.json", sample=samples.index),
		expand("output/{sample}/01_reads/orig_reads/qc-report.tgz", sample=samples.index),
		expand("output/{sample}/01_reads/host-mapping/host_mapping.json", sample=samples.index),
		expand("output/{sample}/02_assembly/unicycler/unicycler-subsampled-contigs-report.json", sample=samples.index),
		expand("output/{sample}/02_assembly/unicycler/checkv/process-files-report.json", sample=samples.index)



rule download_reads:
	output:
		fwd=temp("scratch/{sample}/fwd.fq.gz"),
		rev=temp("scratch/{sample}/rev.fq.gz")
	params:
		fwd=get_fwd_reads,
		rev=get_rev_reads,
		store=get_read_server
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/
			rsync -avh --progress  {params.store}:{params.fwd} {output.fwd}
			rsync -avh --progress  {params.store}:{params.rev} {output.rev}
		"""

rule download_host:
	output:
		"output/{sample}/01_reads/host-mapping/host.fa"
	params:
		host=get_host,
		store=get_read_server
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/
			rsync -avh --progress  {params.store}:{params.host} {output}
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
		fwd="output/reads/{sample}/orig_reads/{sample}-fwd.qc.fq.gz",
		rev="output/reads/{sample}/orig_reads/{sample}-rev.qc.fq.gz",
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

rule parse_qc_output:
	input:
		report=rules.qc_reads.output.report
	output:
		"output/{sample}/01_reads/orig_reads/stats.json"
	params:
		date=get_sequence_date,
		centre=get_sequencing_centre,
		type=get_sequencing_type
	shell:
		"""
			python scripts/parse_qc_output.py \
			--report {input.report} \
			--output {output} \
			--date {params.date} \
			--centre "{params.centre}" \
			--type "{params.type}"
			
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
		unmapped=temp("scratch/{sample}/unmapped-reads.bam"),
		mapped_fwd="output/reads/{sample}/host-mapping/host-mapped-reads-fwd.fq.gz",
		mapped_rev="output/reads/{sample}/host-mapping/host-mapped-reads-rev.fq.gz",
		mapped_single=temp("scratch/{sample}/mapped-reads-single.fq.gz"),
		unmapped_fwd="output/reads/{sample}/host-mapping/unmapped-reads-fwd.fq.gz",
		unmapped_rev="output/reads/{sample}/host-mapping/unmapped-reads-rev.fq.gz",
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
		"../envs/unicycler.yml"
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
		viruses="scratch/{sample}/02_assembly/unicycler/checkv/viruses.fna",
		quality_summary="scratch/{sample}/02_assembly/unicycler/checkv/quality_summary.tsv"
	conda:
		"../envs/checkv.yml"
	threads: 16
	log: "logs/{sample}/checkv.log"
	shell:
		"""
		checkv end_to_end {input.assembly} \
		scratch/{wildcards.sample}/checkv \
		-t {threads} \
		-d {input.db}/checkv-db-v* \
		--remove_tmp 2>&1 | tee {log}

		mv scratch/{wildcards.sample}/checkv/quality_summary.tsv {output.quality_summary}
		mv scratch/{wildcards.sample}/checkv/viruses.fna {output.viruses}
		"""

rule process_checkv_output:
	input:
		viruses=rules.run_checkv_on_unicycler.output.viruses,
		quality_summary=rules.run_checkv_on_unicycler.output.quality_summary
	output:
		report="output/{sample}/02_assembly/unicycler/checkv/process-files-report.json",
		store=directory("output/{sample}/02_assembly/unicycler/checkv/contig-store")
	shell:
		"""
			python scripts/parse_checkv.py \
			--viruses {input.viruses} \
			--quality {input.quality_summary} \
			--output {output.report} \
			--contig_store {output.store}
		"""
