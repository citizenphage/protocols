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
    return f'{row.host_fasta}'

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
		expand("output/reads/{sample}/orig_reads/{sample}-fwd.qc.fq.gz", sample=samples.index),
		expand("output/reads/{sample}/orig_reads/{sample}-rev.qc.fq.gz", sample=samples.index),
		expand("output/{sample}/01_reads/host-mapping/host.fa", sample=samples.index),
		expand("output/{sample}/01_reads/orig_reads/stats.json", sample=samples.index),
		expand("output/{sample}/01_reads/orig_reads/qc-report.tgz", sample=samples.index),
		expand("output/{sample}/01_reads/host-mapping/host_mapping.json", sample=samples.index),
		expand("output/{sample}/02_assembly/unicycler/shovill-reads/{sample}-shovill-bandage.png", sample=samples.index),
		expand("output/{sample}/02_assembly/unicycler/1pc/{sample}-1pc-bandage.png", sample=samples.index),
		expand("output/reports/{sample}/assembly-report.json", sample=samples.index),
		expand("output/{sample}/01_reads/orig_reads/stats-processed.json", sample=samples.index),
		expand("output/{sample}/01_reads/host-mapping/host_mapping-processed.json", sample=samples.index),
		expand("output/{sample}/02_assembly/unicycler/shovill-reads/report-processed.json", sample=samples.index)



rule download_reads:
	output:
		fwd=temp("scratch/{sample}/all-fwd.fq.gz"),
		rev=temp("scratch/{sample}/all-rev.fq.gz")
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

			mv scratch/{wildcards.sample}/fastqc-out/all-fwd_fastqc.html {output.fwd}
			mv scratch/{wildcards.sample}/fastqc-out/all-rev_fastqc.html {output.rev}

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

rule update_db_reads:
	input:
		report=rules.parse_qc_output.output,
		fwd_qc=rules.fastqc_reads.output.fwd,
		rev_qc=rules.fastqc_reads.output.rev
	output:
		"output/{sample}/01_reads/orig_reads/stats-processed.json"
	shell:
		"""
		python scripts/data/update_db_reads.py \
		--report {input.report} \
		--sample {wildcards.sample} \
		--fwd {input.fwd_qc} \
		--rev {input.rev_qc} \
		--output {output}
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


rule update_mapped_reads:
	input:
		report= rules.parse_mapped_reads.output,
		bam_file = rules.map_reads_against_host.output.mapped,
	output:
		"output/{sample}/01_reads/host-mapping/host_mapping-processed.json"
	shell:
		"""
			python scripts/data/update_db_mapped_reads.py \
			--report {input.report} \
			--sample {wildcards.sample} \
			--bam_file {input.bam_file} \
			--output {output}
		"""


rule subsample_reads:
	input:
		fwd=rules.map_reads_against_host.output.unmapped_fwd,
		rev=rules.map_reads_against_host.output.unmapped_rev
	output:
		sub_fwd = "scratch/{sample}/1pc/fwd.fq.gz",
		sub_rev = "scratch/{sample}/1pc/rev.fq.gz"
	conda:
		"../envs/assembly.yml"
	shell:
		"""
			seqtk sample {input.fwd} 0.01  | pigz --fast -c -p 16 > {output.sub_fwd}
			seqtk sample {input.rev} 0.01  | pigz --fast -c -p 16 > {output.sub_rev}
		"""

rule assemble_subsample:
	input:
		fwd= "scratch/{sample}/1pc/fwd.fq.gz",
		rev="scratch/{sample}/1pc/rev.fq.gz"
	output:
		contigs = "output/{sample}/02_assembly/unicycler/1pc/contigs.fa",
		report="output/{sample}/02_assembly/unicycler/1pc/report.json",
		graph = "output/{sample}/02_assembly/unicycler/1pc/bandage.gfa.gz"
	conda:
		"../envs/unicycler.yml"
	threads: 16
	shell:
		"""
		unicycler \
		-1 {input.fwd} \
		-2 {input.rev} \
		--min_fasta_length 1000 \
		--threads {threads} \
		--out scratch/{wildcards.sample}/unicycler/1pc/
	
		python scripts/rename_contigs.py \
		--input scratch/{wildcards.sample}/unicycler/1pc/assembly.fasta \
		--outfile scratch/{wildcards.sample}/unicycler/1pc/renamed.fa \
		--report {output.report} \
		--prefix {wildcards.sample}-1pc \
		--assembly_method unicycler-1pc
	
		cp scratch/{wildcards.sample}/unicycler/1pc/renamed.fa {output.contigs}
		pigz -c scratch/{wildcards.sample}/unicycler/1pc/assembly.gfa > {output.graph}
	
	
		"""

rule visualise_subsample_assembly:
	input:
		graph=rules.assemble_subsample.output.graph
	output:
		img="output/{sample}/02_assembly/unicycler/1pc/{sample}-1pc-bandage.png"
	conda:
		"../envs/bandage.yml"
	shell:
		"""
			unpigz -c {input.graph} > scratch/{wildcards.sample}/unicycler/1pc/bandage.gfa
			Bandage image scratch/{wildcards.sample}/unicycler/1pc/bandage.gfa {output.img}
			rm scratch/{wildcards.sample}/unicycler/1pc/bandage.gfa

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
		contigs="output/{sample}/02_assembly/unicycler/shovill-reads/contigs.fa",
		report ="output/{sample}/02_assembly/unicycler/shovill-reads/report.json",
		graph="output/{sample}/02_assembly/unicycler/shovill-reads/bandage.gfa.gz"
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
		--min_fasta_length 1000 \
		--threads {threads} \
		--out scratch/{wildcards.sample}/unicycler-subsample 2>&1 | tee {log}

		python scripts/rename_contigs.py \
		--input scratch/{wildcards.sample}/unicycler-subsample/assembly.fasta \
		--outfile scratch/{wildcards.sample}/unicycler-subsample/renamed.fa \
		--report {output.report} \
		--prefix {wildcards.sample} \
		--assembly_method unicycler-shovill-reads

		cp scratch/{wildcards.sample}/unicycler-subsample/renamed.fa {output.contigs}
		pigz -c scratch/{wildcards.sample}/unicycler-subsample/assembly.gfa > {output.graph}


		"""

rule visualise_unicycler_assembly:
	input:
		graph=rules.reassemble_shovill_reads_with_unicycler.output.graph
	output:
		img="output/{sample}/02_assembly/unicycler/shovill-reads/{sample}-shovill-bandage.png"
	conda:
		"../envs/bandage.yml"
	shell:
		"""
			unpigz -c {input.graph} > scratch/{wildcards.sample}/unicycler-subsample/bandage.gfa
			Bandage image scratch/{wildcards.sample}/unicycler-subsample/bandage.gfa {output.img}
			rm scratch/{wildcards.sample}/unicycler-subsample/bandage.gfa
			
			mkdir -p output/{wildcards.sample}/03_selected_contigs/
			
		"""


rule update_db_unicycler_shovill:
	input:
		contigs=rules.reassemble_shovill_reads_with_unicycler.output.contigs,
		graph=rules.reassemble_shovill_reads_with_unicycler.output.graph,
		graph_img=rules.visualise_unicycler_assembly.output.img
	output:
		"output/{sample}/02_assembly/unicycler/shovill-reads/report-processed.json"
	shell:
		"""
			mkdir -p scratch/{wildcards.sample}/unicycler/shovill-reads/
			pigz -c {input.contigs} > scratch/{wildcards.sample}/unicycler/shovill-reads/contigs.fa.gz
			
			python scripts/data/update_db_assembly.py \
			--contigs scratch/{wildcards.sample}/unicycler/shovill-reads/contigs.fa.gz \
			--assembler unicycler \
			--assembler_version "0.5.0" \
			--command "--min_fasta_length 1000" \
			--graph {input.graph} \
			--graph_img {input.graph_img} \
			--subsampled_coverage shovill-reads \
			--sample {wildcards.sample} \
			--output {output}
			
			rm scratch/{wildcards.sample}/unicycler/shovill-reads/contigs.fa.gz
		"""


rule generate_assembly_report_json:
	input:
		qc_report="output/{sample}/01_reads/orig_reads/stats.json",
		mapping_report="output/{sample}/01_reads/host-mapping/host_mapping.json",
		shovill_report="output/{sample}/02_assembly/shovill/shovill-subsampled-contigs-report.json",
		unicycler_report="output/{sample}/02_assembly/unicycler/shovill-reads/report.json",
		onepct_report="output/{sample}/02_assembly/unicycler/1pc/report.json",
		qc_archive=rules.archive_qc_data.output,
		bandage=rules.visualise_subsample_assembly.output.img,
		bandage2=rules.visualise_unicycler_assembly.output.img

	output:
		"output/reports/{sample}/assembly-report.json"
	params:
		config='resources/config.json'
	shell:
		"""
		python scripts/generate_assembly_report.py \
		--output {output} \
		--config {params.config} \
		--sample {wildcards.sample} \
		--qc_report {input.qc_report} \
		--mapping_report {input.mapping_report} \
		--shovill_report {input.shovill_report} \
		--unicycler_report {input.unicycler_report} \
		--onepct_report {input.onepct_report}
		"""
