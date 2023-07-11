import pandas as pd
samples = pd.read_csv('../resources/metadata.tbl', index_col=0, comment='#')

def get_reads(wildcards):
    row = samples.loc[wildcards.sample]
    rtnValue = [f'{row.short_fwd}',
            f'{row.short_rev}']
    return rtnValue

def get_reads(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.host}'

rule all:
	input:
		expand("output/{phage}/fastqc-fwd.html", phage=samples.index)


rule download_reads:
	input:
		get_reads
	output:
		fwd="scratch/{phage}/fwd.fq.gz",
		rev="scratch/{phage}/rev.fq.gz"
	log:
	benchmark:
	shell:
		"""
			cp {input[0]} {output.fwd}
			cp {input[1]} {output.rev}
		"""

rule download_host:
	input:
		get_host
	output:
		"scratch/{phage}/host.fa"
	shell:
		"""
			cp {input} {output}
		"""

rule fastqc_reads:
	input:
		fwd=rules.download_reads.fwd,
		rev=rules.download_reads.rev
	output:
		fwd="output/{phage}/fastqc-fwd.html",
		rev="output/{phage}/fastqc-rev.html"
	conda:
		"../envs/qc.yml"
	log:
		"../logs/{phage}/fastqc.log"
	benchmark:
	threads: 16
	shell:
		"""
			mkdir -p scratch/{wildcards.phage}/fastqc-out
			fastqc \
			--threads {threads} \
			--outdir scratch/{wildcards.phage}/fastqc-out \
			{input.fwd} \
			{input.rev} 2>&1 | tee {log}

			mv scratch/{wildcards.phage}/fastqc-out/fwd_fastqc.html {output.fwd}
			mv scratch/{wildcards.phage}/fastqc-out/rev_fastqc.html {output.rev}

			rm -rf scratch/{wildcards.phage}/fastqc-out
		"""


rule qc_reads:
	input:
		fwd=rules.download_reads.fwd,
		rev=rules.download_reads.rev
	output:
		fwd="scratch/{phage}/fwd.qc.fq.gz",
		rev="scratch/{phage}/rev.qc.fq.gz",
		merged="scratch/{phage}/merged.qc.fq.gz",
		report="scratch/{phage}/qc-report.json"
	conda:
		"../envs/qc.yml"
	threads: 16
	log:
	benchmark:
	shell:
		"""
			fastp \
			--in1 {input.fwd} \
			--in2 {input.rev} \
			--out1 {output.fwd} \
			--out2 {output.rev} \
			--merge \
			--merged_out {output.merged} \
			--dedup \
			--dup_calc_accuracy 6 \
			--length_required 30 \
			--correction \
			--json {output.report} \
			--thread {threads} 2>&1 | tee {log}


		"""

rule assemble_reads:
	input:
		fwd=rules.qc_reads.fwd,
		rev=rules.qc_reads.rev
	output:
		contigs="output/{phage}/assembly/contigs.fa"
	conda:
		"../envs/assembly.yml"
	log:
	benchmark:
	shell:
		"""
		"""

