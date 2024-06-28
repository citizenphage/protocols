import pandas as pd
import json

with open('resources/config.json', 'r') as handle:
	config = json.load(handle)

samples = pd.read_csv(config['metadata_file'], index_col=0, comment='#')

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
		assembly="output/{sample}/03_selected_contigs/{contig}.fa",
		db=rules.setup_checkv.output
	output:
		quality_summary="scratch/checkv/quality_summary.tsv"
	conda:
		"../envs/checkv.yml"
	threads: 16
	log: "logs/{sample}/checkv.log"
	shell:
		"""
		checkv end_to_end {input.assembly} \
		scratch/{wildcards.sample}/{wildcards.contig}/checkv \
		-t {threads} \
		-d {input.db}/checkv-db-v* \
		--remove_tmp 2>&1 | tee {log}

		mv scratch/{wildcards.sample}/{wildcards.contig}/checkv/quality_summary.tsv {output.quality_summary}
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




rule setup_pharokka:
    output:
        directory("scratch/pharokka-db")
    conda: "../envs/pharokka.yml"
    log: "logs/setup-pharokka.log"
    shell:
        """
        install_databases.py -o {output} 2>&1 | tee {log}
        """
def get_title(wildcards):
	row = samples.loc[wildcards.sample]
	return f'{row.host_genus} phage {row.given_name if row.given_name else ""} {wildcards.sample}'

rule annotate_genome:
    input:
        contigs="scratch/{sample}/checkv/contig-store/{contig}.fa",
        db=rules.setup_pharokka.output
    output:
        gbk="output/{sample}/03_annotation/{contig}/pharokka.gbk",
        inphared="output/{sample}/03_annotation/{contig}/top_hits_mash_inphared.tsv"
    conda:
        "../envs/pharokka.yml"
    params:
        title=get_title
    threads: 16
    log: "logs/{sample}/{contig}/pharokka-unicycler.log"
    shell:
        """
            mkdir -p scratch/{wildcards.sample}/{wildcards.contig}
            pharokka.py \
            -i {input.contigs} \
            -o scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler \
            -d {input.db} \
            -t {threads} \
            -l {wildcards.contig} \
            -g prodigal \
            --force \
            --dnaapler \
            --prefix {wildcards.contig} 2>&1 | tee {log}

            cp scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/{wildcards.contig}.gbk {output.gbk}
            cp scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/{wildcards.contig}_top_hits_mash_inphared.tsv {output.inphared}


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
        plot_png="output/{sample}/03_annotation/{contig}/phold.png",
        plot_svg="output/{sample}/03_annotation/{contig}/phold.svg",
        gbk="output/{sample}/03_annotation/{contig}/phold.gbk"
    conda:
        "../envs/phold.yml"
    threads: 16
    log: "logs/{sample}/{contig}/pharokka-phold.log"
    shell:
        """
            mkdir -p scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold
            phold run -i {input.gbk} \
            -p {wildcards.contig} \
            -o scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold/{wildcards.contig} \
            -t {threads} \
            --database {input.db} \
            --force

            cp scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold/{wildcards.contig}/{wildcards.contig}.gbk {output.gbk}

            phold plot -i {output.gbk} \
            -p {wildcards.contig} \
            -o scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold/plots \
            -t {wildcards.contig} \
            --force

            cp scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold/plots/{wildcards.contig}.svg {output.plot_svg}
            cp scratch/{wildcards.sample}/{wildcards.contig}/pharokka-unicycler/phold/plots/{wildcards.contig}.png {output.plot_png}


        """


rule get_nearest_ncbi_hit:
    input:
        rules.annotate_genome.output.inphared
    output:
        genome="output/{sample}/03_annotation/{contig}/nearest-neighbour.fa",
        report="output/{sample}/03_annotation/{contig}/nearest-neighbour-report.json"
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
        "output/{sample}/03_annotation/{contig}/genome.fa"
    shell:
        """
        python scripts/extract_fasta_from_genbank.py \
        --gbk {input} \
        --outfile {output} \
        --name {wildcards.contig}

        """

rule map_reads_against_assembly:
    input:
        fwd="output/{sample}/01_reads/host-mapping/unmapped-reads-fwd.fq.gz",
        rev="output/{sample}/01_reads/host-mapping/unmapped-reads-rev.fq.gz",
        assembly=rules.extract_fasta_from_genbank.output
    output:
        mapped="output/{sample}/03_annotation/{contig}/assembly-mapped-reads.bam"
    conda:
        "../envs/map-reads.yml"
    threads: 16
    shell:
        """

        minimap2 \
        -a -x sr -t {threads} \
        -o scratch/{wildcards.sample}/{wildcards.contig}/assembly-mapping.sam \
        {input.assembly} \
        {input.fwd} \
        {input.rev}

        samtools view -bS -F 4 \
        -@ {threads} scratch/{wildcards.sample}/{wildcards.contig}/assembly-mapping.sam | 
        samtools sort -@ {threads} -o {output.mapped}

        samtools index {output.mapped}


        rm scratch/{wildcards.sample}/{wildcards.contig}/assembly-mapping.sam
        """

rule plot_coverage:
    input:
        bam=rules.map_reads_against_assembly.output.mapped,
        genome=rules.extract_fasta_from_genbank.output
    output:
        json="output/{sample}/03_annotation/{contig}/coverage.json"
    shell:
        """
        python scripts/plot_genome_coverage.py \
        --bam {input.bam} \
        --prefix output/{wildcards.sample}/03_annotation/{wildcards.contig}/coverage \
        --genome {input.genome} \
        --json {output.json}
        """

rule check_for_variants:
    input:
        gbk=rules.run_phold.output.gbk,
        fwd="output/{sample}/01_reads/host-mapping/unmapped-reads-fwd.fq.gz",
        rev="output/{sample}/01_reads/host-mapping/unmapped-reads-rev.fq.gz",
    output:
        report="output/{sample}/04_variants/{contig}/snps.html",
        bam="output/{sample}/04_variants/{contig}/snps.bam"
    conda:
        "../envs/snippy.yml"
    threads: 16
    shell:
        """
            snippy --cpus {threads} --outdir scratch/{wildcards.sample}/{wildcards.contig}/snippy/ \
            --ref {input.gbk} \
            --R1 {input.fwd} \
            --R2 {input.rev} \
            --force 

            cp scratch/{wildcards.sample}/{wildcards.contig}/snippy/snps.html {output.report}
            cp scratch/{wildcards.sample}/{wildcards.contig}/snippy/snps.bam {output.bam}
        """
