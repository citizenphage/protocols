import pandas as pd
samples = pd.read_csv('resources/metadata.csv', index_col=0, comment='#')

def get_fwd_reads(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.short_fwd}'

def get_rev_reads(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.short_rev}'

def get_host(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.host_genome}'

def get_phage_gbk(wildcards):
    row = samples.loc[wildcards.sample]
    return f'{row.phage_gbk}'

rule all:
    input:
        expand("output/{sample}/fastqc-fwd.html", sample=samples.index),
        expand("output/{sample}/fastqc-rev.html", sample=samples.index),
        expand("output/{sample}/qc-report.json", sample=samples.index),
        expand("output/{sample}/mapping/{sample}-mapped-reads.bam", sample=samples.index),
        expand("output/{sample}/assembly/unmapped-contigs.fa.gz", sample=samples.index),
        expand("output/{sample}/assembly/non-phage-contigs.fa.gz",sample=samples.index),
        expand("output/{sample}/assembly/{sample}-non-phage-host-mapping.bam", sample=samples.index),
        expand("output/{sample}/mapping/{sample}-coverage.json", sample=samples.index),
        expand("output/{sample}/assembly/checkv-summary-non-phage-contigs.tsv", sample=samples.index)




rule fastqc_reads:
    input:
        fwd="scratch/{sample}/fwd.fq.gz",
        rev="scratch/{sample}/rev.fq.gz"
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
        fwd="scratch/{sample}/fwd.fq.gz",
        rev="scratch/{sample}/rev.fq.gz"
    output:
        fwd="scratch/{sample}/fwd.qc.fq.gz",
        rev="scratch/{sample}/rev.qc.fq.gz",
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
            --length_required 30 \
            --correction \
            --json {output.report} \
            --thread {threads} 2>&1 | tee {log}


        """

rule create_fasta_for_mapping:
    input:
        phage="scratch/{sample}/phage.gbk",
        host="scratch/{sample}/host.fa"
    output:
        combined="output/{sample}/{sample}-combined-genomes.fa",
        phage="scratch/{sample}/phage.fa"
    shell:
        """
        mkdir -p scratch/{wildcards.sample}
        python scripts/extract_fasta_from_genbank.py \
        --gbk {input.phage} \
        --outfile {output.phage}

        cat {output.phage} {input.host} > {output.combined}

        """


rule map_reads:
    input:
        fwd=rules.qc_reads.output.fwd,
        rev=rules.qc_reads.output.rev,
        phage=rules.create_fasta_for_mapping.output.phage,
        target=rules.create_fasta_for_mapping.output.combined
    output:
        mapped="output/{sample}/mapping/{sample}-mapped-reads.bam",
        combined="output/{sample}/mapping/{sample}-all-reads.bam",
        unmapped=temp("scratch/{sample}/unmapped-reads.bam"),
        mapped_fwd=temp("scratch/{sample}/mapped-reads-fwd.fq.gz"),
        mapped_rev=temp("scratch/{sample}/mapped-reads-rev.fq.gz"),
        mapped_single=temp("scratch/{sample}/mapped-reads-single.fq.gz"),
        unmapped_fwd=temp("scratch/{sample}/unmapped-reads-fwd.fq.gz"),
        unmapped_rev=temp("scratch/{sample}/unmapped-reads-rev.fq.gz"),
        unmapped_single=temp("scratch/{sample}/unmapped-reads-single.fq.gz"),
        non_phage_fwd=temp("scratch/{sample}/non-phage-reads-fwd.fq.gz"),
        non_phage_rev=temp("scratch/{sample}/non-phage-reads-rev.fq.gz"),
        non_phage_single=temp("scratch/{sample}/non-phage-reads-single.fq.gz")

    conda:
        "../envs/map-reads.yml"
    threads: 16
    log: "logs/{sample}/map-reads.log"
    shell:
        """
        minimap2 \
        -a -x sr -t {threads} \
        -o scratch/{wildcards.sample}/mapping.sam \
        {input.target} \
        {input.fwd} \
        {input.rev} 2>&1 | tee {log}

        samtools view -bS \
        -@ {threads} scratch/{wildcards.sample}/mapping.sam | samtools sort -@ {threads} -o {output.combined}
        samtools index {output.combined}

        samtools view -bS -F 4 \
        -@ {threads} scratch/{wildcards.sample}/mapping.sam | samtools sort -@ {threads} -o {output.mapped}

        samtools view -bS -f 4 \
        -@ {threads} scratch/{wildcards.sample}/mapping.sam | samtools sort -@ {threads} -o {output.unmapped}

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


        minimap2 \
        -a -x sr -t {threads} \
        -o scratch/{wildcards.sample}/phage-mapping.sam \
        {input.phage} \
        {input.fwd} \
        {input.rev} 2>&1 | tee {log}

        samtools view -bS -f 4 \
        -@ {threads} scratch/{wildcards.sample}/phage-mapping.sam | samtools sort -@ {threads} -o scratch/{wildcards.sample}/phage-mapping.bam

        samtools fastq \
        -0 /dev/null \
        -1 {output.non_phage_fwd} \
        -2 {output.non_phage_rev} \
        -s {output.non_phage_single} \
        -@ {threads} \
        -n scratch/{wildcards.sample}/phage-mapping.bam

        rm scratch/{wildcards.sample}/phage-mapping.sam
        rm scratch/{wildcards.sample}/phage-mapping.bam


        """

rule assemble_unmapped_reads:
    input:
        fwd=rules.map_reads.output.unmapped_fwd,
        rev=rules.map_reads.output.unmapped_rev
    output:
        contigs="output/{sample}/assembly/unmapped-contigs.fa.gz",
        graph="output/{sample}/assembly/{sample}-unmapped-contigs.gfa.gz"
    conda:
        "../envs/assembly.yml"
    log: "logs/{sample}/unmapped-assembly.log"
    params:
        min_length=3000
    threads: 16
    shell:
        """
            spades.py \
            -1 {input.fwd} \
            -2 {input.rev} \
            -t {threads} \
            --only-assembler \
            --meta \
            -o scratch/{wildcards.sample}/unmapped-spades 2>&1 | tee {log}

            pigz -c scratch/{wildcards.sample}/unmapped-spades/contigs.fasta > {output.contigs}
            pigz -c scratch/{wildcards.sample}/unmapped-spades/assembly_graph.fastg > {output.graph}

            rm -rf scratch/{wildcards.sample}/unmapped-spades
        """


rule assemble_non_phage_reads:
    input:
        fwd=rules.map_reads.output.non_phage_fwd,
        rev=rules.map_reads.output.non_phage_rev,
        single=rules.map_reads.output.non_phage_single
    output:
        contigs="output/{sample}/assembly/non-phage-contigs.fa.gz",
        graph="output/{sample}/assembly/{sample}-non-phage-contigs.gfa.gz"
    conda:
        "../envs/assembly.yml"
    log: "logs/{sample}/non-phage-assembly.log"
    threads: 16
    shell:
        """
            spades.py \
            -1 {input.fwd} \
            -2 {input.rev} \
            -s {input.single} \
            -t {threads} \
            --only-assembler \
            --meta \
            -o scratch/{wildcards.sample}/non-phage-spades 2>&1 | tee {log}



            pigz -c scratch/{wildcards.sample}/non-phage-spades/contigs.fasta > {output.contigs}
            pigz -c scratch/{wildcards.sample}/non-phage-spades/assembly_graph.fastg > {output.graph}

            rm -rf scratch/{wildcards.sample}/non-phage-spades
        """

rule check_assembled_contigs_for_prophage:
    input:
        rules.assemble_unmapped_reads.output.contigs,
    output:
        quality="output/{sample}/assembly/checkv-summary-non-phage-contigs.tsv",
        contigs="output/{sample}/assembly/checkv-summary-non-phage-viral-contigs.fa"
    conda:
        "../envs/checkv.yml"
    log: "logs/{sample}/non-phage-assembly-checkv.log"
    params:
        checkv_db="/home/artic/bens_toys/dbs/checkv-db-v1.1"
    threads: 16
    shell:
        """
            checkv end_to_end \
            {input} \
            scratch/{wildcards.sample}/checkv \
            -t {threads} \
            -d {params.checkv_db}

            cp scratch/{wildcards.sample}/checkv/quality_summary.tsv {output.quality}
            cp scratch/{wildcards.sample}/checkv/viruses.fna {output.contigs}

            rm -rf scratch/{wildcards.sample}/checkv
            
        """

rule map_non_phage_assembly:
    input:
        contigs=rules.assemble_non_phage_reads.output.contigs,
        host="scratch/{sample}/host.fa"
    output:
        "output/{sample}/assembly/{sample}-non-phage-host-mapping.bam"
    conda:
        "../envs/map-reads.yml"
    threads: 16
    log: "logs/{sample}/map-non-host-contigs.log"
    shell:
        """
            minimap2 \
            -L -a -x asm5 -t {threads} \
            -o scratch/{wildcards.sample}/non_phage_assembly.sam \
            {input.host} \
            {input.contigs} 2>&1 | tee {log}

            samtools view -bS scratch/{wildcards.sample}/non_phage_assembly.sam | samtools sort -o {output}

            rm scratch/{wildcards.sample}/non_phage_assembly.sam

        """


rule evaluate_coverage:
    input:
        phage=rules.create_fasta_for_mapping.output.phage,
        host="scratch/{sample}/host.fa",
        bam=rules.map_reads.output.combined
    output:
        "output/{sample}/mapping/{sample}-coverage.json"
    shell:
        """
            python scripts/evaluate_bam_for_coverage.py \
            --bam {input.bam} \
            --phage_fasta {input.phage} \
            --host_fasta {input.host} \
            --output {output}
        """
