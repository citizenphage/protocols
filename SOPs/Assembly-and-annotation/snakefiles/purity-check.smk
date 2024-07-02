

rule map_for_purity:
    input:
        host="output/{sample}/01_reads/host-mapping/host.fa",
        viral_contigs="output/{sample}/03_selected_contigs/",
        fwd="output/reads/{sample}/orig_reads/{sample}-fwd.qc.fq.gz",
        rev="output/reads/{sample}/orig_reads/{sample}-rev.qc.fq.gz"
    output:
        unmapped_fwd="scratch/{sample}/03_selected_contigs/unmapped.fwd.fq.gz",
        unmapped_rev="scratch/{sample}/03_selected_contigs/unmapped.rev.fq.gz",
        unmapped_single="scratch/{sample}/03_selected_contigs/unmapped.single.fq.gz",
        mapped_fwd="scratch/{sample}/03_selected_contigs/mapped.fwd.fq.gz",
        mapped_rev="scratch/{sample}/03_selected_contigs/mapped.rev.fq.gz",
        mapped_single="scratch/{sample}/03_selected_contigs/mapped.single.fq.gz"
    threads: 16
    conda:
        "../envs/map-reads.yml"
    shell:
        """
            mkdir -p scratch/{wildcards.sample}/03_selected_contigs
            cat {input.host} {input.viral_contigs}/*.fa > scratch/{wildcards.sample}/03_selected_contigs/genomes.fa
            
            minimap2 \
            -a -x sr -t {threads} \
            -o scratch/{wildcards.sample}/03_selected_contigs/mapping.sam \
            scratch/{wildcards.sample}/03_selected_contigs/genomes.fa \
            {input.fwd} \
            {input.rev}
            
            samtools view -bS -F 4 \
            -@ {threads} scratch/{wildcards.sample}/03_selected_contigs/mapping.sam | 
            samtools sort -@ {threads} -o scratch/{wildcards.sample}/03_selected_contigs/mapped.bam
    
            samtools view -bS -f 4 \
            -@ {threads} scratch/{wildcards.sample}/03_selected_contigs/mapping.sam | 
            samtools sort -@ {threads} -o scratch/{wildcards.sample}/03_selected_contigs/unmapped.bam
            
            rm scratch/{wildcards.sample}/03_selected_contigs/mapping.sam
            
            samtools fastq \
            -0 /dev/null \
            -1 {output.unmapped_fwd} \
            -2 {output.unmapped_rev} \
            -s {output.unmapped_single} \
            -@ {threads} \
            -n scratch/{wildcards.sample}/03_selected_contigs/unmapped.bam
            
            samtools fastq \
            -0 /dev/null \
            -1 {output.mapped_fwd} \
            -2 {output.mapped_rev} \
            -s {output.mapped_single} \
            -@ {threads} \
            -n scratch/{wildcards.sample}/03_selected_contigs/mapped.bam
            
            rm scratch/{wildcards.sample}/03_selected_contigs/mapped.bam
            rm scratch/{wildcards.sample}/03_selected_contigs/unmapped.bam
            
            
            
        """

rule parse_mapping:
    input:
        unmapped_fwd = rules.map_for_purity.output.unmapped_fwd,
        unmapped_rev = rules.map_for_purity.output.unmapped_rev,
        unmapped_single = rules.map_for_purity.output.unmapped_single,
        mapped_fwd = rules.map_for_purity.output.mapped_fwd,
        mapped_rev = rules.map_for_purity.output.mapped_rev,
        mapped_single = rules.map_for_purity.output.mapped_single
    output:
        "output/{sample}/03_selected_contigs/purity-mapping.json"
    shell:
        """
            python scripts/parse_mapped_reads.py \
            --mapped_fwd {input.mapped_fwd} \
            --mapped_rev {input.mapped_rev} \
             --mapped_single {input.mapped_single} \
            --unmapped_fwd {input.unmapped_fwd} \
            --unmapped_rev {input.unmapped_rev} \
            --unmapped_single {input.unmapped_single} \
            --output {output}
        """

rule assemble_unmapped:
    input:
        fwd=rules.map_for_purity.output.unmapped_fwd,
        rev=rules.map_for_purity.output.unmapped_rev,
        single=rules.map_for_purity.output.unmapped_single
    output:
        contigs="output/{sample}/03_selected_contigs/unmapped-contigs.fa.gz",
        report="scratch/{sample}/03_selected_contigs/unmapped-contigs-report.json"
    threads: 16
    params:
        min_length=1000
    conda:
        "../envs/assembly.yml"
    shell:
        """
            spades.py -1 {input.fwd} -2 {input.rev} -s {input.single} --only-assembler --threads {threads} -o scratch/{wildcards.sample}/03_selected_contigs/spades
            
            python scripts/rename_contigs.py \
            --input scratch/{wildcards.sample}/03_selected_contigs/spades/contigs.fasta \
            --outfile scratch/{wildcards.sample}/03_selected_contigs/spades/renamed.fa \
            --min_length {params.min_length} \
            --report {output.report} \
            --prefix unassembled \
            --assembly_method "spades meta"
            
            pigz -c scratch/{wildcards.sample}/03_selected_contigs/spades/renamed.fa > {output.contigs}
            
        """

