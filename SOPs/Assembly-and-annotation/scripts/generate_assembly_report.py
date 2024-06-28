import argparse
import pandas as pd
import json
import sys


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument("--qc_report", type=str, required=True)
    parser.add_argument("--mapping_report", type=str, required=True)
    parser.add_argument("--shovill_report", type=str, required=True)
    parser.add_argument("--unicycler_report", type=str, required=True)
    parser.add_argument("--onepct_report", type=str, required=True)

    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    results = {}
    with open(args.config, 'r') as handle:
        config = json.load(handle)

    metadata = pd.read_csv(config['metadata_file'], index_col=0, comment='#')
    row = metadata.loc[args.sample]

    with open(args.qc_report, 'r') as handle:
        results['sequencing'] = json.load(handle)

    results['sequencing']['reads'] = {"read_url": [f'reads/{args.sample}/orig_reads/{args.sample}-fwd.qc.fq.gz',
                                                   f'reads/{args.sample}/orig_reads/{args.sample}-fwd.qc.fq.gz'],
                                      'read_qc_url': f"{args.sample}/01_reads/orig_reads/qc-report.tgz"}

    results['sequencing']['host_removal'] = {"host_genome": row['host_name'],
                                             "genome_url": f"{args.sample}/01_reads/host-mapping/host.fa"}

    with open(args.mapping_report, 'r') as handle:
        results['sequencing']['host_removal']['mapping'] = json.load(handle)

    results['sequencing']['host_removal']['mapping'][
        "mapped_reads_bam_url"] = f"{args.sample}/01_reads/host-mapping/host-mapped-reads.bam"


    results['assemblies'] = []

    with open(args.shovill_report, 'r') as handle:
        results['assemblies'].append({
            "assembler": "shovill",
            "assembler_version": "1.1.0",
            "command": "--minlen 10000 --depth 500 --mincov 20 --keepfiles --cpus {threads} --force --noreadcorr",
            "subsampled_coverage": "500",
            "contigs_file" : f"{args.sample}/02_assembly/shovill/{args.sample}-shovill-contigs.fa.gz",
            "assembly_graph" : f"{args.sample}/02_assembly/shovill/{args.sample}-shovill-contigs.gfa.gz",
            "contigs": json.load(handle)['contigs']
        })
        with open(args.unicycler_report, 'r') as handle:
            results['assemblies'].append({
                "assembler": "unicycler",
                "assembler_version": "0.5.0",
                "command": "--min_fasta_length 1000",
                "subsampled_coverage": "500",
                "subsampled_method": {
                    "program_name": "shovill",
                    "program_version": "1.1.0",
                    "command": "--minlen 10000 --depth 500 --mincov 20 --keepfiles --cpus {threads} --force --noreadcorr"
                },
                "contigs_file": f"{args.sample}/02_assembly/unicycler/shovill-reads/contigs.fa",
                "bandage_png": f"{args.sample}/02_assembly/unicycler/shovill-reads/{args.sample}-shovill-bandage.png",
                "assembly_graph": f"{args.sample}/02_assembly/unicycler/shovill-reads/bandage.gfa.gz",
                "contigs": json.load(handle)['contigs']
            })

        with open(args.onepct_report, 'r') as handle:
            results['assemblies'].append({
                "assembler": "unicycler",
                "assembler_version": "0.5.0",
                "command": "--min_fasta_length 1000",
                "subsampled_coverage": "1pc",
                "subsampled_method": {
                    "program_name": "seqtk",
                    "program_version": "1.3",
                    "command": "seqtk sample {input.fwd} 0.01  | pigz --fast -c -p 16 > {output.sub_fwd}"
                },
                "contigs_file": f"{args.sample}/02_assembly/unicycler/1pc/contigs.fa",
                "bandage_png": f"{args.sample}/02_assembly/unicycler/1pc/{args.sample}-1pc-bandage.png",
                "assembly_graph": f"{args.sample}/02_assembly/unicycler/1pc/bandage.gfa.gz",
                "contigs": json.load(handle)['contigs']
            })

    with open(args.output, 'w') as handle:
        json.dump(results, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    main(args)
