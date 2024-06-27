import argparse
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import json


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--viruses", type=str, required=True)
    parser.add_argument("--quality", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument('--contig_store', type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    results = {}
    results['version'] = "1.0.1"
    results['command'] = "checkv end_to_end {input.assembly} scratch/{wildcards.sample}/shovill-checkv -t {threads} -d {input.db}/checkv-db-v* --remove_tmp"
    results['contigs'] = []
    checkv_df = pd.read_csv(args.quality, sep='\t')
    genomes = SeqIO.index(args.viruses, 'fasta')
    Path(args.contig_store).mkdir(parents=True, exist_ok=True)

    for index, row in checkv_df.iterrows():
        if row['checkv_quality'] in ['Complete', 'High-quality', 'Medium-quality']:
            results['contigs'].append({"name": row['contig_id'], "quality": row['checkv_quality'],"estimated_completeness": f"{row['completeness']:.2f}"})
            with open(f"{args.contig_store}/{row['contig_id']}.fa", 'w') as handle:
                SeqIO.write([genomes[row['contig_id']]], handle, 'fasta')

    with open(args.output, 'w') as handle:
        json.dump(results, handle, indent=4)



if __name__ == "__main__":
    args = parse_args()
    main(args)
