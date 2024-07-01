import argparse
import pandas as pd
import json


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--vfdb", type=str, required=True)
    parser.add_argument("--card", type=str, required=True)
    parser.add_argument("--gbk", type=str, required=True) # might need this in the future
    parser.add_argument("--output", type=str, required=True)

    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    results = {"predicted_genes_of_concern": [],
               "prediction_methods": []}

    vfdb = pd.read_csv(args.vfdb, sep='\t')
    for index, row in vfdb.iterrows():
        results['predicted_genes_of_concern'].append({
            "gene": row['gene'],
            "hit": row['vfdb_hit'],
            "method": 'vfdb'
        })

    results['prediction_methods'].append({
        "method": "vfdb",
        "command": "pharokka.py -i {input.contigs} -o scratch/{wildcards.sample}/pharokka-unicycler -d {input.db} -t {threads} -l {wildcards.sample} -g prodigal --force --dnaapler --prefix {wildcards.sample}"
    })

    card = pd.read_csv(args.card, sep='\t')
    for index, row in card.iterrows():
        results['predicted_genes_of_concern'].append({
            "gene": row['gene'],
            "hit": row['card_hit'],
            "method": 'card'
        })

    results['prediction_methods'].append({
        "method": "card",
        "command": "pharokka.py -i {input.contigs} -o scratch/{wildcards.sample}/pharokka-unicycler -d {input.db} -t {threads} -l {wildcards.sample} -g prodigal --force --dnaapler --prefix {wildcards.sample}"
    })

    with open(args.output, 'w') as handle:
        json.dump(results, handle, indent=4)

if __name__ == "__main__":
    args = parse_args()
    main(args)
