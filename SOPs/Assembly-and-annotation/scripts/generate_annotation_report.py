import argparse
from Bio import SeqIO
import pandas as pd
import json
import os



def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--fasta", type=str, required=True)
    parser.add_argument("--checkv", type=str, required=True)
    parser.add_argument("--coverage", type=str, required=True)
    parser.add_argument("--genes_of_concern", type=str, required=True)
    parser.add_argument("--sample_name", type=str, required=True)
    parser.add_argument("--variants", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--config", type=str, default='resources/config.json')

    # Parse arguments
    args = parser.parse_args()
    return args


def check_file(file_path, local_store, backblaze_store):
    if os.path.isfile(f'{local_store}/{file_path}'):
        return f'{backblaze_store}/{file_path}'
    print(f'Could not find file: [{file_path}]')
    return 'FILE_NOT_FOUND'


def main(args):
    with open(args.config, 'r') as handle:
        config = json.load(handle)
        local_store = config['local_data_store']
        backblaze_store = config['data_store']


    results = {}

    with open(args.fasta, 'r') as handle:
        seq = [x for x in SeqIO.parse(handle, 'fasta')][0]
        results['selected_contig_name'] = seq.id
        results['selected_contig_length'] = len(seq)

    with open(args.checkv, 'r') as handle:
        results['checkv'] = json.load(handle)

    results['gene_calling'] = {}

    results['gene_calling']['program_used'] = "pharokka + phold"
    results['gene_calling']['program_version'] = "1.7.1"
    results['gene_calling']['command'] = '-g prodigal --force --dnaapler'
    results['gene_calling'][
        "genbank_url"] = check_file(f"{args.sample_name}/05_annotation/{results['selected_contig_name']}/phold.gbk", local_store, backblaze_store)
    results['gene_calling'][
        "plot_url"] = check_file(f"{args.sample_name}/05_annotation/{results['selected_contig_name']}/phold.png", local_store, backblaze_store)

    results['lifestyle'] = {}
    results['lifestyle']['predicted_lifestyle'] = "TO_COMPLETE"
    results['lifestyle']['prediction_methods'] = []
    results['lifestyle']['prediction_methods'].append({
        "program_name": "PhageLeads",
        "program_version": "https://doi.org/10.3390/v14020342",
        "command": "https://phageleads.ku.dk/",
        "result": "TO_COMPLETE"
    })

    with open(args.genes_of_concern, 'r') as handle:
        results['genes_of_concern'] = json.load(handle)

    results['genome_coverage'] = {}
    with open(args.coverage, 'r') as handle:
        results['genome_coverage']["assembly_average_coverage"] = json.load(handle)['coverage']
    results['genome_coverage'][
        "assembly_bam_file"] = check_file(f"{args.sample_name}/06_coverage/{results['selected_contig_name']}/assembly-mapped-reads.bam", local_store,backblaze_store)
    results['genome_coverage'][
        "assembly_coverage_plot_file"] = check_file(f"{args.sample_name}/06_coverage/{results['selected_contig_name']}/coverage-{results['selected_contig_name']}.png", local_store, backblaze_store)

    results['variants'] = {
        "program_name": "snippy",
        "program_version": "4.6.0",
        "command": "snippy --cpus 16"
    }

    results['variants']["full_read_bam_file"] = check_file(f"{args.sample_name}/07_variants/{results['selected_contig_name']}/snps.bam", local_store, backblaze_store)

    snps = pd.read_csv(args.variants).to_dict('records')

    results['variants']['variants'] = snps

    with open(args.output, 'w') as handle:
        json.dump(results, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    main(args)
