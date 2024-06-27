import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import json
import gzip
import re
def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--unicycler_log", type=str, required=True)
    parser.add_argument("--shovill_contigs", type=str, required=True)
    parser.add_argument("--checkv_log", type=str, required=True)
    parser.add_argument("--inphared_log", type=str, required=True)
    parser.add_argument("--pharokka_log", type=str, required=True)
    parser.add_argument("--coverage_log", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def parse_unicycler(logfile):
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    # Initialize empty lists to store data for the first DataFrame
    data_1 = []

    # Initialize empty lists to store data for the second DataFrame

    results = {}
    in_first_dataframe = False
    in_second_dataframe = False
    with open(args.unicycler_log, 'r') as file:
        for line in file:
            line = line.strip()
            line = ansi_escape.sub('', line)
            if line.startswith("Component"):
                parts = [x.replace('\x1b[0m', '') for x in re.split(r'\s+', next(file)) if x and x not in ['\x1b[32m', '\x1b[31m']]
                parts = [x.replace(',', '') for x in parts]
                data_1.append(parts)
                continue
            if line.startswith("Segment"):
                parts = [x for x in re.split(r'\s{2,}', next(file)) if x]
                
                results['terminase_found'] = not parts[3] =='\x1b[31mnone found\x1b[0m'
                continue

    df1 = pd.DataFrame(data_1, columns=["Component", "Segments", "Links", "Length", "N50", "Longest_segment", "Status"])
    df1['Component'] = pd.to_numeric(df1['Component'])
    df1['Segments'] = pd.to_numeric(df1['Segments'])
    df1['Links'] = pd.to_numeric(df1['Links'])
    df1['N50'] = pd.to_numeric(df1['N50'])
    df1['Longest_segment'] = pd.to_numeric(df1['Longest_segment'])
    df1['Length'] = pd.to_numeric(df1['Length'])

    
    results['assembly'] = df1.to_dict('records')
    return results

def parse_shovill(infile):
    results = {}
    contigs = {}
    with gzip.open(infile, 'rt') as handle:
        for seq in SeqIO.parse(handle, 'fasta'):
            contigs[seq.id] = len(seq)
    
    results['num_contigs'] = len(contigs)

    results['contig_lengths'] = contigs
    return results

def parse_pharokka(infile):
    results = {}

    with open(args.pharokka_log, 'r') as file:
        virulence_factors = np.NaN
        amr_genes = np.NaN
        for line in file:
            line = line.strip()
            if line.endswith('VFDB virulence factors identified.'):
                virulence_factors = int(line.split()[0])
                continue
            if line.endswith('CARD AMR genes identified.'):
                amr_genes = int(line.split()[0])
                continue

    results['virulence_factors'] = virulence_factors
    results['amr_genes'] = amr_genes
    return results

def parse_checkv(logfile):
    results = {}
    df = pd.read_csv(logfile, sep='\t')
    for index, row in df.iterrows():
        results[row['contig_id']] = {'contig_length': row['contig_length'], 'quality': row['checkv_quality']}
    return results

def parse_inphared(logfile):
    results = {}
    df = pd.read_csv(logfile, sep='\t')
    for index, row in df.iterrows():
        results[row['contig']] = {'best_hit_accession': row['Accession'],
                                    'best_hit_name': row['Description'],
                                    'best_hit_length': row['Genome_Length_(bp)'],
                                    'taxonomy': {
                                    'genus': row['Genus'],
                                    'sub-family': row['Sub-family'],
                                    'family': row['Family'],
                                    'order': row['Order'],
                                    'class': row['Class'],
                                    'phylum': row['Phylum'],
                                    'kingdom': row['Kingdom']
                                    },
                                    'mash_similarity': f'{(1 - row["mash_distance"]) * 100:.2f} %'}
    return results

def parse_coverage(logfile):
    with open(logfile, 'r') as handle:
        results = json.load(handle)
    return float(results['coverage'])



def main(args):
    results = {}
    results['unicycler'] = parse_unicycler(args.unicycler_log)
    results['shovill'] = parse_shovill(args.shovill_contigs)
    results['pharokka'] = parse_pharokka(args.pharokka_log)
    results['checkv'] = parse_checkv(args.checkv_log)
    results['nearest_inphared_hit'] = parse_inphared(args.inphared_log)
    results['coverage'] = parse_coverage(args.coverage_log)

    with open(args.output, 'w') as handle:
        json_str = json.dump(results, handle, indent=4)
    

if __name__=="__main__":
    args = parse_args()
    main(args)




