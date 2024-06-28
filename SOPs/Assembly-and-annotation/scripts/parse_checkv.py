import argparse
import pandas as pd
import json
import shutil


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--contig", type=str, required=True)
    parser.add_argument("--quality", type=str, required=True)
    parser.add_argument("--report", type=str, required=True)
    parser.add_argument("--output_contig", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    results = {}
    results['version'] = "1.0.1"
    results['command'] = "checkv end_to_end {input.assembly} scratch/{wildcards.sample}/shovill-checkv -t {threads} -d {input.db}/checkv-db-v* --remove_tmp"
    results['contigs'] = []
    checkv_df = pd.read_csv(args.quality, sep='\t')

    for index, row in checkv_df.iterrows():
        results['contigs'].append({"name": row['contig_id'], "quality": row['checkv_quality'],"estimated_completeness": f"{row['completeness']:.2f}"})
        if row['checkv_quality'] in ['Complete', 'High-quality', 'Medium-quality']:
            shutil.copy(args.contig, args.output_contig)


    with open(args.report, 'w') as handle:
        json.dump(results, handle, indent=4)



if __name__ == "__main__":
    args = parse_args()
    main(args)
