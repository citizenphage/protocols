import argparse
import numpy as np
import pandas as pd
import json
import os

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--reports", type=str, nargs='+')
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    results = []
    for f in args.reports:
        with open(f, 'r') as handle:
            try:
                phage_id = os.path.basename(f).split('-report.json')[0]
                d = json.load(handle)
                unicycler_subsample_num_segments = d['unicycler']['assembly'][0]['Segments']
                unicycler_subsample_status = d['unicycler']['assembly'][0]['Status']
                unicycler_subsample_length = d['unicycler']['assembly'][0]['Length']
                unicycler_found_terminase = d['unicycler']['terminase_found']
                
                shovill_num_contigs = d['shovill']['num_contigs']
                shovill_length = d['shovill']['contig_lengths'][phage_id]
                checkv_status = d['checkv'][phage_id]['quality']

                closest_hit_accession = d['nearest_inphared_hit'][phage_id]['best_hit_accession']
                closest_hit_name = d['nearest_inphared_hit'][phage_id]['best_hit_name']
                closest_hit_length = d['nearest_inphared_hit'][phage_id]['best_hit_length']
                closest_hit_similarity = d['nearest_inphared_hit'][phage_id]['mash_similarity']

                taxonomy=[d['nearest_inphared_hit'][phage_id]['taxonomy']['genus'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['sub-family'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['family'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['order'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['class'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['phylum'],
                            d['nearest_inphared_hit'][phage_id]['taxonomy']['kingdom']]
                
                taxonomy_str = '; '.join(taxonomy[::-1])

                proposed_length = 0
                if unicycler_subsample_status == 'complete':
                    proposed_length = unicycler_subsample_length
                elif shovill_num_contigs == 1:
                    proposed_length = shovill_length



                difference_in_length= float(f'{proposed_length / closest_hit_length * 100:.2f}')
                
                results.append([phage_id,
                                unicycler_subsample_num_segments,
                                unicycler_subsample_status,
                                unicycler_subsample_length,
                                unicycler_found_terminase,
                                shovill_num_contigs,
                                shovill_length,
                                checkv_status,
                                proposed_length,
                                closest_hit_accession,
                                closest_hit_name,
                                closest_hit_length,
                                float(closest_hit_similarity.split()[0]),
                                difference_in_length,
                                taxonomy_str])
            except KeyError:
                print(f'An error occured when parsing file {f}')

    df = pd.DataFrame(results, columns=['phage_id',
                                        'unicycler_num_segments',
                                        'unicycler_status',
                                        'unicycler_length',
                                        'unicycler_found_terminase',
                                        'shovill_num_contigs',
                                        'shovill_length',
                                        'checkv_status',
                                        'proposed_length',
                                        'closest_hit_accession',
                                        'closest_hit_name',
                                        'closest_hit_length',
                                        'closest_hit_similarity',
                                        'difference_in_length',
                                        'closest_hit_taxonomy'])
    print(df)

    df.to_csv(args.output, sep='\t', index=False)


    

if __name__=="__main__":
    args = parse_args()
    main(args)