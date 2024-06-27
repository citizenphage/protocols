from Bio import SeqIO
import argparse
import json


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--outfile", type=str, required=True)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--report", type=str, required=True)
    parser.add_argument("--assembly_method", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args

def main(args):
    
    reads = []
    old_names = []
    new_names = []
    lengths = []
    with open(args.input) as handle:
        reads = [x for x in SeqIO.parse(handle, 'fasta')]
        old_names = [x.id for x in reads]

    if len(reads) == 1:
        reads[0].id = args.prefix
        reads[0].description = f'{reads[0].description} {args.assembly_method}'
    else:
        counter = 1
        for r in reads:
            r.id = f'{args.prefix}_{counter:02d}'
            r.description = f'{r.description} {args.assembly_method}'
            counter +=1

    new_names = [x.id for x in reads]
    lengths = [len(x) for x in reads]

    with open(args.outfile, 'w') as handle:
        SeqIO.write(reads, handle, 'fasta')

    report = {"contigs": []}
    for i in range(len(old_names)):
        report['contigs'].append({'old_name': old_names[i], 'new_name': new_names[i], 'length': f'{lengths[i]}'})

    with open(args.report, 'w') as handle:
        json.dump(report, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    main(args)