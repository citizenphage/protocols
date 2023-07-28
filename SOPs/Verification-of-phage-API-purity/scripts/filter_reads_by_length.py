import argparse
from Bio import SeqIO

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--length", type=int, default=3000)
    # Parse arguments
    args = parser.parse_args()
    return args

def main(args):

	contigs = [x for x in SeqIO.parse(args.input, 'fasta') if len(x.seq) >= args.length]
	with open(args.output, 'w') as handle:
		SeqIO.write(contigs, handle, 'fasta')

if __name__ == "__main__":
    args = parse_args()
    main(args)
