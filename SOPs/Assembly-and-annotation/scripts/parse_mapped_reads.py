import argparse
from Bio import SeqIO
import json
import gzip
import numpy as np


def parse_args():
	# Create argument parser
	parser = argparse.ArgumentParser()

	# Positional mandatory arguments
	parser.add_argument("--mapped_fwd", type=str, required=True)
	parser.add_argument("--mapped_rev", type=str, required=True)
	parser.add_argument("--mapped_single", type=str)
	parser.add_argument("--unmapped_fwd", type=str, required=True)
	parser.add_argument("--unmapped_rev", type=str, required=True)
	parser.add_argument("--unmapped_single", type=str)
	parser.add_argument("--output", type=str, required=True)
	# Parse arguments
	args = parser.parse_args()
	return args


def read_file(file):
	with gzip.open(file, 'rt') as handle:
		reads = [len(x) for x in SeqIO.parse(handle, 'fastq')]
	return len(reads), np.sum(reads)


def main(args):
	mapped_fwd_count, mapped_fwd_bp = read_file(args.mapped_fwd)
	mapped_rev_count, mapped_rev_bp = read_file(args.mapped_rev)
	unmapped_fwd_count, unmapped_fwd_bp = read_file(args.unmapped_fwd)
	unmapped_rev_count, unmapped_rev_bp = read_file(args.unmapped_rev)

	total_mapped_count = mapped_fwd_count + mapped_rev_count
	total_unmapped_count = unmapped_fwd_count + unmapped_rev_count
	total_mapped_bp = mapped_fwd_bp + mapped_rev_bp
	total_unmapped_bp = unmapped_fwd_bp + unmapped_rev_bp

	if args.mapped_single:
		mapped_single_count, mapped_single_bp = read_file(args.mapped_single)
		total_mapped_count += mapped_single_count
		total_mapped_bp += mapped_single_bp

	if args.unmapped_single:
		unmapped_single_count, unmapped_single_bp = read_file(args.unmapped_single)
		total_unmapped_count += unmapped_single_count
		total_unmapped_bp += unmapped_single_bp

	pct_mapped_count = total_mapped_count / (total_unmapped_count + total_mapped_count) * 100
	pct_mapped_bp = total_mapped_bp / (total_mapped_bp + total_unmapped_bp) * 100

	with open(args.output, 'w') as handle:
			json.dump({'total_mapped_count': f'{total_mapped_count}',
						'total_unmapped_count' : f'{total_unmapped_count}',
						'total_mapped_bp' : f'{total_mapped_bp}',
						'total_unmapped_bp' : f'{total_unmapped_bp}',
						'pct_mapped_count' : f'{pct_mapped_count:.3f}',
						'pct_mapped_bp' : f'{pct_mapped_bp:.3f}'}, handle, indent=4)


if __name__=="__main__":
	args = parse_args()
	main(args)
