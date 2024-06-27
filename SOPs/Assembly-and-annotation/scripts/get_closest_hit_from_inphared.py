from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import argparse
import json

def parse_args():
	# Create argument parser
	parser = argparse.ArgumentParser()

	# Positional mandatory arguments
	parser.add_argument("--inphared_file", type=str, required=True)
	parser.add_argument("--outfile", type=str, required=True)
	parser.add_argument("--report", type=str, required=True)
	parser.add_argument("--email", type=str, default='admin@citizenphage.com')
	parser.add_argument("--inphared_version", type=str, default='1Aug2023')
	# Parse arguments
	args = parser.parse_args()
	return args

def main(args):
	df = pd.read_csv(args.inphared_file, sep='\t')
	accession = df.iloc[0]['Accession']
	print(f'Looking for accession number {accession}')
	seqs = []
	report_obj = {'estimated_distance': df.iloc[0]['mash_distance'], 
					"estimated_distance_method": "mash",
					"database_name": "inphared",
					"database_version": args.inphared_version }

	if accession:
		Entrez.email = args.email
		with Entrez.efetch(db='nuccore', id = accession, rettype="fasta", retmode="text") as handle:
			for record in SeqIO.parse(handle, 'fasta'):
				seqs.append(record)
			report_obj['name'] = seqs[0].description
			report_obj['accession_number'] = accession

		with Entrez.efetch(db='nuccore', id = accession, rettype="gb", retmode="text") as handle:
			x = SeqIO.read(handle, 'genbank')
			taxonomy = x.annotations['taxonomy']
			if taxonomy:
				report_obj['taxonomy'] = '; '.join(taxonomy)


	with open(args.outfile, 'w') as handle:
		SeqIO.write(seqs, handle, 'fasta')

	with open(args.report, 'w') as handle:
		json.dump(report_obj, handle, indent=4)


if __name__=="__main__":
	args = parse_args()
	main(args)
