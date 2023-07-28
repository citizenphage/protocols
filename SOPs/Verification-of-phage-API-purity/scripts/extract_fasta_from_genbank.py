from Bio import SeqIO
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--gbk", type=str, required=True)
    parser.add_argument("--outfile", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args

def extract_fasta_from_genbank(genbank_file):
    fasta_records = []
    with open(genbank_file) as handle:
	    for record in SeqIO.parse(handle, "genbank"):
	        # Check if the record contains a sequence
	        if record.seq:
	            fasta_records.append(record)
    
    return fasta_records

def main(args):
	records = extract_fasta_from_genbank(args.gbk)
	if records:
		count = SeqIO.write(records, args.outfile, 'fasta')
		print(f'Printed {count} records')
	else:
		print('No records were found')


if __name__=="__main__":
    args = parse_args()
    main(args)