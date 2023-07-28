import pysam
import argparse
from Bio import SeqIO
import numpy as np
import json


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--phage_fasta", type=str, required=True)
    parser.add_argument("--host_fasta", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    bam_file = pysam.AlignmentFile(args.bam, "rb")

    phage_contigs = [x.id for x in SeqIO.parse(args.phage_fasta, 'fasta')]
    host_contigs = [x.id for x in SeqIO.parse(args.host_fasta, 'fasta')]

    read_counts = {}

    phage_read_count = 0
    unmapped_read_count = 0


    bacterial_read_count = 0
    
    phage_total_coverage = 0
    
    # Iterate through each read in the BAM file
    for read in bam_file.fetch():
        # Check if the read is mapped (not unmapped)
        if read.is_unmapped:
            unmapped_read_count += 1
        else:
            if read.reference_name in phage_contigs:
                phage_total_coverage += read.query_alignment_length
        try:
            read_counts[read.reference_name] += 1
        except KeyError:
            read_counts[read.reference_name] = 1
    
    for k, v in read_counts.items():
        if k in phage_contigs:
            phage_read_count += v
        else:
            bacterial_read_count += v
    
    total_reads = phage_read_count + unmapped_read_count + bacterial_read_count
    
    phage_genome_length = np.sum([bam_file.get_reference_length(x) for x in phage_contigs])
    
    results = {
        'total_reads': total_reads,
        'bacterial_reads': bacterial_read_count,
        'bacterial_pct': np.around(bacterial_read_count / total_reads * 100, decimals=2),
        'phage_reads': phage_read_count,
        'phage_pct': np.around(phage_read_count / total_reads * 100, decimals=2),
        'unmapped_reads': unmapped_read_count,
        'unmapped_pct': np.around(unmapped_read_count / total_reads * 100, decimals=2),
        'phage_coverage': np.around(phage_total_coverage / phage_genome_length, decimals=2),
        'read_counts': read_counts
    
    }
    json_object = json.dumps(results, indent=4)

    
    with open(args.output, "w") as outfile:
        outfile.write(json_object)

if __name__ == "__main__":
    args = parse_args()
    main(args)
