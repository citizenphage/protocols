import pysam
import matplotlib.pyplot as plt
import argparse
from Bio import SeqIO
import seaborn as sns
import numpy as np
import json

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--genome", type=str, required=True)
    parser.add_argument("--window_size", type=int, default=1000)
    parser.add_argument("--json", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args

def calculate_sliding_window_coverage(coverage, window_size):
    coverage_mean = []
    for i in range(len(coverage) - window_size + 1):
        window_coverage = coverage[i:i + window_size]
        window_mean = np.mean(window_coverage)
        coverage_mean.append(window_mean)
    return coverage_mean

def plot_genome_coverage(genome, bam_file, prefix, window_size):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    contigs = [record.id for record in SeqIO.parse(genome, 'fasta')]
    
    for c in contigs:
    # Retrieve the coverage information
        coverage = bam.count_coverage(contig = c)

        # Calculate the total coverage across all positions
        sliding_window_coverage = calculate_sliding_window_coverage(coverage[0], window_size)

        total_coverage = np.mean(sliding_window_coverage)

        # Generate the x-axis positions for plotting
        positions = list(range(1, len(sliding_window_coverage) + 1))

        ## Set seaborn style
        sns.set(style="whitegrid")

        # Plot the coverage
        plt.plot(positions, sliding_window_coverage, linewidth=1.5, color='#1f77b4')

        # Format the plot
        plt.xlabel("Locus", fontsize=12)
        plt.ylabel("Coverage", fontsize=12)
        plt.title(f"Genome Coverage (total = {total_coverage:.1f}x)", fontsize=14)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout()

        # Save the plot as a PNG file with 300 dpi
        plt.savefig(f'{prefix}-{c}.png', dpi=300)

        # Display a confirmation message
        print(f"Genome coverage plot saved as {prefix}-{c}.png")

        with open(args.json, 'w') as handle:
            json.dump({'coverage': f'{total_coverage:.1f}'}, handle)

def main(args):
    plot_genome_coverage(args.genome, args.bam, args.prefix, args.window_size)

if __name__=="__main__":
    args = parse_args()
    main(args)
