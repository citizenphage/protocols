import argparse
import json
import uuid
from Bio import SeqIO
from pymongo import MongoClient
from dotenv import load_dotenv
import os
from utils import upload_file_to_s3
import gzip

def get_database():
    load_dotenv()
    client = MongoClient(os.getenv('MONGO_DB'))
    return client[os.getenv('DB_NAME')]


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--contigs", type=str, required=True)
    parser.add_argument("--assembler", type=str, required=True)
    parser.add_argument("--assembler_version", type=str, required=True)
    parser.add_argument("--command", type=str, required=True)
    parser.add_argument("--subsampled_coverage", type=str, required=True)
    parser.add_argument("--graph", type=str, required=True)
    parser.add_argument("--graph_img", type=str, required=True)
    parser.add_argument("--sample", type=str, required=True)
    parser.add_argument("--version", type=str, default='1.0.0')
    parser.add_argument("--output", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def parse_contigs(file):
    contigs = []
    with gzip.open(file, 'rt') as handle:
        for seq in SeqIO.parse(handle, 'fasta'):
            contigs.append({"name": seq.id, "length": len(seq.seq)})
    return contigs


def main(args, db_client):
    phages = db_client['phages']

    assembly = {
        "_id": str(uuid.uuid4()),
        "assembler": args.assembler,
        "assembler_version": args.assembler_version,
        "command": args.command,
        "subsampled_coverage": args.subsampled_coverage,
        "assembly_graph": upload_file_to_s3(args.graph,
                                            f'phages/{args.sample}/assembly/{args.assembler}/{args.subsampled_coverage}/assembly.gfa.gz'),
        "assembly_graph_img": upload_file_to_s3(args.graph_img, f'phages/{args.sample}/assembly/{args.assembler}/{args.subsampled_coverage}/bandage.png'),
        "contigs_file": upload_file_to_s3(args.graph,
                                            f'phages/{args.sample}/assembly/{args.assembler}/{args.subsampled_coverage}/contigs.fa.gz'),
        "contigs": parse_contigs(args.contigs)

    }

    phages.update_one({"short_name": args.sample,
                       "version": args.version
                       },
                      {
                        "$set": {'status': 'assembled'},
                        "$push": {'sequencing.assembly': assembly}
                      })

    with open(args.output, 'w') as handle:
        json.dump(assembly, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    db_client = get_database()
    main(args, db_client)
