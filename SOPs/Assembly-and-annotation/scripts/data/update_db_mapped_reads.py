import argparse
import json
from pymongo import MongoClient
from dotenv import load_dotenv
import os
from utils import upload_file_to_s3


def get_database():
    load_dotenv()
    client = MongoClient(os.getenv('MONGO_DB'))
    return client[os.getenv('DB_NAME')]


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--report", type=str, required=True)
    parser.add_argument("--bam_file", type=str, required=True)
    parser.add_argument("--sample", type=str, required=True)
    parser.add_argument("--version", type=str, default='1.0.0')
    parser.add_argument("--output", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args, db_client):
    phages = db_client['phages']

    with open(args.report, 'r') as handle:
        mapping = json.load(handle)

    mapping = {k: float(v) for k, v in mapping.items()}

    mapping['mapped_reads_bam'] = upload_file_to_s3(args.bam_file,
                                                   f'phages/{args.sample}/reads/host-mapped-reads.bam')

    phages.update_one({"short_name": args.sample, "version": args.version},
                          {"$set": {'sequencing.host_removal': mapping}})

    with open(args.output, 'w') as handle:
        json.dump(mapping, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    db_client = get_database()
    main(args, db_client)
