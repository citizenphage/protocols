import argparse
import json


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--report", type=str, required=True)
    parser.add_argument("--date", type=str, required=True)
    parser.add_argument("--centre", type=str, required=True)
    parser.add_argument("--type", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    # Parse arguments
    args = parser.parse_args()
    return args


def main(args):
    output = {}
    with open(args.report, 'r') as handle:
        report_data = json.load(handle)

    output["date_sequenced"] = args.date
    output["sequencing_center"] = args.centre
    output["sequencing_type"] = args.type
    output["number_of_reads"] = report_data['summary']['after_filtering']['total_reads']
    output["total_bp"] = report_data['summary']['after_filtering']['total_bases']

    with open(args.output, 'w') as handle:
        json.dump(output, handle, indent=4)


if __name__ == "__main__":
    args = parse_args()
    main(args)
