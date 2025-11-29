#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to replace FASTA Ids in a multi-FASTA files"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i", "--input_fasta", action="store", dest="input_fasta",
        help="Input multi-FASTA file",
        required=True, default=True, metavar="INPUT")
    parser.add_argument(
        "-t", "--ids_table", action="store", dest="ids_table",
        help="Two column table with current (first column) and new (second column) to be replaced",
        required=True, default=True, metavar="IDS_TABLE")
    parser.add_argument(
        "-o", "--output_fasta", action="store", dest="output_fasta",
        help="Output FASTA with Ids replaced",
        required=True, default=True, metavar="OUTPUT")

    return parser.parse_args()


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------

def _main():

    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Dictionary to save old ids as keys and new Ids as values
    replace_ids = dict()

    num_input_ids = 0  # number of input FASTA sequences
    num_ids_replaced = 0  # number of ids replaced

    if not os.path.isfile(args.ids_table):
        logging.error(f"Input Ids tables {args.ids_table} not found")
        exit(-1)
    else:
        with open(args.ids_table, "r") as id_file:
            for line in id_file:
                ids = line.strip().split('\t')
                replace_ids[ids[0]] = ids[1]

    if not os.path.isfile(args.input_fasta):
        logging.error(f"Input FASTA file {args.input_fasta} not found")
        exit(-1)
    else:
        logging.info(f"Reading input FASTA file {args.input_fasta}")
        output_records = []
        for record in SeqIO.parse(args.input_fasta, "fasta"):
            num_input_ids += 1
            if record.id in replace_ids:
                num_ids_replaced += 1
                record.id = record.description = replace_ids[record.id]
            else:
                logging.warning(f"Id in input FATSA {record.id} not found")
            output_records.append(record)
        logging.info(f"Saving output FASTA file {args.output_fasta}")
        SeqIO.write(output_records, args.output_fasta, "fasta")

    logging.info(f"{num_ids_replaced}/{num_input_ids} input Ids replaced")


if __name__ == "__main__":
    _main()
