#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import statistics

# ------------------------------------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------------------------------------

# Motivation and functionality:

# This script is used to summarise Pyseer annotated kmers per locus
# First, Pyseer associated kmers were extracted using script filterPyseerOutput.r
# Second, Pyseer associated kmers were mapped to a set of reference genomes using annotate_hits_pyseer-runner.py
# Third, Pyseer mapped kmers were annotated using custom script annotate_kmers_pyseer.py to extract locus_tag, gene
#   and product information from embl files, including for intergenic regions
# And forth, this script is used to summarise Pyseer annotated kmers per locus


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to summarise Pyseer annotated kmers per locus.\n" \
                  "This script expects the output of script annotate_kmers_pyseer.py.\n"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-k", "--pyseer_kmers_ann_file", action="store", dest="pyseer_kmers_ann_file",
        help="Pyseer mapped and annotated kmers file obtained from running annotate_hits_pyseer-runner.py followed "
             "by annotate_kmers_pyseer.py",
        required=True, metavar="KMERS")
    group.add_argument(
        "-o", "--output_file", action="store", dest="output_file",
        help="Name of output file with Pyseer annotated kmers summarised by locus (including intergenic regions)",
        required=True, metavar="OUT")

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

    # Making sure input files exist
    logging.info(f'Making sure input files exist...')
    if not os.path.isfile(args.pyseer_kmers_ann_file):
        logging.error(f'Input file {args.pyseer_kmers_ann_file} not found!')
        sys.exit(-1)

    # symbols used in args.pyseer_kmers_ann_file to represent not annotated/found locus_tag
    locus_tag_na = ['-', '']

    # Saving annotated kmers along with mapped genome id and locus_tag
    logging.info(f'Savings kmers from {args.pyseer_kmers_ann_file}...')
    kmers_ref_locus = dict()  # dictionary to save all annotated kmers by reference genome and locus_tag
    # kmers_ref_locus{ref_id}{locus_tag} = [kmer_ids]
    all_kmers = dict()  # dictionary to save all kmer ids and original kmer line
    #     # all_kmers{kmer_id}{kmer_mapping} = kmer_line
    locus_gene = dict()  # dictionary to the gene names of locus tags
    locus_product = dict()  # dictionary to the products of locus tags
    with open(args.pyseer_kmers_ann_file, "r") as file:
        next(file)
        for line in file:
            kmer_line_items = line.strip().split('\t')
            kmer_id = kmer_line_items[0]
            kmer_mapping = kmer_line_items[-2]
            if kmer_id not in all_kmers:
                all_kmers[kmer_id] = dict()
                all_kmers[kmer_id][kmer_mapping] = line
            else:
                all_kmers[kmer_id][kmer_mapping] = line
            # extracting reference genome id
            # NOTE: although the same kmer can map to multiple locations on the same reference, each kmer mapping is
            # expected in a separate line
            ref_id = kmer_line_items[-2].split(":")[0]  # reference genome id
            # Extracting kmer annotations
            # NOTE: a kmer can have multiple locus annotations (separated by comma) because:
            #   - it maps to multiple locations on the reference genome
            #   - it spans multiple adjacent loci on the reference genome
            kmer_anns = kmer_line_items[-1].split(",")  # information on kmer locus annotation
            for kmer_ann in kmer_anns:
                kmer_ann_items = kmer_ann.split(";")
                locus_tag = kmer_ann_items[0]
                # if locus_tag == "|":
                #     print('locus_tag ' + locus_tag)
                #     print(line)
                #     exit(-1)
                if locus_tag not in locus_tag_na:
                    # saving gene and product information for locus_tag
                    gene = kmer_ann_items[1]
                    product = kmer_ann_items[2]
                    if ref_id not in locus_gene:
                        locus_gene[ref_id] = dict()
                        locus_product[ref_id] = dict()
                    locus_gene[ref_id][locus_tag] = gene
                    locus_product[ref_id][locus_tag] = product
                    # saving kmers for locus_tag
                    if ref_id not in kmers_ref_locus:
                        kmers_ref_locus[ref_id] = dict()
                    if locus_tag not in kmers_ref_locus[ref_id]:
                        kmers_ref_locus[ref_id][locus_tag] = list()
                    kmers_ref_locus[ref_id][locus_tag].append(kmer_id + "__" + kmer_mapping)
    print(str(len(all_kmers)) + " kmers saved in " + str(len(kmers_ref_locus)) +
          " different contigs (from assemblies)/chromosomes (from reference genomes)")
    logging.info(f'Savings kmers from {args.pyseer_kmers_ann_file}. DONE.')

    # Summarising annotated kmers per locus_tag
    header_items = ["ref_id", "median_pos", "min_pos", "max_pos", "locus_tag", "gene", "product", "num_kmers",
                    "min_pval", "max_pval", "mean_af", "min_af", "max_af", "mean_beta", "min_beta", "max_beta",
                    "lineages", "repeat_kmers", "num_repeat_kmers"]

    logging.info(f'Saving summarised kmers per locus...')
    header = '\t'.join(header_items) + '\n'
    output = open(args.output_file, 'w')
    output.write(header)
    for ref_id in kmers_ref_locus:
        for locus_tag in kmers_ref_locus[ref_id]:
            kmers_pos = []  # list to save the nucleotide coordinates of mapped kmers > to derive median
            gene = locus_gene[ref_id][locus_tag]
            if gene == "|":
                gene = "-"  # intergenic regions without gene names
            product = locus_product[ref_id][locus_tag]
            num_kmers = str(len(kmers_ref_locus[ref_id][locus_tag]))
            kmers_pval = []  # list to save the Pyseer lrt-pvalue of all locus kmers
            kmers_af = []  # list to save the Pyseer kmer frequence of all locus kmers
            kmers_beta = []  # list to save the Pyseer beta of all locus kmers
            kmers_lineages = []  # list to save the Pyseer lineages of all locus kmers
            repeat_kmers_dict = dict()  # dictionary with kmers in locus mapping to multiple locations in genome
            repeat_kmers = "no"
            num_repeat_kmers = "0"
            for kmer_mapping_id in kmers_ref_locus[ref_id][locus_tag]:
                (kmer_id, kmer_mapping) = kmer_mapping_id.split("__")
                kmer_line = all_kmers[kmer_id][kmer_mapping]
                # variant	af	filter-pvalue	lrt-pvalue	beta	beta-std-err	variant_h2	lineage	notes
                (var, af, fil_pval, lrt_pval, beta, beta_se, h2, lin, note, map, *_) = kmer_line.strip().split('\t')
                # NOTE: line below edited to accomodate other Pyseer outputs
                # (var, af, fil_pval, lrt_pval, beta, beta_se, h2, lin, map, *_) = kmer_line.strip().split('\t')
                kmers_pval.append(float(lrt_pval))
                kmers_af.append(float(af))
                kmers_beta.append(float(beta))
                kmers_lineages.append(str(lin))
                # adding kmer positions: e.g. MRSA252:2196930-2197029;;;
                # NOTE: only one kmer mapping expected per line
                print(map)
                (pos_from, pos_to) = map.strip().split(";")[0].split(":")[1].split("-")
                kmers_pos.extend(list(range(int(pos_from), int(pos_to), 1)))
                # if kmer has multiple mappings, then save as repeat_kmer
                if len(all_kmers[kmer_id].keys()) > 1:
                    repeat_kmers_dict[kmer_id] = "yes"
            kmers_pos = list(set(kmers_pos))
            median_pos = str(statistics.median(kmers_pos))
            min_pos = str(min(kmers_pos))
            max_pos = str(max(kmers_pos))
            min_pval = str(min(kmers_pval))
            max_pval = str(max(kmers_pval))
            mean_af = str(statistics.mean(kmers_af))
            min_af = str(min(kmers_af))
            max_af = str(max(kmers_af))
            mean_beta = str(statistics.mean(kmers_beta))
            min_beta = str(min(kmers_beta))
            max_beta = str(max(kmers_beta))
            lineages = '-'.join(list(set(kmers_lineages)))
            num_repeat_kmers = str(len(repeat_kmers_dict.keys()))
            if len(repeat_kmers_dict.keys()) > 0:
                repeat_kmers = "yes"
            newline_items = [ref_id, median_pos, min_pos, max_pos, locus_tag, gene, product, num_kmers, min_pval,
                             max_pval, mean_af, min_af, max_af, mean_beta, min_beta, max_beta, lineages, repeat_kmers,
                             num_repeat_kmers]
            print('\t'.join(newline_items))
            output.write('\t'.join(newline_items) + "\n")
    logging.info(f'Saving summarised kmers per locus. DONE.')
    output.close()
    # to do: saved not annotated kmers
    

if __name__ == "__main__":
    _main()