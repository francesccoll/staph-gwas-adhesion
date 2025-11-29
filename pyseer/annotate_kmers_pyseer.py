#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ------------------------------------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------------------------------------

# Motivation and functionality:

# This script is used to extract gene annotation for Pyseer associated kmers. Pyseer's script
# annotate_hits_pyseer-runner.py needs to be run first to map kmers to a set of reference genomes and assemblies.
# annotate_hits_pyseer-runner.py only outputs chromosome id and kmer coordinates. This script is used to extract kmer
# gene annotation from this output file and assign his to kmers.
#
# NOTE: two types of mapped genomes are allowed (reference and assemblies).
# Only Bakta-generated .embl formatted genome annotations supported.

# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to extract gene annotation for Pyseer associated kmers.\n" \
                  "Pyseer's script annotate_hits_pyseer-runner.py needs to be run first.\n"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-k", "--pyseer_kmers_file", action="store", dest="pyseer_kmers_file",
        help="Pyseer kmers file obtained from annotate_hits_pyseer-runner.py with kmers mapped",
        required=True, metavar="KMERS")
    group.add_argument(
        "-g", "--map_genomes", action="store", dest="map_genomes",
        help="Tab-delimited file with absolute paths to FASTA and annotated genomes used for mapping kmers "
             "(i.e. 'references' positional argument in script annotate_hits_pyseer-runner.py)",
        required=True, metavar="GENOMES")
    group.add_argument(
        "-o", "--output_file", action="store", dest="output_file",
        help="Name of output with Pyseer kmers file extended with gene function information",
        required=True, metavar="OUT")
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-c", "--min_contig_length", action="store", dest="min_contig_length",
        help="Minimum length of contigs to keep mapped kmers (default: 1000)",
        required=False, metavar="CONTIG", type=float, default=1000)

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

    # Other arguments
    feature_types = ['CDS', 'inter']
    # Using 'ignore_features' is probably more inclusive
    ignore_features = ['source', 'gene']
    min_contig_length = args.min_contig_length

    # Making sure input files exist
    logging.info(f'Making sure input files exist...')
    input_files = [args.pyseer_kmers_file, args.map_genomes]
    for input_file in input_files:
        if not os.path.isfile(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    map_genome_ann = dict()  # dictionary to save genome annotation file for each genome used for kmer mapping
    map_genome_type = dict()  # dictionary to save type of genome used for kmer mapping
    map_genome_len = dict()  # dictionary to save length of chromosome/contig in reference FASTA files
    for line in open(args.map_genomes, "r"):
        (fasta_file, ann_file, genome_type) = line.strip().split('\t')
        if not os.path.isfile(fasta_file):
            logging.error(f'FASTA file {fasta_file} not found!')
            sys.exit(-1)
        if not os.path.isfile(ann_file):
            logging.error(f'Annotation file {ann_file} not found!')
            sys.exit(-1)
        # Saving FASTA header (the one in pyseer_kmers_file) along with
        input_records = SeqIO.parse(fasta_file, "fasta")
        for record in input_records:
            # print(record.id)
            map_genome_ann[record.id] = ann_file
            map_genome_type[record.id] = genome_type
            map_genome_len[record.id] = len(record.seq)
    logging.info(f'Making sure input files exist. DONE')

    # Saving mapped kmers along with mapped genome id and position
    logging.info(f'Reading kmers from {args.pyseer_kmers_file}...')
    kmers_ref_pos = dict()  # dictionary to save all kmer coordinates by reference genome
    # kmers_ref_pos{ref_id}{pos} = [kmer_ids]
    all_kmers = dict()  # dictionary to save all kmer ids and original kmer line
    # all_kmers{kmer_id}{mapping} = kmer_line
    # NOTE: each kmer can have multiple 'mappings', that is, it can map to multiple locations in the genome;
    #   and each kmer mapping multiple kmer locus annotations (as one single kmer can span multiple adjacent features)
    for line in open(args.pyseer_kmers_file, "r"):
        kmer_line_items = line.strip().split('\t')
        kmer_id = kmer_line_items[0]
        all_kmers[kmer_id] = dict()
        # here ref_id refers to the FASTA header(s) in assemblies/complete genome (e.g. contig ids)
        ref_id = line.strip().split('\t')[-1].split(":")[0]
        # if contig length longer than chosen, then keep kmers
        if map_genome_len[ref_id] >= min_contig_length:
            kmer_mappings = kmer_line_items[-1].split(",")
            for kmer_mapping in kmer_mappings:
                kmer_mapping_items = kmer_mapping.split(";")
                # Saving kmer mappings
                kmer_line_mapping = '\t'.join(kmer_line_items[:-1]) + '\t' + kmer_mapping
                all_kmers[kmer_id][kmer_mapping] = kmer_line_mapping
                (pos_from, pos_to) = kmer_mapping_items[0].split(":")[1].split("-")
                for pos in range(int(pos_from), int(pos_to)):
                    if ref_id not in kmers_ref_pos:
                        kmers_ref_pos[ref_id] = dict()
                    if pos not in kmers_ref_pos[ref_id]:
                        kmers_ref_pos[ref_id][pos] = list()
                    kmers_ref_pos[ref_id][pos].append(kmer_id)
                    # print("kmers_ref_pos["+ref_id+"]["+str(pos)+"] > " + str(kmers_ref_pos[ref_id][pos]))
    print(str(len(all_kmers)) + " kmers saved in " + str(len(kmers_ref_pos)) +
          " different contigs (from assemblies)/chromosomes (from reference genomes)")
    print("Kmers mapped on: " + '\t\n'.join(kmers_ref_pos.keys()))  # prints keys
    logging.info(f'Reading kmers from {args.pyseer_kmers_file}. DONE.')

    # Extracting gene annotation information for kmers
    # Loop structure below:
    #   For each reference/contig
    #       For each annotated feature
    #           > extract feature annotation: locus_tag;gene;product
    #           For all nucleotide positions in annotated feature
    #               > Extract all kmers mapped onto these positions and assign feature annotation to kmer ids
    logging.info(f'Extracting gene annotation information for kmers...')
    kmers_ann = dict()  # dictionary to save kmer annotations, list of locus_tag;gene;product
    # kmers_ann{kmer_id}{mapping_id} = kmer_mapping_ann
    for ref_id in kmers_ref_pos:
        print("ref_id: " + ref_id)
        if ref_id not in map_genome_ann:
            logging.error(f'{ref_id} not found in {args.map_genomes}, please check.')
            exit(-1)
        ann_file = map_genome_ann[ref_id]
        input_records = SeqIO.parse(ann_file, "embl")
        for record in input_records:
            use_record = False
            # NOTE: if mapped genome is a reference genome, then use single expected record/chromosome
            if map_genome_type[ref_id] == "ref":
                use_record = True
            if map_genome_type[ref_id] == "draft":
                # expected contig ids of draft assemblies on pyseer_kmers_file: sample.contigX
                kmer_contig = ref_id.split(".")[1]
                # Bakta-annotated genome files my contain '_' in contig ids
                record_contig = record.id.replace('_', '')
                if kmer_contig == record_contig:
                    use_record = True
            if use_record:
                print(record.id)
                for feature in record.features:
                    # if feature.type in feature_types:
                    if feature.type not in ignore_features:
                        # NOTE: neither of the CDS feature attributes extracted can contain ';'
                        # as this is used as a separator; ',' is also used a separator as the same kmer
                        # can span multiple features
                        locus_tag = "-"
                        if "locus_tag" in feature.qualifiers:
                            locus_tag = str(feature.qualifiers["locus_tag"][0]).replace(";", "_").replace(",", " ")
                        gene_name = "-"
                        if "gene" in feature.qualifiers:
                            gene_name = str(feature.qualifiers["gene"][0]).replace(";", "_").replace(",", " ")
                        product = "-"
                        if "product" in feature.qualifiers:
                            product = str(feature.qualifiers["product"][0]).replace(";", "_").replace(",", " ")
                        gene_ann = locus_tag + ";" + gene_name + ";" + product
                        if locus_tag == "|":
                            print('locus_tag ' + locus_tag)
                            print(feature)
                            # exit(-1)
                        # NOTE: a feature can contain a CompoundLocation, i.e. multiple FeatureLocations
                        for feature_location in feature.location.parts:
                            for pos in range(int(feature_location.start), int(feature_location.end)):
                                if pos in kmers_ref_pos[ref_id]:
                                    kmer_ids = kmers_ref_pos[ref_id][pos]
                                    for kmer_id in kmer_ids:
                                        if kmer_id not in kmers_ann:
                                            kmers_ann[kmer_id] = dict()
                                        # making sure locus annotation is assigned to correct kmer_mapping
                                        kmer_mappings_tmp = list()
                                        for kmer_mapp in all_kmers[kmer_id]:
                                            #print("kmer_mapp " + kmer_mapp)
                                            # example of kmer_mapping: MRSA252:644058-644118;;;
                                            (kmer_mapp_l, kmer_mapp_r) = kmer_mapp.split(";")[0].split(":")[1].split("-")
                                            kmer_mapp_range = list(range(int(kmer_mapp_l), int(kmer_mapp_r)))
                                            #print("kmer_mapp_range " + str(kmer_mapp_range))
                                            if pos in kmer_mapp_range:
                                                kmer_mappings_tmp.append(kmer_mapp)
                                        #print("kmer_mappings_tmp " + str(kmer_mappings_tmp))
                                        kmer_mapping = str(list(set(kmer_mappings_tmp))[0])
                                        #print("kmer_mapping " + str(kmer_mapping))
                                        if kmer_mapping not in kmers_ann[kmer_id]:
                                            kmers_ann[kmer_id][kmer_mapping] = list()
                                        else:
                                            kmers_ann[kmer_id][kmer_mapping].append(gene_ann)
                                        # getting unique list of annotations
                                        kmers_ann[kmer_id][kmer_mapping] = list(set(kmers_ann[kmer_id][kmer_mapping]))
                                        # print("kmers_ann["+kmer_id+"]["+kmer_mapping+"] > " + str(kmers_ann[kmer_id][kmer_mapping]))
    logging.info(f'Extracting gene annotation information for kmers. DONE.')

    logging.info(f'Saving gene annotation for kmers...')
    output = open(args.output_file, 'w')
    header_items = ["variant", "af", "filter.pvalue", "lrt.pvalue", "beta", "beta.std.err", "variant_h2", "lineage",
                    "notes", "kmer_map", "kmer_ann"]
    header = '\t'.join(header_items) + '\n'
    output.write(header)
    for kmer_id in all_kmers:
        for kmer_mapping in all_kmers[kmer_id]:
            line = all_kmers[kmer_id][kmer_mapping]
            ann = "-"
            if kmer_id in kmers_ann:
                if kmer_mapping in kmers_ann[kmer_id]:
                    ann = ','.join(kmers_ann[kmer_id][kmer_mapping])
            output_line = line.strip() + "\t" + ann + "\n"
            output.write(output_line)
    logging.info(f'Saving gene annotation for kmers. DONE.')
    output.close()


if __name__ == "__main__":
    _main()