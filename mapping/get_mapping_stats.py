#!/usr/bin/env python3

import string
import argparse
import logging
import os
import subprocess
import sys
import statistics

# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------
# Script tested with the following dependency versions:
#   samtools v1.9 (using htslib 1.9)
#   bcftools v1.9 (using htslib 1.9)
# Tested with snps.raw.vcf produced by snippy v4.6.0
# Tested with snps.bam produced by snippy v4.6.0
#
# To do:
#   make sure samtools v1.9 and bcftools v1.9 are locally installed, otherwise num_ht_sites may not be calculated

# ------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------

_DEPENDENCIES = ['samtools', 'bcftools']

# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to obtain short-read mapping statistics from BAM and VCF files for QC purposes"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-b", "--bam_file", action="store", dest="bam_file",
        help="BAM file obtained from short-read mapping pipeline", required=True, metavar="BAM_FILE")
    group.add_argument(
        "-v", "--vcf_file", action="store", dest="vcf_file",
        help="VCF file (uncompressed) obtained from short-read mapping pipeline", required=True, metavar="VCF_FILE")
    group.add_argument(
        "-s", "--sample_id", action="store", dest="sample_id",
        help="sample id used as prefix to name output files", required=True, metavar="SAMPLE_ID")

    group = parser.add_argument_group('optional arguments')
    group.add_argument(
        "-d", "--delete_tmp", action="store", dest="delete_tmp",
        help="delete temporary files", required=False, default=True, metavar="DELETE_TMP")
    group.add_argument(
        '--version', action='version', version='%(prog)s 1.0')

    return parser.parse_args()


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    output = subprocess.check_output(['which', executable_name]).strip()
    if output:
        found = True
    return found


def check_file_exist(file, file_tag):
    if not os.path.isfile(file):
        logging.error(f'{file_tag} {file} not found!')
        sys.exit(-1)
    else:
        logging.info(f'{file_tag} {file} found!')


def run_command_shell_string(command_line_string):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line_string: it must be a string not a list
    """
    print('\tRunning: ' + command_line_string)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)


def run_command_string(command_line_string):
    """
    This function executes a command line, check for execution errors and returns stdout
    :param command_line_string: it must be a string
    :return: stdout
    """
    print('\tRunning: ' + command_line_string)
    try:
        process_completed = subprocess.run(
            command_line_string,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            shell=True,
        )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)
    return process_completed.stdout.decode('utf-8')


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

    # Making sure dependencies are installed
    logging.info('Making sure dependencies exist...')
    for dependency in _DEPENDENCIES:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure required input files exist
    check_file_exist(args.bam_file, "Input BAM file")
    check_file_exist(args.vcf_file, "Input VCF file")

    # Extracting mapping stats from bam file
    output_file1 = args.sample_id + '.samtools_stats.txt'
    run_command_shell_string(''.join(['samtools stats ', args.bam_file, ' | grep ^SN | cut -f 2- > ', output_file1]))
    check_file_exist(output_file1, "samtools stats output")

    # Parsing samtools stats output file to save stats
    samtools_stats = dict()
    for line in open(output_file1, 'r'):
        items = line.strip().split('\t')
        variable = items[0].replace(':', '')
        value = items[1]
        samtools_stats[variable] = value

    # Extracting depth of coverage stats from BAM file
    output_file2 = args.sample_id + '.samtools_depth.txt'
    run_command_shell_string(''.join(['samtools depth ', args.bam_file, ' > ', output_file2]))
    check_file_exist(output_file2, "samtools depth output")

    # Extracting chromosome length from VCF header
    # e.g. ##contig=<ID=BX571856,length=2902619>
    vcf_line = run_command_string(''.join(['cat ', args.vcf_file, ' | grep "^##contig"']))
    chr_length = vcf_line.strip().split(',')[1].replace('length=', '').replace('>', '')
    logging.info(f'chr_length extracted: {chr_length}')

    # Calculating average and stdev depth of coverage
    chr_depth = []
    for line in open(output_file2, 'r'):
        n_depth = int(line.strip().split('\t')[2])
        chr_depth.append(n_depth)
    avg_depth = sum(chr_depth)/int(chr_length)
    logging.info(f'avg_depth extracted: {str(avg_depth)}')
    stdev_depth = statistics.stdev(chr_depth)
    logging.info(f'stdev_depth extracted: {str(stdev_depth)}')

    # Extracting number of heterozygous sites from VCF file
    num_ht_sites = run_command_string(''.join(["bcftools view ",
                                               args.vcf_file,
                                               " | awk -F'\t' '{ print $8}' | "
                                               "awk -F';' '{for (i=1;i<=NF;i++) {if ($i ~ /^AF=/) {print $i}}}' | "
                                               "sed 's/AF=//g' | "
                                               "awk -F'\t' '{ if($1>0.2 && $1<0.8) print $0}' | wc -l"]))
    logging.info(f'num_ht_sites extracted: {str(num_ht_sites)}')

    # Saving output file
    final_output_file = args.sample_id + '.mapping_qc_stats.csv'
    logging.info(f'Writing mapping statistics in {final_output_file}')
    header = 'sample_id'
    newline = args.sample_id
    for variable in samtools_stats:
        header += '\t' + variable
        newline += '\t' + samtools_stats[variable]
    header += '\t' + 'avg_depth' + '\t' + 'stdev_depth' + '\t' + 'num_ht_sites' + '\n'
    newline += '\t' + str(avg_depth) + '\t' + str(stdev_depth) + '\t' + str(num_ht_sites) + '\n'
    output = open(final_output_file, 'w')
    output.write(header)
    output.write(newline)
    output.close()

    # Removing temporary files
    logging.info('Removing temporary files.')
    rm_command = ['rm', output_file1, output_file2]
    run_command_shell_string(' '.join(rm_command))


if __name__ == "__main__":
    _main()
