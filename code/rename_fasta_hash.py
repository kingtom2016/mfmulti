#!/usr/bin/env python3

"""
Author: Matthew Kweskin, github: @mkweskin
This script converts the sequence names in a FASTA file
to the md5 hash of the sequence.
Notes:
- Gaps will automatically be removed by biopython
"""

import hashlib
import argparse
from os import path
import sys
import re
from datetime import datetime

if (sys.version_info < (3, 0)):
    raise Exception("This script requires python3. One way to install this is with miniconda3: https://docs.conda.io/en/latest/miniconda.html")

try:
    from Bio import SeqIO
except ImportError:
    raise Exception("This script requires Biopython. One way to install it is with conda: conda install -c conda-forge biopython")

def get_args():
    parser = argparse.ArgumentParser(description="Convert names of a FASTA file to the md5 of their sequence", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    currdatetime = datetime.now()
    parser.add_argument("INPUT", help = 'Input FASTA file')
    parser.add_argument("OUTPUT", help = 'Output FASTA file')
    parser.add_argument("--overwrite",  action="store_true", help = 'Automatically overwrite the output file and log file if they\'re already present')
    parser.add_argument("--case",  choices=['upper', 'keep'], help = 'Convert sequence to all upper case (\'upper\') or maintain case of the original file (\'keep\')', default = 'upper')
    parser.add_argument("--outfmt", choices=['multi', 'single'], help = 'Format of output FASTA file: \'single\' for sequence to be on one line, \'multi\' for a multi-line FASTA file.', default = 'single')
    parser.add_argument("--log", help = 'Log file of a tab separated list of old and sequence names.',default = 'header_hash.log') #default = currdatetime.strftime("%Y%m%dT%H%M%S.%f") + '.log')
    parser.add_argument("--nolog", action="store_true", help = 'Disable the log file')
    return parser.parse_args()

def file_checks(args):
    """
    Some checks on the input/output/log files
    """
    if not path.exists(args.INPUT):
        raise Exception("Input file does not exist")

    if path.exists(args.OUTPUT) and not args.overwrite:
        answer = input("Output file exists. (use --overwrite to Automatically overwrite the output file)\nOverwrite [y/n]? ")
        if answer.lower() != 'y':
            print ("Exiting script.")
            sys.exit()

    if (not args.nolog) and path.exists(args.log) and not args.overwrite:
        answer = input("Log file exists. (use --overwrite to Automatically overwrite the output file)\nOverwrite [y/n]? ")
        if answer.lower() != 'y':
            print ("Exiting script.")
            sys.exit()

    logdir = path.dirname(args.log)
    if not path.exists(logdir) and logdir != "":
        raise Exception("Log directory does not exist")

    outputdir = path.dirname(args.OUTPUT)
    if not path.exists(outputdir) and outputdir != "":
        raise Exception("Output directory does not exist")


def convert_md5(original, modified, log, args):
    """
    Does the conversion of sequnce names to md5 hashes
    """
    if args.outfmt == 'single' or args.outfmt is None:
        outformat = 'fasta-2line'
    else:
        outformat = 'fasta'

    records = SeqIO.parse(original, "fasta")
    for record in records:
        seq = record.seq
        if args.case == 'upper':
            seq = seq.upper()
        seq_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
        orig_description = record.description
        record.id = seq_md5
        record.description = seq_md5
        if log != None:
            log.write("\t".join([orig_description, seq_md5]) + "\n")
        SeqIO.write(record, modified, outformat)

def main():
    args = get_args()
    file_checks(args)

    if args.nolog:
        with open(args.INPUT) as original, open(args.OUTPUT, 'w') as modified:
            convert_md5(original, modified, None, args)
    else:
        with open(args.INPUT) as original, open(args.OUTPUT, 'w') as modified, open(args.log, 'w') as log:
            convert_md5(original, modified, log, args)

if __name__ == '__main__':
    main()