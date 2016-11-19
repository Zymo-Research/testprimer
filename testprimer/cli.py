"""Command line interface(CLI).

Command line interface for testprimer.
"""

from __future__ import print_function, absolute_import

import os
import sys
import argparse

from testprimer.pcr import pcr
from testprimer.report import report


def main():
    parser = argparse.ArgumentParser(
        description="Zymo's internal utility to evaluate the performance of \
                     primer pools by running in silico PCR on the given \
                     microbial database."
    )
    subparsers = parser.add_subparsers(
        dest='subcommand',
        help='sub-command help'
    )
 
    pcr_parser = subparsers.add_parser(
        'pcr',
        help='Run in silico PCR on sequences in FASTA format with forward \
              primer pool and reverse primer pool.'
    )
    pcr_parser.add_argument(
        '-a', '--fasta',
        dest='fasta_path',
        metavar='',
        required=True,
        help='Path to the FASTA file containing aligned microbial genomes.'
    )
    pcr_parser.add_argument(
        '-f', '--forward',
        dest='fw_path',
        metavar='',
        required=True,
        help='Path to the file containing forward primer pool of which \
              primers have the same coordinates.'
    )
    pcr_parser.add_argument(
        '-r', '--reverse',
        dest='rv_path',
        metavar='',
        required=True,
        help='Path to the file containing reverse primer pool of which \
              primers have the same coordinates.'
    )
    pcr_parser.add_argument(
        '-n', '--name',
        dest='filename',
        metavar='',
        default=None,
        help='Output SQL file name. Please include extension. \
              [default: <fasta>.<fw>.<rv>.sql]'
    )
    pcr_parser.add_argument(
        '-o',
        dest='out_dir',
        metavar='',
        default=os.getcwd(),
        help='Directory file outputs to. [default: <cwd>]'
    )
    pcr_parser.add_argument(
        '--coverage',
        dest='taxa_coverage',
        action='store_true',
        help='Perform taxa coverage analysis jointly after in silico PCR is \
              complete and output result in Excel format.'
    )

    report_parser = subparsers.add_parser(
        'report',
        help='Perform analysis and generate report in various forms based on \
              in silico PCR result.'
    )
    report_parser.add_argument(
        '-s', '--sql',
        dest='sql_path',
        metavar='',
        required=True,
        help='Path to the SQL file containing in silico PCR result.'
    )
    report_parser.add_argument(
        '-o',
        dest='out_dir',
        metavar='',
        default=os.getcwd(),
        help='Directory file outputs to. [default: cwd]'
    )
    report_parser.add_argument(
        '--coverage',
        dest='taxa_coverage',
        action='store_true',
        help='Flag to perform taxa coverage analysis and output result in \
              Excel format.'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    if len(sys.argv) == 2:
        if sys.argv[1] == 'pcr':
            pcr_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == 'report':
            report_parser.print_help()
            sys.exit(0)
        else:
            # Python automatically prints argparse's built-in error message
            # for wrong subcommands.
            pass

    args = parser.parse_args()

    if args.subcommand == 'pcr':
        if not os.path.isfile(args.fasta_path):
            print("Error: FASTA file '%s' not found." % args.fasta_path)
            sys.exit(1)
        if not os.path.isfile(args.fw_path):
            print("Error: forward primer pool '%s' not found." % args.fw_path)
            sys.exit(1)
        if not os.path.isfile(args.rv_path):
            print("Error: reverse primer pool '%s' not found." % args.rv_path)
            sys.exit(1)

        if os.path.splitext(args.fasta_path)[-1] == '.gz':
            print("Error: Please unzip FASTA file first.")
            sys.exit(1)

        if not args.filename:
            filename = '%s.%s.%s.sql' % (
                os.path.splitext(os.path.basename(args.fasta_path))[0],
                os.path.splitext(os.path.basename(args.fw_path))[0],
                os.path.splitext(os.path.basename(args.rv_path))[0],
            )
        else:
            filename = args.filename
            
        try:
            pcr(args.fasta_path, args.fw_path, args.rv_path, filename, args.out_dir)
            report(os.path.join(args.out_dir, filename), args.out_dir, args.taxa_coverage)
        except KeyboardInterrupt:
            print("Error: Execution aborted.")
            sys.exit(1)

    if args.subcommand == 'report':
        if not os.path.isfile(args.sql_path):
            print("Error: SQL file '%s' does not exist." % args.sql_path)
            sys.exit(1)

        try:
            report(args.sql_path, args.out_dir, args.taxa_coverage)
        except KeyboardInterrupt as e:
            print("Error: Execution aborted.")
            sys.exit(1)
        

if __name__ == '__main__':
    main()
