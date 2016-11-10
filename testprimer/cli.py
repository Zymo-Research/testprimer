# from __future__ import absolute_import

# from pcr import PCRArray
# from report import Analysis, TaxaCoverage


# def main(fasta_path, fw_path, rv_path, taxa_coverage, out_dir):
    # pcrarray = PCRArray(fasta_path, fw_path, rv_path)
    # pcrarray.to_sql('filename.sql', out_dir)

    # if taxa_coverage:
        # taxa_analyzer = TaxaAnalyzer()
        # analysis = Analysis(sql_path, taxa_analyzer, out_dir)
        # analysis.execute()


if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Zymo's internal utility to evaluate the performance of \
                     primer pools by running in silico PCR on the given \
                     microbial database."
    )
    subparsers = parser.add_subparsers(help='sub-command help')

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
        help='Path to the file containing forward primers that have the same \
              coordinates.'
    )
    pcr_parser.add_argument(
        '-r', '--reverse',
        dest='rv_path',
        metavar='',
        required=True,
        help='Path to the file containing reverse primers that have the same \
              coordinates.'
    )
    pcr_parser.add_argument(
        '--coverage',
        dest='taxa_coverage',
        # metavar='',
        action='store_true',
        default=False,
        help='Perform taxa coverage analysis jointly after in silico PCR is \
              complete and output result in Excel format.'
    )
    pcr_parser.add_argument(
        '-o',
        dest='out_dir',
        metavar='',
        help='Directory file outputs to.'
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
        '--coverage',
        dest='analyzer',
        # metavar='',
        action='store_true',
        default=False,
        help='Flag to perform taxa coverage analysis and output result in \
              Excel format.'
    )
    report_parser.add_argument(
        '-o',
        dest='out_dir',
        metavar='',
        help='Directory file outputs to.'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if sys.argv[1] == pcr_parser.prog.split()[1]:
            pcr_parser.print_help()
            sys.exit(1)
        elif sys.argv[1] == report_parser.prog.split()[1]:
            report_parser.print_help()
            sys.exit(1)
        else:
            pass
    else:
        pass

    args = parser.parse_args()
