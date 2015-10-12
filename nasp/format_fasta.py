#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1.0.0"
__email__ = "dsmith@tgen.org"

import logging


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Reformats a fasta to be split 80 characters per line, with system line-endings.")
    parser.add_argument("--inputfasta", required=True, help="Path to input fasta.")
    parser.add_argument("--outputfasta", required=True, help="Path to output fasta.")
    return parser.parse_args()

def format_fasta(inputfasta, outputfasta):
    from nasp.nasp_objects import Genome
    fasta_data = Genome()
    fasta_data.import_fasta_file(inputfasta)
    fasta_data.write_to_fasta_file(outputfasta)


def main():
    commandline_args = _parse_args()
    format_fasta(commandline_args.inputfasta, commandline_args.outputfasta)


if __name__ == "__main__":
    main()


