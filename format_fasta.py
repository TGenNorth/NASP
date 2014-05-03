#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1"
__email__ = "dsmith@tgen.org"


def _parse_args():
    import argparse
    parser = argparse.ArgumentParser( description="Reformats a fasta to be split 80 characters per line, with system line-endings." )
    parser.add_argument( "--inputfasta", required=True, help="Path to input fasta." )
    parser.add_argument( "--outputfasta", required=True, help="Path to output fasta." )
    return parser.parse_args()

def main():
    from nasp_objects import Genome
    commandline_args = _parse_args()
    fasta_data = Genome()
    fasta_data.import_fasta_file( commandline_args.inputfasta )
    fasta_data.write_to_fasta_file( commandline_args.outputfasta )


if __name__ == "__main__": main()


