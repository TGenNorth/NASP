#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1.0.0"
__email__ = "dsmith@tgen.org"

import logging


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--nucmerpath", default="nucmer", help="Path to the 'nucmer' executable.")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file.")
    return parser.parse_args()


# TODO: This should eventually be moved to the main job manager section
def run_nucmer_on_reference(nucmer_path, reference_path):
    import subprocess
    import sys

    return_code = subprocess.call(
        [nucmer_path, "--prefix=reference", "--maxmatch", "--nosimplify", reference_path, reference_path])
    if return_code > 0:
        sys.stderr.write(
            "NASP WARNING: nucmer may have encountered errors during reference duplicates checking, proceeding anyway\n")


def _parse_delta_line(line_from_delta_file, dups_data, current_contigs):
    import re

    line_match = re.match(r'^>([^ ]+) ([^ ]+) (\d+) (\d+)\s*$', line_from_delta_file)
    if line_match:
        current_contigs = ( line_match.group(1), line_match.group(2) )
        contig_0_end = int(line_match.group(3))
        contig_1_end = int(line_match.group(4))
        dups_data.add_contig(current_contigs[0])
        dups_data.add_contig(current_contigs[1])
        if contig_0_end > dups_data.get_contig_length(current_contigs[0]):
            dups_data.extend_contig(contig_0_end, "0", current_contigs[0])
        if contig_0_end > dups_data.get_contig_length(current_contigs[1]):
            dups_data.extend_contig(contig_1_end, "0", current_contigs[1])
    else:
        line_match = re.match(r'^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$', line_from_delta_file)
        if line_match:
            contig_0_start = int(line_match.group(1))
            contig_0_end = int(line_match.group(2))
            contig_1_start = int(line_match.group(3))
            contig_1_end = int(line_match.group(4))
            if ( current_contigs[0] != current_contigs[1] ) or (
                ( contig_0_start != contig_1_start ) and ( contig_0_end != contig_1_end ) ):
                if contig_0_end < contig_0_start:
                    contig_0_end, contig_0_start = contig_0_start, contig_0_end
                if contig_1_end < contig_1_start:
                    contig_1_end, contig_1_start = contig_1_start, contig_1_end
                dups_data.set_value(( "1" * ( contig_0_end - contig_0_start + 1 ) ), contig_0_start, "!",
                                    current_contigs[0])
                dups_data.set_value(( "1" * ( contig_1_end - contig_1_start + 1 ) ), contig_1_start, "!",
                                    current_contigs[1])
    return current_contigs


def parse_delta_file(delta_filename, dups_data):
    current_contigs = ( "", "" )
    delta_handle = open(delta_filename, 'r')
    for line_from_delta_file in delta_handle:
        current_contigs = _parse_delta_line(line_from_delta_file, dups_data, current_contigs)
    delta_handle.close()


def main():
    from nasp.nasp_objects import GenomeStatus

    commandline_args = _parse_args()
    run_nucmer_on_reference(commandline_args.nucmerpath, commandline_args.reference)
    dups_data = GenomeStatus()
    parse_delta_file("reference.delta", dups_data)
    dups_data.write_to_fasta_file("duplicates.txt")


if __name__ == "__main__":
    main()


