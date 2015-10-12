#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1.0.0"
__email__ = "dsmith@tgen.org"

import logging


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--nucmerpath", default="nucmer", help="Path to the 'nucmer' executable.")
    parser.add_argument("--nucmerargs", default="", help="Optional arguments to pass to the 'nucmer' executable.")
    parser.add_argument("--deltafilterpath", default="delta-filter", help="Path to the 'delta-filter' executable.")
    parser.add_argument("--deltafilterargs", default="", help="Optional arguments to pass to the 'delta-filter' executable.")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file.")
    parser.add_argument("--external", required=True, help="Path to the external genome fasta file.")
    parser.add_argument("--name", default="", help="Name of this external genome.")
    return parser.parse_args()


# This should eventually be moved to the main job manager section
def generate_delta_file(nucmer_path, nucmer_args, delta_filter_path, delta_filter_args, external_nickname, reference_path, external_path):
    import subprocess
    import sys

    return_code = subprocess.call(
        [nucmer_path, ( "--prefix=" + external_nickname ), reference_path, external_path] + nucmer_args.split())
    if return_code > 0:
        sys.stderr.write(
            "NASP WARNING: nucmer may have encountered errors while processing external genome '" + external_nickname + "', proceeding anyway\n")
    filtered_delta_handle = open(( external_nickname + ".filtered.delta" ), "w")

    return_code = subprocess.call([delta_filter_path, "-q", "-r", "-o", "100"] + delta_filter_args.split() + [(external_nickname + ".delta")],
                                  stdout=filtered_delta_handle)
    if return_code > 0:
        sys.stderr.write(
            "NASP WARNING: delta-filter may have encountered errors while processing external genome '" + external_nickname + "', proceeding anyway\n")
    filtered_delta_handle.close()


def _update_genome_from_delta_data(franken_genome, external_genome, parser_state, distance_covered, is_external_insert):
    from nasp.nasp_objects import Genome

    if distance_covered == -1:
        distance_covered = parser_state['final_pos'] - parser_state['reference_pos'] + 1
        is_external_insert = True
    if distance_covered > 0:
        if parser_state['external_is_reversed']:
            matching_segment = Genome.reverse_complement(''.join(
                external_genome.get_call(( parser_state['external_pos'] - distance_covered + 1 ),
                                         parser_state['external_pos'])))
        else:
            matching_segment = ''.join(external_genome.get_call(parser_state['external_pos'], (
                parser_state['external_pos'] + distance_covered - 1 )))
        franken_genome.set_call(matching_segment, parser_state['reference_pos'], 'X')
    parser_state['reference_pos'] = parser_state['reference_pos'] + distance_covered
    parser_state['external_pos'] = parser_state['external_pos'] + (
        -distance_covered if parser_state['external_is_reversed'] else distance_covered )
    if is_external_insert:
        parser_state['external_pos'] += -1 if parser_state['external_is_reversed'] else 1
    else:
        franken_genome.set_call('.', parser_state['reference_pos'], '!')
        parser_state['reference_pos'] += 1
    return parser_state


def _parse_delta_line(line_from_delta_file, franken_genome, external_genome, parser_state):
    import re

    line_match = re.match(r'^>([^ ]+) ([^ ]+) (\d+) \d+\s*$', line_from_delta_file)
    if line_match:
        current_contig = line_match.group(1)
        external_genome.set_current_contig(line_match.group(2))
        parser_state['contig_sizes'][current_contig] = int(line_match.group(3))
        franken_genome.add_contig(current_contig)
    else:
        line_match = re.match(r'^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$', line_from_delta_file)
        if line_match:
            parser_state['reference_pos'] = int(line_match.group(1))
            parser_state['final_pos'] = int(line_match.group(2))
            parser_state['external_pos'] = int(line_match.group(3))
            parser_state['external_is_reversed'] = (
                True if ( parser_state['external_pos'] > int(line_match.group(4)) ) else False )
        else:
            line_match = re.match(r'^(\-?)(\d+)\s*$', line_from_delta_file)
            if line_match:
                distance_covered = int(line_match.group(2)) - 1
                is_external_insert = ( True if ( line_match.group(1) == '-' ) else False )
                parser_state = _update_genome_from_delta_data(franken_genome, external_genome, parser_state,
                                                              distance_covered, is_external_insert)
    return parser_state


def parse_delta_file(delta_filename, franken_genome, external_genome):
    parser_state = dict(zip(['contig_sizes', 'reference_pos', 'external_pos', 'final_pos', 'external_is_reversed'],
                            [dict(), None, None, None, None]))
    delta_handle = open(delta_filename, 'r')
    for line_from_delta_file in delta_handle:
        parser_state = _parse_delta_line(line_from_delta_file, franken_genome, external_genome, parser_state)
    delta_handle.close()
    for current_contig in franken_genome.get_contigs():
        franken_genome.extend_contig(parser_state['contig_sizes'][current_contig], 'X', current_contig)


def main():
    import subprocess
    from nasp.nasp_objects import Genome, GenomeMeta

    commandline_args = _parse_args()
    external_nickname = commandline_args.name if commandline_args.name else GenomeMeta.generate_nickname_from_filename(
        commandline_args.external)
    #external_genome = Genome()
    #external_genome.import_fasta_file(commandline_args.external)
    generate_delta_file(commandline_args.nucmerpath, commandline_args.nucmerargs, commandline_args.deltafilterpath,
                        commandline_args.deltafilterargs, external_nickname, commandline_args.reference,
                        commandline_args.external)
    #franken_genome = Genome()
    #parse_delta_file(( external_nickname + ".filtered.delta" ), franken_genome, external_genome)
    #franken_genome.write_to_fasta_file(external_nickname + ".frankenfasta", "franken::")

    import pkg_resources
    nasptool_path = pkg_resources.resource_filename('nasp', 'nasptool_linux_64')
    with open(external_nickname + ".frankenfasta", 'w') as handle:
        return_code = subprocess.call([nasptool_path, "frankenfasta", external_nickname + ".filtered.delta"], stdout=handle)


if __name__ == "__main__":
    main()


