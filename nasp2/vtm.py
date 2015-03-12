#!/usr/bin/env python3

__author__ = 'jtravis'

import sys

import itertools
import sys
from nasp2.matrix_DTO import parse_dto


def explain(matrix_parameters, sample_analyses):
    sample_groups = tuple(tuple(v) for _, v in itertools.groupby(sorted(sample_analyses), lambda x: x.name))

    from nasp2.analyze import sample_positions, analyze_position

    while True:
        try:
            print('\a')
            contig_name, _, position = input("Enter LocusID: ").partition('::')
            index = int(position) - 1
        except ValueError:
            print('LocusID is <contig name>::<position number> such as 500WT1::42')
            continue
        reference_contig = matrix_parameters.reference_fasta.get_contig(contig_name)
        dups_contig = matrix_parameters.reference_dups.get_contig(contig_name)
        for idx, row in enumerate(zip(reference_contig.positions, dups_contig.positions, sample_positions(reference_contig.name, sample_groups))):
            if index == idx:
                print('\a')
                print('Positions: ', row, '\n')
                print('Position Analysis:')
                print(analyze_position(row[0], row[1], row[2]))
                break


def main():
    from nasp2.matrix_DTO import parse_dto
    print("Building contig indices...")
    matrix_parameters, sample_analyses = parse_dto(sys.argv[1])

    if len(sys.argv) > 2 and sys.argv[2] == "explain":
        return explain(matrix_parameters, sample_analyses)

    from nasp2.analyze import analyze_samples

    analyze_samples(matrix_parameters.reference_fasta, matrix_parameters.reference_dups, sample_analyses)


if __name__ == '__main__':
    main()
