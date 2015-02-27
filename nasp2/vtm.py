#!/usr/bin/env python3

__author__ = 'jtravis'


def main():
    from nasp2.matrix_DTO import parse_dto
    matrix_parameters, sample_analyses = parse_dto("/Users/jtravis/NASP/results/matrix_dto.xml")
    # matrix_parameters, sample_analyses = parse_dto("/Users/jtravis/test_time/matrix_dto.xml")

    from nasp2.analyze import analyze_samples

    analyze_samples(matrix_parameters.reference_fasta, matrix_parameters.reference_dups, sample_analyses)


if __name__ == '__main__':
    main()