#!/usr/bin/env python3

__author__ = 'jtravis'

import itertools
from concurrent.futures import ProcessPoolExecutor
from nasp2.matrix_DTO import parse_dto
from nasp2.parse import Vcf, Fasta
from nasp2.analyze import analyze_samples


# def explain(matrix_parameters, sample_groups):
#     from nasp2.analyze import sample_positions, analyze_position
#
#     while True:
#         try:
#             print('\a')
#             contig_name, _, position = input("Enter LocusID: ").partition('::')
#             index = int(position) - 1
#         except ValueError:
#             print('LocusID is <contig name>::<position number> such as 500WT1::42')
#             continue
#         reference_contig = matrix_parameters.reference_fasta.get_contig(contig_name)
#         dups_contig = matrix_parameters.reference_dups.get_contig(contig_name)
#         print("Contig Objects:")
#         print(reference_contig)
#         print(dups_contig)
#         for sample in sample_groups:
#             print(sample[0].name)
#             for analysis in sample:
#                 print('\t', analysis.get_contig(reference_contig.name))
#
#         print("Scanning files...")
#         for idx, row in enumerate(zip(reference_contig.positions, dups_contig.positions, sample_positions(reference_contig.name, sample_groups))):
#             if index == idx:
#                 print('\a')
#                 #print('Positions:', "\n".join(str(row)), '\n')
#                 print('Position Analysis:')
#                 print(analyze_position(row[0], row[1], row[2]))
#                 break


def _parse_args():
    """
    Parse command line arguments and return an object whose attributes are the
    settings provided, as per standard argparse behavior.
    For backward compatibility when we were changing dispatchers, the options
    could be directly provided on the command line, or a single argument could
    specify the path to the XML file with the remaining arguments in it.  As
    the old dispatcher is deprecated, so perhaps is every option below except
    '--dto-file', which would then be made required.

    Returns:
        argparse.Namespace: The properties of this object are the parsed command line arguments
    """
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    # parser.add_argument("--mode", required=True, choices=['commandline', 'xml'],
    #                     help="Data passing mode, must be set to 'commandline' or 'xml'.")
    parser.add_argument("--mode", choices=['commandline', 'xml'],
                        help="Data passing mode, must be set to 'commandline' or 'xml'.")
    parser.add_argument("--reference-fasta", help="Path to input reference fasta file.")
    parser.add_argument("--reference-dups", help="Path to input reference dups file.")
    parser.add_argument("--input-files", nargs="+", help="Path to input VCF/fasta files for matrix conversion.")
    parser.add_argument("--matrix-folder", default="matrices", help="Name of folder to write output matrices to.")
    parser.add_argument("--stats-folder", default="statistics", help="Name of folder to write statistics files to.")
    parser.add_argument("--minimum-coverage", type=int, default=10, help="Minimum coverage depth at a position.")
    parser.add_argument("--minimum-proportion", type=float, default=0.9,
                        help="Minimum proportion of reads that must match the call at a position.")
    # parser.add_argument("--num-threads", type=int, default=1, help="Number of threads to use when processing input.")
    # A default of None here means the max number of worker processes will be the number of CPUs on the machine.
    parser.add_argument("--num-threads", type=int, default=None, help="Number of threads to use when processing input.")
    parser.add_argument("--dto-file", required=True, help="Path to a matrix_dto XML file that defines all the parameters.")
    return parser.parse_args()


def main():
    # TODO: remove
    OUTPUT_DIR = './'

    arguments = _parse_args()
    matrix_params = parse_dto(arguments.dto_file)

    matrix_params['reference_fasta'] = matrix_params.get('reference_fasta', arguments.reference_fasta)
    matrix_params['reference_dups'] = matrix_params.get('reference_dups', arguments.reference_dups)
    matrix_params['matrix_folder'] = matrix_params.get('matrix_folder', arguments.matrix_folder)
    matrix_params['stats_folder'] = matrix_params.get('stats_folder', arguments.stats_folder)
    matrix_params['minimum_coverage'] = float(matrix_params.get('minimum_coverage', arguments.minimum_coverage))
    matrix_params['minimum_proportion'] = float(matrix_params.get('minimum_proportion', arguments.minimum_proportion))

    futures = []
    with ProcessPoolExecutor() as executor:
        print("Building contig indices...")

        reference_fasta = executor.submit(Fasta, matrix_params['reference_fasta'], 'reference', None, True)
        # TODO: handle undefined reference_dups
        reference_dups = executor.submit(Fasta, matrix_params['reference_dups'], 'reference', None)

        # Index Vcf and Frankenfastas in parallel.
        for frankenfasta in matrix_params['frankenfasta']:
            futures.append(executor.submit(Fasta, frankenfasta.path, frankenfasta.name, frankenfasta.aligner))
        for vcf in matrix_params['vcf']:
            futures.append(executor.submit(Vcf, vcf.path, vcf.name, vcf.aligner, vcf.snpcaller))

        # Return when the indexing processes complete.
        # TODO: try/catch futures exception to return a failed genome object and write the error to the parse log.
        # Group the analyses by sample name in order to collect sample-level statistics.
        # The SampleAnalyses are sorted before grouping because groups are determined by when the key changes.
        # See analyse.sample_positions() for more details regarding the structure of sample_groups
        sample_groups = tuple(tuple(v) for _, v in itertools.groupby(sorted(future.result() for future in futures), lambda x: x.name))
        # return MatrixParameters(**matrix_params), sample_analyses

        # if len(sys.argv) > 2 and sys.argv[2] == "explain":
        #     return explain(matrix_parameters, sample_groups)

        print("Starting analysis...")
        analyze_samples(OUTPUT_DIR, reference_fasta.result(), reference_dups.result(), sample_groups, arguments.num_threads)


if __name__ == '__main__':
    main()
