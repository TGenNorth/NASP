#!/usr/bin/env python3

__author__ = 'jtravis'

import itertools
from concurrent.futures import ProcessPoolExecutor

from nasp.vtm.matrix_DTO import parse_dto
from nasp.vtm.parse import Vcf, Fasta
from nasp.vtm.analyze import GenomeAnalysis
from nasp.vtm.write_matrix import analyze_samples


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
    # help="Data passing mode, must be set to 'commandline' or 'xml'.")
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
    parser.add_argument("--dto-file", required=True,
                        help="Path to a matrix_dto XML file that defines all the parameters.")
    return parser.parse_args()


def _index_contigs(reference_fasta, reference_dups, frankenfastas, vcfs, num_workers=None):
    """
    At initialization, each SampleAnalysis will scan their file to build an index of contig file positions.
    This allows contigs be processed in parallel without assuming they appear in any order or exist in any given
    sample. _index_contigs is a helper function that uses a ProcessPoolExecutor to allow this initialization to run in
    parallel.

    Args:
        reference_fasta (str): Path to the reference fasta
        reference_dups (str):
        frankenfastas (tuple of NaspFile): Normalized external Fasta sample analyses
        vcfs (tuple of NaspFile): VCF sample analyses
        num_workers (int or None): The max number of process pool workers or None to default to the number of available
        processors.

    Return:
        Fasta, Fasta, tuple of SampleAnalysis: The given files wrapped in objects that will uniformly parse the files
        using the interface defined by SampleAnalysis. The tuple of SampleAnalysis are sorted lexically by identifier.
    """
    futures = []
    with ProcessPoolExecutor(num_workers) as executor:
        reference_fasta = executor.submit(Fasta, reference_fasta, 'reference', None, True)
        # TODO: handle undefined reference_dups
        reference_dups = executor.submit(Fasta, reference_dups, 'reference', None)

        # FIXME: frankenfasta and vcfs are not read from the commandline
        # Index Vcf and Frankenfastas in parallel.
        for frankenfasta in frankenfastas:
            futures.append(executor.submit(Fasta, frankenfasta.path, frankenfasta.name, frankenfasta.aligner))
        for vcf in vcfs:
            futures.append(executor.submit(Vcf, vcf.path, vcf.name, vcf.aligner, vcf.snpcaller))

        # Return when the indexing processes complete.
        # TODO: try/catch futures exception to return a failed genome object and write the error to the parse log.
        # Group the analyses by sample name in order to collect sample-level statistics.
        # The SampleAnalyses are sorted before grouping because groups are determined by when the key changes.
        # See analyse.sample_positions() for more details regarding the structure of sample_groups
        sample_groups = tuple(
            tuple(v) for _, v in itertools.groupby(sorted(future.result() for future in futures), lambda x: x.name)
        )

        executor.shutdown()

    return reference_fasta.result(), reference_dups.result(), sample_groups


def main():
    arguments = _parse_args()
    matrix_params = parse_dto(arguments.dto_file)

    reference_fasta = matrix_params.get('reference_fasta', arguments.reference_fasta)
    reference_dups = matrix_params.get('reference_dups', arguments.reference_dups)
    matrix_dir = matrix_params.get('matrix_folder', arguments.matrix_folder)
    stats_dir = matrix_params.get('stats_folder', arguments.stats_folder)
    # TODO: handle cast to float exception
    coverage_threshold = float(matrix_params.get('minimum_coverage', arguments.minimum_coverage))
    proportion_threshold = float(matrix_params.get('minimum_proportion', arguments.minimum_proportion))

    print("Building contig indices...")
    (reference, dups, sample_groups) = _index_contigs(reference_fasta, reference_dups, matrix_params['frankenfasta'],
                                                      matrix_params['vcf'], arguments.num_threads)

    print("Starting analysis...")
    genome_analysis = GenomeAnalysis(coverage_threshold, proportion_threshold)
    analyze_samples(matrix_dir, stats_dir, genome_analysis, reference, dups, sample_groups, arguments.num_threads)


if __name__ == '__main__':
    main()
