"""
write_matrix handles all disk IO operations.

TODO: Most of the functions are written as coroutines that could be used improve performance by using asyncio to
continue processing instead of waiting for data to be read/written to/from disk.
"""
__author__ = 'jtravis'
__version__ = "0.9.8"

from nasp import __version__ as nasp_version

import os
import csv
from collections import Counter
from contextlib import ExitStack
from tempfile import TemporaryDirectory
import itertools
import functools


def get_vcf_metadata(nasp_version, identifiers, contigs):
    """
    Args:
        nasp_version (str): The current nasp release version.
        identifiers (tuple): contig_name::aligner,snpcaller
        contigs (tuple of Contig):

    Return:
        str: VCF file metadata
    """
    vcf_metadata = "##fileFormat=VCFv4.2\n##source=NASPv{0}\n".format(nasp_version)
    vcf_metadata += "\n".join(
        "##contig=<ID=\"{0}\",length={1}>".format(contig.name, len(contig)) for contig in contigs)
    vcf_metadata += "\n".join(
        "##SAMPLE=<ID=\"{0}\",Genomes=\"{0}\",Mixture=1.0>".format(identifier) for identifier in identifiers)
    vcf_metadata += ("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
                     "##FILTER=<ID=NoCall,Description=\"No call for this sample at this position\">\n"
                     "##FILTER=<ID=CovFail,Description=\"Insufficient depth of coverage for this sample at this position\">\n"
                     "##FILTER=<ID=PropFail,Description=\"Insufficient proportion of reads were variant for this sample at this position\">\n"
                     "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                     "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filters that failed for this sample at this position\">\n")
    return vcf_metadata


def get_header(type, identifiers):
    """
    Args:
        type (str): 'master', 'missingdata', 'best_snp', 'vcf'
        identifiers (tuple): A list of sample_name::aligner,snpcaller or just sample_name if type is 'best_snp'

    Returns:
        list: Header columns for the requested file type.
    """
    if type == 'vcf':
        return ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT') + identifiers

    matrix_header = ('LocusID', 'Reference') + identifiers

    # TODO: bestsnp consensus
    if type == 'best_snp':
        # header += consensus column/sample_name
        pass

    matrix_header += (
        '#SNPcall', '#Indelcall', '#Refcall', '#CallWasMade', '#PassedDepthFilter', '#PassedProportionFilter', '#A',
        '#C', '#G', '#T', '#Indel', '#NXdegen', 'Contig', 'Position', 'InDupRegion', 'SampleConsensus'
    )

    # all callable or missing data
    if type in ['master', 'missingdata']:
        matrix_header += ('CallWasMade', 'PassedDepthFilter', 'PassedProportionFilter')

    matrix_header += ('Pattern', 'Pattern#')

    return matrix_header


def _vcf_filter_column(coverage_threshold, proportion_threshold, is_all_pass_coverage, is_all_pass_proportion):
    """
    Return:
        str: 'PASS' if all filters passed, or an error code.

    Example:
        'PASS' means all filters passed.
        'c10' means the coverage was below the user defined threshold of 10 for at least one analysis.
        'p0.9' means the proportion was below the user defined threshold of 0.9 for at least one analysis.
        'c10;p0.9' means at least one analysis was below either the coverage or proportion threshold.
    """
    filter = []
    if not is_all_pass_coverage:
        filter.append('c{0}'.format(coverage_threshold))
    if not is_all_pass_proportion:
        filter.append('p{0}'.format(proportion_threshold))

    if filter:
        return ';'.join(filter)

    return 'PASS'


def _vcf_analysis_column(pattern, analysis_stats):
    """
    Args:
        pattern (str):
        analysis_stats:

    Return:
        str: The format string for each analysis column.
    """
    for pattern_num, analysis_stat in zip(pattern, itertools.chain.from_iterable(analysis_stats)):
        try:
            gt = int(pattern_num) - 1
        except ValueError:
            # pattern_num for the analysis was 'N'
            gt = '.'
        if not analysis_stat['was_called']:
            ft = "NoCall"
        elif not analysis_stat['passed_coverage_filter']:
            ft = "CovFail"
        elif not analysis_stat['passed_proportion_filter']:
            ft = "PropFail"
        else:
            ft = "PASS"
        yield '{0}:{1}'.format(gt, ft)


def write_missingdata_vcf(directory, contig_name, identifiers, metadata):
    """
    Args:
        directory (str):
        contig_name (str):
        identifiers:
    """
    coverage_threshold = 10
    proportion_threshold = 0.9

    with open('{0}_missingdata.vcf'.format(os.path.join(directory, contig_name)), 'w') as handle:
        handle.write(metadata)
        writer = csv.DictWriter(handle, fieldnames=get_header('vcf', identifiers), delimiter='\t')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_missing_matrix:
                continue

            ref = row.call_str[0]
            alts = set(row.call_str[1:])
            alts.difference_update(('X', 'N', ref))

            line = {
                '#CHROM': contig_name,
                'POS': position,
                'ID': '.',
                'REF': ref,
                'ALT': ','.join(alts) or '.',
                'QUAL': '.',
                'FILTER': _vcf_filter_column(coverage_threshold, proportion_threshold, row.is_all_passed_consensus,
                                             row.is_all_passed_proportion),
                # TODO: Add #indel stat to NS
                'INFO': 'AN={0};NS={1}'.format(len(alts) + 1, row.called_reference + row.called_snp),
                'FORMAT': 'GT:FT'
            }
            # Match each analysis with its analysis column.
            line.update({k: v for k, v in zip(identifiers, _vcf_analysis_column(row.Pattern, row.all_sample_stats))})
            writer.writerow(line)


def write_bestsnp_vcf(directory, contig_name, identifiers, metadata):
    """
    Args:
        directory (str):
        contig_name (str):
        identifiers:
    """
    coverage_threshold = 10
    proportion_threshold = 0.9

    with open('{0}_bestsnp.vcf'.format(os.path.join(directory, contig_name)), 'w') as handle:
        handle.write(metadata)
        # writer = csv.DictWriter(handle, fieldnames=get_header('vcf', identifiers), delimiter='\t')
        # writer.writeheader()
        # position = 0
        # while True:
        #     row = yield
        #     position += 1
        #
        #     if not row.is_best_snp:
        #         continue
        #
        #     ref = row.call_str[0]
        #     alts = set(row.call_str[1:])
        #     alts.discard(ref)
        #
        #     line = {
        #         '#CHROM': contig_name,
        #         'POS': position,
        #         'ID': '.',
        #         'REF': ref,
        #         'ALT': ','.join(alts) or '.',
        #         'QUAL': '.',
        #         'FILTER': _vcf_filter_column(coverage_threshold, proportion_threshold, row.is_all_passed_consensus,
        #                                      row.is_all_passed_proportion),
        #         # TODO: AN is the number of snps + 1 for the reference.
        #         # TODO: Add #indel stat to NS
        #         'INFO': 'AN={0};NS={1}'.format(len(alts) + 1, row.called_reference + row.called_snp),
        #         'FORMAT': 'GT:FT'
        #     }
        #     # Match each analysis with its analysis column.
        #     line.update({k: v for k, v in zip(identifiers, _vcf_analysis_column(row.Pattern, row.all_sample_stats))})
        #     writer.writerow(line)

        handle.write('{0}\n'.format('\t'.join(get_header('vcf', identifiers))))

        position = 0

        while True:
            row = yield

            if not row.is_best_snp:
                continue

            position += 1

            ref = row.call_str[0]
            alts = set(row.call_str[1:])
            alts.discard(ref)

            handle.write(
                '{0}\t{1}\t.\t{2}\t{3}\t.\t{4}\tGT:FT\n'.format(
                    contig_name,
                    position,
                    ref,
                    ','.join(alts) or '.',
                    '\t'.join(_vcf_analysis_column(row.Pattern, row.all_sample_stats))
                )
            )

def write_sample_stats(filepath, sample_stats, sample_groups, reference_length):
    """
    Args:
        filepath (str):
        sample_stats (list of lists of Counters):
        sample_groups (tuple of tuples of SampleAnalysis):
        reference_length (int): Total number of positions in the reference.

    Example:
        sample_stats = (
            # Any/All Counter summary for the whole genome.
            ({
                 'was_called': 0,
                 'passed_coverage_filter': 0,
                 'passed_proportion_filter': 0,
                 'quality_breadth': 0,
                 'called_reference': 0,
                 'called_snp': 0,
                 'called_degen': 0
             }, {
                 'was_called': 0,
                 'passed_coverage_filter': 0,
                 'passed_proportion_filter': 0,
                 'quality_breadth': 0,
                 'called_reference': 0,
                 'called_snp': 0,
                 'called_degen': 0
             }),
             # Followed by a tuple for each sample
             # Where the first two elements are a Any/All Counter summary for the sample
             # Followed by a Counter for each analysis combination
             ({
                 'was_called': 0,
                 'passed_coverage_filter': 0,
                 'passed_proportion_filter': 0,
                 'quality_breadth': 0,
                 'called_reference': 0,
                 'called_snp': 0,
                 'called_degen': 0
              },
              ...
             ),
        )
    """
    fieldnames = ('Sample', 'Sample::Analysis', 'was_called', 'was_called (%)', 'passed_coverage_filter',
                  'passed_coverage_filter (%)', 'passed_proportion_filter', 'passed_proportion_filter (%)',
                  'quality_breadth', 'quality_breadth (%)', 'called_reference', 'called_reference (%)',
                  'called_snp', 'called_snp (%)', 'called_degen', 'called_degen (%)')

    genome_any = sample_stats[0][0]
    genome_all = sample_stats[0][1]

    # Join identifiers with genome any/all summary data.
    genome_any['Sample'] = '[any]'
    genome_all['Sample'] = '[all]'

    # Calculate genome any/all stat percentages.
    for stat in fieldnames[2::2]:
        genome_any[stat + ' (%)'] = "{0:.2f}%".format(genome_any[stat] / reference_length * 100)
        genome_all[stat + ' (%)'] = "{0:.2f}%".format(genome_all[stat] / reference_length * 100)

    with open(filepath, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        handle.write('\n')
        # Write the whole genome any/all summaries.
        writer.writerow(genome_any)
        writer.writerow(genome_all)

        # Join sample names and analysis identifiers with the analysis stats before writing them to file.
        # sample_stats[1:] skips the overall genome summaries which have no corresponding sample.
        for sample_stat, sample in zip(sample_stats[1:], sample_groups):
            sample_any = sample_stat[0]
            sample_all = sample_stat[1]

            handle.write('\n')
            # Join identifiers with per-analysis any/all summary data.
            sample_any['Sample'] = sample[0].name
            sample_any['Sample::Analysis'] = '[any]'
            sample_all['Sample'] = sample[0].name
            sample_all['Sample::Analysis'] = '[all]'

            # Calculate sample any/all summary stat percentages
            for stat in fieldnames[2::2]:
                sample_any[stat + ' (%)'] = "{0:.2f}%".format(sample_any[stat] / reference_length * 100)
                sample_all[stat + ' (%)'] = "{0:.2f}%".format(sample_all[stat] / reference_length * 100)
            writer.writerow(sample_any)
            writer.writerow(sample_all)
            # Join sample identifiers with sample data then write to file.
            for analysis_stats, identifier in zip(sample_stat[2:], sample):
                # Calculate per-analysis stat percentages
                for stat in fieldnames[2::2]:
                    analysis_stats[stat + ' (%)'] = "{0:.2f}%".format(analysis_stats[stat] / reference_length * 100)
                analysis_stats['Sample'] = identifier.name
                analysis_stats['Sample::Analysis'] = identifier.identifier
                writer.writerow(analysis_stats)


def write_general_stats(filepath, contig_stats):
    """
    Args:
        filepath (str):
        contig_stats (tuple of Counter):

    Example:
        contig_stats = (Counter({
            'Contig': 'contig_name',
            'reference_length': 0,
            'reference_clean': 0,
            'reference_duplicated': 0,
            'all_called': 0,
            'all_passed_coverage': 0,
            'all_passed_proportion': 0,
            'all_passed_consensus': 0,
            'quality_breadth': 0,
            'any_snps': 0,
            'best_snps': 0
        }),)

    Return:
        int: Total reference length.
    """
    # Sum contig stats
    whole_genome_stats = Counter({'Contig': ''})
    for contig_stat in contig_stats:
        whole_genome_stats.update(contig_stat)
    whole_genome_stats['Contig'] = 'Whole Genome'

    # NOTE: Raises ZeroDivisionError if reference_length is 0.
    reference_length = whole_genome_stats['reference_length']

    with open(filepath, 'w') as handle:
        fieldnames = ('Contig', 'reference_length', 'reference_clean', 'reference_clean (%)', 'reference_duplicated',
                      'reference_duplicated (%)', 'all_called', 'all_called (%)', 'all_passed_coverage',
                      'all_passed_coverage (%)', 'all_passed_proportion', 'all_passed_proportion (%)',
                      'all_passed_consensus', 'all_passed_consensus (%)', 'quality_breadth', 'quality_breadth (%)',
                      'any_snps', 'any_snps (%)', 'best_snps', 'best_snps (%)')
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        handle.write('\n')
        # Calculate whole genome stat percentages.
        for stat in fieldnames[2::2]:
            whole_genome_stats[stat + ' (%)'] = "{0:.2f}%".format(whole_genome_stats[stat] / reference_length * 100)
        writer.writerow(whole_genome_stats)
        for contig_stat in contig_stats:
            # Calculate contig stat percentages.
            for stat in fieldnames[2::2]:
                contig_stat[stat + ' (%)'] = "{0:.2f}%".format(
                    contig_stat[stat] / contig_stat['reference_length'] * 100)
            writer.writerow(contig_stat)

    return reference_length


def write_master_matrix(directory, contig_name, identifiers):
    """
    Args:
        filepath (str): Path to the output file.
        contig_name (str): Name
        identifiers (tuple of
    """
    with open('{0}_master.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # writer = csv.DictWriter(handle, fieldnames=get_header('master', identifiers), delimiter='\t', lineterminator='\n')
        # print(get_header('master', identifiers))
        # writer.writeheader()
        # position = 0
        # while True:
        # row = yield
        # position += 1
        #     # num_samples is the number of analyses not including the reference.
        #     num_samples = len(row.call_str) - 1
        #
        #     line = {
        #         'LocusID': "{0}::{1}".format(contig_name, position),
        #         'Reference': row.call_str[0],
        #         '#SNPcall': row.called_snp,
        #         # TODO: replace with n/a
        #         '#Indelcall': '0',
        #         '#Refcall': row.called_reference,
        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
        #         '#A': row.num_A,
        #         '#C': row.num_C,
        #         '#G': row.num_G,
        #         '#T': row.num_T,
        #         # TODO: replace with n/a
        #         '#Indel': '0',
        #         '#NXdegen': row.num_N,
        #         'Contig': contig_name,
        #         'Position': position,
        #         'InDupRegion': row.is_reference_duplicated,
        #         'SampleConsensus': row.is_all_passed_consensus,
        #         'CallWasMade': row.CallWasMade,
        #         'PassedDepthFilter': row.PassedDepthFilter,
        #         'PassedProportionFilter': row.PassedProportionFilter,
        #         'Pattern': "".join(row.Pattern)
        #     }
        #     # Match each base call with its sample analysis column.
        #     line.update({k: v for k, v in zip(identifiers, row.call_str[1:])})
        #
        #     writer.writerow(line)
        handle.write('{0}\n'.format('\t'.join(get_header('master', identifiers))))

        position = 0
        num_samples = len(identifiers)

        while True:
            row = yield

            position += 1

            handle.write(
                '{0}::{1}\t{2}\t{3}\t0\t{4}\t{5:d}/{8:d}\t{6:d}/{8:d}\t{7:d}/{8:d}\t{9}\t{10}\t{11}\t{12}\t0\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\n'.format(
                    contig_name,
                    position,
                    '\t'.join(row.call_str),
                    row.called_snp,
                    # FIXME: #Indelcall unimplemented
                    row.called_reference,
                    num_samples - row.CallWasMade.count('N'),
                    row.passed_coverage_filter,
                    row.passed_proportion_filter,
                    num_samples,
                    row.num_A,
                    row.num_C,
                    row.num_G,
                    row.num_T,
                    # FIXME: #Index unimplemented
                    row.num_N,
                    contig_name,
                    position,
                    row.is_reference_duplicated,
                    row.is_all_passed_consensus,
                    row.CallWasMade,
                    row.PassedDepthFilter,
                    row.PassedProportionFilter,
                    "".join(row.Pattern)
                )
            )


def write_missingdata_matrix(directory, contig_name, identifiers):
    """
    Args:
        directory(str):
        contig_name (str):
        identifiers:
    """
    with open('{0}_missingdata.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # writer = csv.DictWriter(handle, fieldnames=get_header('missingdata', identifiers), delimiter='\t',
        #                         lineterminator='\n')
        # writer.writeheader()
        # position = 0
        # while True:
        #     row = yield
        #     position += 1
        #
        #     if not row.is_missing_matrix:
        #         continue
        #
        #     # num_samples is the number of analyses not including the reference.
        #     num_samples = len(row.masked_call_str) - 1
        #
        #     line = {
        #         'LocusID': "{0}::{1}".format(contig_name, position),
        #         'Reference': row.masked_call_str[0],
        #         '#SNPcall': row.called_snp,
        #         # TODO: replace with n/a
        #         '#Indelcall': '0',
        #         '#Refcall': row.called_reference,
        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
        #         '#A': row.num_A,
        #         '#C': row.num_C,
        #         '#G': row.num_G,
        #         '#T': row.num_T,
        #         # TODO: replace with n/a
        #         '#Indel': '0',
        #         '#NXdegen': row.num_N,
        #         'Contig': contig_name,
        #         'Position': position,
        #         'InDupRegion': row.is_reference_duplicated,
        #         'SampleConsensus': row.is_all_passed_consensus,
        #         'CallWasMade': row.CallWasMade,
        #         'PassedDepthFilter': row.PassedDepthFilter,
        #         'PassedProportionFilter': row.PassedProportionFilter,
        #         'Pattern': "".join(row.Pattern)
        #     }
        #     # Match each base call with its sample analysis column.
        #     line.update({k: v for k, v in zip(identifiers, row.masked_call_str[1:])})
        #     writer.writerow(line)

        handle.write('{0}\n'.join(get_header('missingdata', identifiers)))

        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_missing_matrix:
                continue

            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.masked_call_str) - 1

            handle.write(
                '{0}::{1}\t{2}\t{3}\t0\t{4}\t{5:d}/{8:d}\t{6:d}/{8:d}\t{7:d}/{8:d}\t{9}\t{10}\t{11}\t{12}\t0\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\n'.format(
                    contig_name,
                    position,
                    '\t'.join(row.call_str),
                    row.called_snp,
                    # FIXME: #Indelcall unimplemented
                    row.called_reference,
                    num_samples - row.CallWasMade.count('N'),
                    row.passed_coverage_filter,
                    row.passed_proportion_filter,
                    num_samples,
                    row.num_A,
                    row.num_C,
                    row.num_G,
                    row.num_T,
                    # FIXME: #Index unimplemented
                    row.num_N,
                    contig_name,
                    position,
                    row.is_reference_duplicated,
                    row.is_all_passed_consensus,
                    row.CallWasMade,
                    row.PassedDepthFilter,
                    row.PassedProportionFilter,
                    "".join(row.Pattern)
                )
            )


def write_bestsnp_matrix(directory, contig_name, sample_groups):
    """
    Args:
        directory (str):
        contig_name (str):
        sample_groups
    """
    sample_names = tuple(sample[0].name for sample in sample_groups)
    # first_analysis_index is a list of the index of the first analysis for each sample in the call string
    first_analysis_index = []
    # num_analyses starts at 1 to skip the reference call
    num_analyses = 1
    for sample in sample_groups:
        first_analysis_index.append(num_analyses)
        num_analyses += len(sample)

    with open('{0}_bestsnp.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # writer = csv.DictWriter(handle, fieldnames=get_header('best_snp', sample_names), delimiter='\t',
        #                         lineterminator='\n')
        # writer.writeheader()
        # position = 0
        # while True:
        #     row = yield
        #     position += 1
        #
        #     if not row.is_best_snp:
        #         continue
        #
        #     # num_samples is the number of analyses not including the reference.
        #     num_samples = len(row.call_str) - 1
        #
        #     line = {
        #         'LocusID': "{0}::{1}".format(contig_name, position),
        #         'Reference': row.call_str[0],
        #         '#SNPcall': row.called_snp,
        #         # TODO: replace with n/a
        #         '#Indelcall': '0',
        #         '#Refcall': row.called_reference,
        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
        #         '#A': row.num_A,
        #         '#C': row.num_C,
        #         '#G': row.num_G,
        #         '#T': row.num_T,
        #         # TODO: replace with n/a
        #         '#Indel': '0',
        #         '#NXdegen': row.num_N,
        #         'Contig': contig_name,
        #         'Position': position,
        #         'InDupRegion': row.is_reference_duplicated,
        #         'SampleConsensus': row.is_all_passed_consensus,
        #         'Pattern': "".join(row.Pattern)
        #     }
        #     # Match each base call with its sample analysis column.
        #     line.update(
        #         {sample_name: row.call_str[index] for sample_name, index in zip(sample_names, first_analysis_index)})
        #
        #     writer.writerow(line)

        handle.write('{0}\n'.join(get_header('best_snp', sample_names)))

        position = 0
        num_samples = num_analyses - 1
        while True:
            row = yield
            position += 1

            if not row.is_best_snp:
                continue

            handle.write(
                '{0}::{1}\t{2}\t{3}\t0\t{4}\t{5:d}/{8:d}\t{6:d}/{8:d}\t{7:d}/{8:d}\t{9}\t{10}\t{11}\t{12}\t0\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\n'.format(
                    contig_name,
                    position,
                    '\t'.join(row.call_str),
                    row.called_snp,
                    # FIXME: #Indelcall unimplemented
                    row.called_reference,
                    num_samples - row.CallWasMade.count('N'),
                    row.passed_coverage_filter,
                    row.passed_proportion_filter,
                    num_samples,
                    row.num_A,
                    row.num_C,
                    row.num_G,
                    row.num_T,
                    # FIXME: #Index unimplemented
                    row.num_N,
                    contig_name,
                    position,
                    row.is_reference_duplicated,
                    row.is_all_passed_consensus,
                    "".join(row.Pattern)
                )
            )


def write_withallrefpos_matrix(directory, contig_name, identifiers):
    """
    includeref_matrix is identical to the master matrix except it uses the same call masking as the missingdata matrix
    for low quality positions.

    Args:
        directory (str):
        contig_name (str):
        identifiers (tuple):
    """
    with open('{0}_withallrefpos.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # writer = csv.DictWriter(handle, fieldnames=get_header('best_snp', identifiers), delimiter='\t',
        #                         lineterminator='\n')
        # writer.writeheader()
        # position = 0
        # while True:
        #     row = yield
        #     position += 1
        #
        #     # if not row.is_all_quality_breadth:
        #     # continue
        #
        #     # num_samples is the number of analyses not including the reference.
        #     num_samples = len(row.call_str) - 1
        #
        #     line = {
        #         'LocusID': "{0}::{1}".format(contig_name, position),
        #         'Reference': row.call_str[0],
        #         '#SNPcall': row.called_snp,
        #         # TODO: replace with n/a
        #         '#Indelcall': '0',
        #         '#Refcall': row.called_reference,
        #         '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
        #         '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
        #         '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
        #         '#A': row.num_A,
        #         '#C': row.num_C,
        #         '#G': row.num_G,
        #         '#T': row.num_T,
        #         # TODO: replace with n/a
        #         '#Indel': '0',
        #         '#NXdegen': row.num_N,
        #         'Contig': contig_name,
        #         'Position': position,
        #         'InDupRegion': row.is_reference_duplicated,
        #         'SampleConsensus': row.is_all_passed_consensus,
        #         'Pattern': "".join(row.Pattern)
        #     }
        #     # Match each base call with its sample analysis column.
        #     line.update({k: v for k, v in zip(identifiers, row.call_str[1:])})
        #
        #     writer.writerow(line)

        handle.write('{0}\n'.format('\t'.join(get_header('best_snp', identifiers))))

        position = 0
        num_samples = len(identifiers)

        while True:
            row = yield

            position += 1

            handle.write(
                '{0}::{1}\t{2}\t{3}\t0\t{4}\t{5:d}/{8:d}\t{6:d}/{8:d}\t{7:d}/{8:d}\t{9}\t{10}\t{11}\t{12}\t0\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\n'.format(
                    contig_name,
                    position,
                    '\t'.join(row.masked_call_str),
                    row.called_snp,
                    # FIXME: #Indelcall unimplemented
                    row.called_reference,
                    num_samples - row.CallWasMade.count('N'),
                    row.passed_coverage_filter,
                    row.passed_proportion_filter,
                    num_samples,
                    row.num_A,
                    row.num_C,
                    row.num_G,
                    row.num_T,
                    # FIXME: #Index unimplemented
                    row.num_N,
                    contig_name,
                    position,
                    row.is_reference_duplicated,
                    row.is_all_passed_consensus,
                    row.CallWasMade,
                    row.PassedDepthFilter,
                    row.PassedProportionFilter,
                    "".join(row.Pattern)
                )
            )



def write_bestsnp_snpfasta(directory, contig_name, identifiers):
    """
    Args:
        directory (str):
        contig_name (str):
        identifiers (tuple of Strings):

    Example:
        directory/contig_sample::aligner,snpcaller_bestsnp.fasta
    """
    # All opened files will automatically be closed at the end of
    # the with statement, even if attempts to open files later
    # in the list raise an exception
    with ExitStack() as stack:
        files = tuple(stack.enter_context(
            open('{0}_{1}_bestsnp.fasta'.format(os.path.join(directory, 'fasta_partials', contig_name), identifier),
                 'w')) for identifier in identifiers)

        while True:
            row = yield

            if not row.is_best_snp:
                continue

            for file, call in zip(files, row.call_str):
                file.write(call)


def write_missingdata_snpfasta(directory, contig_name, identifiers):
    """
    Write the calls

    Args:
        directory (str):
        contig_name (str):
        identifiers (tuple of Strings):

    Example:
        directory/contig_sample::aligner,snpcaller_missingdata.fasta
    """
    # All opened files will automatically be closed at the end of
    # the with statement, even if attempts to open files later
    # in the list raise an exception
    with ExitStack() as stack:
        files = tuple(stack.enter_context(
            open('{0}_{1}_missingdata.fasta'.format(os.path.join(directory, 'fasta_partials', contig_name), identifier),
                 'w')) for identifier in identifiers)

        while True:
            row = yield

            if not row.is_missing_matrix:
                continue

            for file, call in zip(files, row.masked_call_str):
                if call != 'X':
                    file.write(call)
                else:
                    file.write('N')


def _concat_matrix(src, dest, offset=0):
    """
    Concat the contig matrix to the final file.

    Args:
        src (str):
        dest (str):
        offset (int): Seek past the metadata
    """
    # logging.info('Started concat {0}'.format(src))
    with open(src) as partial, open(dest, 'a') as complete:
        partial.seek(offset)
        # Discard the header
        partial.readline()
        complete.writelines(partial)
        # logging.info('Completed concat {0}'.format(src))


def _concat_snpfasta_contig(src_dir, contig_name, identifiers, suffix):
    """
    Concat the snpfasta contigs for each analysis into a file.

    Args:
        src_dir (str):
        contig_name (str):
        identifiers (tuple):
        suffix (str):

    Example:
        Concat sampleA contigs:
          src_dir/contigA_sampleA::snpcaller,aligner_bestsnp_matrix.fasta: GATC
        + ...
        + src_dir/contigZ_sampleA::snpcaller,aligner_bestsnp_matrix.fasta: CTAG
        --------------------------------------------------------
        directory/sampleA::snpcaller,aligner_bestsnp_matrix.fasta: GATCCTAG

        Concat sampleB contigs:
          src_dir/contigA_sampleB::snpcaller,aligner_bestsnp_matrix.fasta: ATCG
        + ...
        + src_dir/contigZ_sampleB::snpcaller,aligner_bestsnp_matrix.fasta: TAGC
        --------------------------------------------------------
        src_dir/sampleB::snpcaller,aligner_bestsnp_matrix.fasta: ATCGTAGC

    """
    with ExitStack() as stack:
        analyses = (stack.enter_context(open(os.path.join(src_dir, identifier + suffix), 'a+')) for identifier in
                    identifiers)
        analysis_contigs = (
            stack.enter_context(open(os.path.join(src_dir, 'fasta_partials', contig_name + '_' + identifier + suffix)))
            for
            identifier in identifiers)

        for analysis, analysis_contig in zip(analyses, analysis_contigs):
            analysis.writelines(analysis_contig)


def _concat_snpfasta(dest_dir, src_dir, dest, identifiers, suffix):
    """
    Concat the snpfastas for each analysis into the final file.
    Each analysis is preceded with a description line and the base calls wrap every 80 characters.

    It assumes the partial results are files with a single line of base calls.

    Args:
        dest_dir (str):
        src_dir (str):
        dest (str):
        identifiers (tuple of str):
        suffix (str):

    Example:
        Concat analysis into the final file:
          src_dir/reference_bestsnp_matrix.fasta:GATCGATC
          src_dir/sampleA::snpcaller,aligner_bestsnp_matrix.fasta: GATCCTAG
          ...
        + src_dir/sampleZ::snpcaller,aligner_bestsnp_matrix.fasta: ATCGTAGC
        -------------------------------------------------------------------
        dest_dir/bestsnp_matrix.fasta:
          >Reference
          GATCGATC
          >sampleZ::snpcaller,aligner
          GATCCTAG
          ...
          >sampleZ::snpcaller,aligner.fasta
          ATCGTAGC
    """
    with ExitStack() as stack, open(os.path.join(dest_dir, dest), 'w') as dest:
        for identifier in identifiers:
            analysis = stack.enter_context(open(os.path.join(src_dir, identifier + suffix), 'r'))
            dest.write('>{0}\n'.format(identifier))
            # Wrap lines every 80 characters.
            while True:
                line = analysis.read(80)
                if not line:
                    break
                dest.write('{0}\n'.format(line))


# def _swap_future(executor, task):
# """
# swap_future is a helper to create callback chains of serialized tasks running in parallel.
#     As each future resolves the callback launches the next task in the series and replaces the reference to the
#     completed task with the newly started task so that tasks that follow can be attached in sequence.
#
#     Args:
#         executor (concurrent.futures.Executor):
#         task (callable):
#
#     Return:
#         callable: A function that replaces the completed future with a new future. The argument is a future that
#         resolves to an array and index.
#
#     Example:
#         In addition to concatenating the partial files in parallel:
#         1 |-----------------------> 2 |------------>        3 |-------->
#         1 |------->                 2 |------->             3 |------->
#         1 |------------------->     2 |-------------------> 3 |-------->
#
#         The done callback chain allows the next section to begin as each section completes.
#         1 |-----------------------> 2 |------------> 3 |-------->
#         1 |-------> 2 |-------> 3 |------->
#         1 |-------------------> 2 |-------------------> 3 |-------->
#     """
#     # FIXME: This method of chaining futures is a hack that is probably incorrectly implementing the intended design
#     # pattern of the futures done callback.
#     def callback(future):
#         result = future.result()
#         result[0][result[1]] = executor.submit(task)
#         return result
#     return callback


def _get_write_coroutines(tempdirname, identifiers, sample_groups, vcf_metadata, contig_name):
    """
    Args:
        tempdirname (str):
        identifiers (tuple):
        sample_groups (tuple of SampleAnalysis tuples):
        vcf_metadata (str):
        contig_name (str):
    Return:
        tuple of coroutines: Each coroutine receives Position tuples
    """
    # Exclude the reference identifier
    sample_identifiers = identifiers[1:]

    return (
        write_master_matrix(tempdirname, contig_name, sample_identifiers),
        write_bestsnp_matrix(tempdirname, contig_name, sample_groups),
        write_missingdata_matrix(tempdirname, contig_name, sample_identifiers),
        write_withallrefpos_matrix(tempdirname, contig_name, sample_identifiers),
        write_bestsnp_vcf(tempdirname, contig_name, sample_identifiers, vcf_metadata),
        write_missingdata_vcf(tempdirname, contig_name, sample_identifiers, vcf_metadata),
        write_missingdata_snpfasta(tempdirname, contig_name, identifiers),
        write_bestsnp_snpfasta(tempdirname, contig_name, identifiers)
    )


from multiprocessing import Pool


def analyze_samples(matrix_dir, stats_dir, genome_analysis, reference_fasta, reference_dups, sample_groups,
                    max_workers=None):
    """
    analyze_samples uses a ProcessPool to read contigs across all the sample analyses in parallel,
    write the partial results to a temporary directory, aggregate the results for each final file
    in parallel, and cleanup the temporary directory.

    Note:
        It is assumed the reference has at least one contig.

    Args:
        matrix_dir (str):
        stats_dir (str):
        reference_fasta (Fasta):
        reference_dups (Fasta):
        sample_groups (tuple of tuple of SampleAnalysis): The sample analysis grouped by Sample.
        max_workers (int or None):
    """
    # TODO: Use consistent ordering with what is passed to analyze contig.
    contig_stats = []
    # NOTE: if the analysis loop does not run, None will be passed to write_general_stats.
    # Sample stats
    sample_stats = None
    is_first_contig = True
    # An identifier is a sample_name::aligner,snpcaller header to identify each analysis file.
    # The Reference identifier is used by the snpfasta writers/concat methods, but skipped by the contig analysis.
    identifiers = ('Reference',) + tuple(
        analysis.identifier for analysis in itertools.chain.from_iterable(sample_groups))

    matrices = ('master.tsv', 'bestsnp.tsv', 'missingdata.tsv', 'withallrefpos.tsv', 'bestsnp.vcf', 'missingdata.vcf')
    # futures = [Future()] * len(matrices)

    # Analyze the contigs in parallel. The partial files leading up to the final result will be written in a
    # temporary directory which is deleted automatically.
    # with ProcessPoolExecutor(max_workers=max_workers) as executor, TemporaryDirectory(dir=matrix_dir) as tempdirname:
    with Pool(processes=max_workers) as pool, TemporaryDirectory(dir=matrix_dir) as tempdirname:

        os.makedirs(os.path.join(tempdirname, 'fasta_partials'))

        vcf_metadata = get_vcf_metadata(nasp_version, identifiers, reference_fasta.contigs)
        vcf_metadata_len = len(vcf_metadata)
        coroutine_partial = functools.partial(_get_write_coroutines, tempdirname, identifiers, sample_groups,
                                              vcf_metadata)

        # Only the reference contig is changing, bind the other parameters to the function.
        analyze = functools.partial(genome_analysis.analyze_contig, coroutine_partial, sample_groups, reference_dups)
        for sample_stat, contig_stat in pool.map(analyze, reference_fasta.contigs):
            contig_stats.append(contig_stat)
            contig_name = contig_stat['Contig']

            # Concatenate the contig matrices.
            for matrix in matrices:
                # logging.info('Scheduled concat contig {0}'.format(contig_name))

                # Path to a contig matrix.
                partial = os.path.join(tempdirname, '{0}_{1}'.format(contig_name, matrix))
                # Path to the final matrix where all the contigs will be concatenated.
                complete = os.path.join(matrix_dir, matrix)

                # The first contig is simply renamed to avoid unnecessary copies with the remaining contigs appended.
                if is_first_contig:
                    # TODO: Replace with shutil.move which can handle cross-filesystem moves.
                    os.rename(partial, complete)
                    #futures[index].set_result((futures, index))
                else:
                    # Schedule a process to append the next contig as soon as the previous contig is done.
                    # if matrix.endswith('.vcf'):
                    #     task = functools.partial(_concat_matrix, partial, complete, vcf_metadata_len)
                    # else:
                    #     task = functools.partial(_concat_matrix, partial, complete, 0)
                    # futures[index].add_done_callback(_swap_future(executor, task))

                    if matrix.endswith('.vcf'):
                        _concat_matrix(partial, complete, vcf_metadata_len)
                    else:
                        _concat_matrix(partial, complete, 0)
            is_first_contig = False

            # TODO: Is there enough work here to run in a separate process?
            _concat_snpfasta_contig(tempdirname, contig_name, identifiers, '_missingdata.fasta')
            _concat_snpfasta_contig(tempdirname, contig_name, identifiers, '_bestsnp.fasta')

            # Sum the SampleAnalysis stats preserving the sample_groups order.
            if sample_stats is None:
                sample_stats = sample_stat
            else:
                for sum, analysis in zip(itertools.chain.from_iterable(sample_stats),
                                         itertools.chain.from_iterable(sample_stat)):
                    sum.update(analysis)

        _concat_snpfasta(matrix_dir, tempdirname, 'missingdata.fasta', identifiers, '_missingdata.fasta')
        _concat_snpfasta(matrix_dir, tempdirname, 'bestsnp.fasta', identifiers, '_bestsnp.fasta')

        # TODO: If reference length is returned by matrix_dto, the stats files can be written in parallel
        reference_length = write_general_stats(os.path.join(stats_dir, 'general_stats.tsv'), contig_stats)
        write_sample_stats(os.path.join(stats_dir, 'sample_stats.tsv'), sample_stats, sample_groups, reference_length)

        # The Python documentation says shutdown() is not explicitly needed inside a context manager, but the snpfasta
        # files were not always complete.
        # https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.Executor.shutdown
        # executor.shutdown()
