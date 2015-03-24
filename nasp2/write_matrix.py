"""
write_matrix handles all disk IO operations.

TODO: Most of the functions are written as coroutines that could be used improve performance by using asyncio to
continue processing instead of waiting for data to be read/written to/from disk.
"""
__author__ = 'jtravis'

import os
import csv
from collections import Counter
from contextlib import ExitStack


def get_vcf_metadata(nasp_version, identifiers, contigs):
    """
    Args:
        nasp_version (str): v1.0.0
        identifiers (tuple): contig_name::aligner,snpcaller
        contigs (tuple of Contig):

    Return:
        str: VCF file metadata
    """
    vcf_metadata = "##fileFormat=VCFv4.2\n##source=NASPv{0}\n".format(nasp_version)
    vcf_metadata += "\n".join(
        "##contig=<ID=\"{0}\",length={1}>\n".format(contig['name'], contig['length']) for contig in contigs)
    vcf_metadata += "\n".join(
        "##SAMPLE=<ID=\"{0}\",Genomes=\"{0}\",Mixture=1.0>\n".format(identifier) for identifier in identifiers)
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
        type (str): 'all_callable', 'missing_data', 'best_snp', 'vcf'
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
    if type in ['all_callable', 'missing_data']:
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
    import itertools
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


def write_missingdata_vcf(directory, contig_name, identifiers):
    """
    Args:
        directory (str):
        contig_name (str):
        identifiers:
    """
    coverage_threshold = 10
    proportion_threshold = 0.9

    with open('{0}_missingdata.vcf'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # handle.write(get_vcf_metadata(version, identifiers, contigs))
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
                'FILTER': _vcf_filter_column(coverage_threshold, proportion_threshold, row.is_all_passed_consensus, row.is_all_passed_proportion),
                # TODO: AN is the number of snps + 1 for the reference.
                # TODO: Add #indel stat to NS
                'INFO': 'AN={0};NS={1}'.format(len(alts) + 1, row.called_reference + row.called_snp),
                'FORMAT': 'GT:FT'
            }
            # Match each analysis with its analysis column.
            line.update({k: v for k, v in zip(identifiers, _vcf_analysis_column(row.Pattern, row.all_sample_stats))})
            writer.writerow(line)


def write_bestsnp_vcf(directory, contig_name, identifiers):
    """
    Args:
        directory (str):
        contig_name (str):
        identifiers:
    """
    coverage_threshold = 10
    proportion_threshold = 0.9

    with open('{0}_bestsnp.vcf'.format(os.path.join(directory, contig_name)), 'w') as handle:
        # handle.write(get_vcf_metadata(version, identifiers, contigs))
        writer = csv.DictWriter(handle, fieldnames=get_header('vcf', identifiers), delimiter='\t')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_best_snp:
                continue

            ref = row.call_str[0]
            alts = set(row.call_str[1:])
            alts.discard(ref)

            line = {
                '#CHROM': contig_name,
                'POS': position,
                'ID': '.',
                'REF': ref,
                'ALT': ','.join(alts) or '.',
                'QUAL': '.',
                'FILTER': _vcf_filter_column(coverage_threshold, proportion_threshold, row.is_all_passed_consensus, row.is_all_passed_proportion),
                # TODO: AN is the number of snps + 1 for the reference.
                # TODO: Add #indel stat to NS
                'INFO': 'AN={0};NS={1}'.format(len(alts) + 1, row.called_reference + row.called_snp),
                'FORMAT': 'GT:FT'
            }
            # Match each analysis with its analysis column.
            line.update({k: v for k, v in zip(identifiers, _vcf_analysis_column(row.Pattern, row.all_sample_stats))})
            writer.writerow(line)


def write_sample_stats(filepath, sample_stats, sample_groups, reference_length):
    """
    Args:
        filepath (str):
        sample_stats (list of lists of Counters):
        sample_groups (tuple of tuples of SampleAnalysis):
        reference_length (int): Total number of positions in the reference.
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

    Return:
        int: Total reference length.
    """
    # Sum contig stats
    whole_genome_stats = Counter({'Contig': ''})
    for contig_stat in contig_stats:
        whole_genome_stats.update(contig_stat)
    whole_genome_stats['Contig'] = 'Whole Genome'

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
                contig_stat[stat + ' (%)'] = "{0:.2f}%".format(contig_stat[stat] / contig_stat['reference_length'] * 100)
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
        writer = csv.DictWriter(handle, fieldnames=get_header('all_callable', identifiers), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1
            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': row.call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
                '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
                '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
                '#A': row.num_A,
                '#C': row.num_C,
                '#G': row.num_G,
                '#T': row.num_T,
                # TODO: replace with n/a
                '#Indel': '0',
                '#NXdegen': row.num_N,
                'Contig': contig_name,
                'Position': position,
                'InDupRegion': row.is_reference_duplicated,
                'SampleConsensus': row.is_all_passed_consensus,
                'CallWasMade': row.CallWasMade,
                'PassedDepthFilter': row.PassedDepthFilter,
                'PassedProportionFilter': row.PassedProportionFilter,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update({k: v for k, v in zip(identifiers, row.call_str[1:])})
            writer.writerow(line)


def write_missingdata_matrix(directory, contig_name, identifiers):
    """
    Args:
        directory(str):
        contig_name (str):
        identifiers:
    """
    with open('{0}_missingdata.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=get_header('missing_data', identifiers), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_missing_matrix:
                continue

            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.masked_call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': row.masked_call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
                '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
                '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
                '#A': row.num_A,
                '#C': row.num_C,
                '#G': row.num_G,
                '#T': row.num_T,
                # TODO: replace with n/a
                '#Indel': '0',
                '#NXdegen': row.num_N,
                'Contig': contig_name,
                'Position': position,
                'InDupRegion': row.is_reference_duplicated,
                'SampleConsensus': row.is_all_passed_consensus,
                'CallWasMade': row.CallWasMade,
                'PassedDepthFilter': row.PassedDepthFilter,
                'PassedProportionFilter': row.PassedProportionFilter,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update({k: v for k, v in zip(identifiers, row.masked_call_str[1:])})
            writer.writerow(line)


def write_bestsnp_matrix(directory, contig_name, sample_groups):
    """
    Args:
        directory (str):
        contig_name (str):
        sample_groups
    """
    sample_names = tuple(sample[0].name for sample in sample_groups)

    # first_analysis_index is a list of the index of the first analysis for each sample in the call string
    # Unlike the
    first_analysis_index = []
    # num_analyses starts at 1 to skip the index call
    num_analyses = 1
    for sample in sample_groups:
        first_analysis_index.append(num_analyses)
        num_analyses += len(sample)

    with open('{0}_bestsnp.tsv'.format(os.path.join(directory, contig_name)), 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=get_header('best_snp', sample_names), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_best_snp:
                continue

            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': row.call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
                '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
                '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
                '#A': row.num_A,
                '#C': row.num_C,
                '#G': row.num_G,
                '#T': row.num_T,
                # TODO: replace with n/a
                '#Indel': '0',
                '#NXdegen': row.num_N,
                'Contig': contig_name,
                'Position': position,
                'InDupRegion': row.is_reference_duplicated,
                'SampleConsensus': row.is_all_passed_consensus,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update(
                {sample_name: row.call_str[index] for sample_name, index in zip(sample_names, first_analysis_index)})

            writer.writerow(line)


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
        writer = csv.DictWriter(handle, fieldnames=get_header('best_snp', identifiers), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            # if not row.is_all_quality_breadth:
            #     continue

            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': row.call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - row.CallWasMade.count('N'), num_samples),
                '#PassedDepthFilter': "{0:d}/{1:d}".format(row.passed_coverage_filter, num_samples),
                '#PassedProportionFilter': "{0:d}/{1:d}".format(row.passed_proportion_filter, num_samples),
                '#A': row.num_A,
                '#C': row.num_C,
                '#G': row.num_G,
                '#T': row.num_T,
                # TODO: replace with n/a
                '#Indel': '0',
                '#NXdegen': row.num_N,
                'Contig': contig_name,
                'Position': position,
                'InDupRegion': row.is_reference_duplicated,
                'SampleConsensus': row.is_all_passed_consensus,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update({k: v for k, v in zip(identifiers, row.call_str[1:])})

            writer.writerow(line)


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
        files = tuple(stack.enter_context(open('{0}_{1}_bestsnp.fasta'.format(os.path.join(directory, contig_name), identifier), 'w')) for identifier in identifiers)

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
        files = tuple(stack.enter_context(open('{0}_{1}_missingdata.fasta'.format(os.path.join(directory, contig_name), identifier), 'w')) for identifier in identifiers)

        while True:
            row = yield

            if not row.is_missing_matrix:
                continue

            for file, call in zip(files, row.masked_call_str):
                if call != 'X':
                    file.write(call)
                else:
                    file.write('N')