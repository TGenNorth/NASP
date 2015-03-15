__author__ = 'jtravis'

import csv
import functools
from collections import Counter


def get_vcf_metadata(nasp_version, identifiers, contigs):
    """
    Args:
        nasp_version (str): v1.0.0
        identifiers (tuple): contig_name::aligner,snpcaller
        contigs (tuple of Contig):
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


# def write_sample_stats(filepath):
# sample_stats = None
#
# try:
# while True:
# all_sample_stats = yield
#
# if sample_stats is None:
# sample_stats = all_sample_stats
#             for i, sample in enumerate(all_sample_stats):
#                 for j, analysis in enumerate(sample):
#                     sample_stats[i][j].update(analysis)
#
#     finally:
#         print(sample_stats)
#         pass


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

    any = Counter({'Sample': '[any]'})
    all = Counter({'Sample': '[all]'})

    for sample_stat, sample in zip(sample_stats, sample_groups):
        any.update(sample_stat[0])
        all.update(sample_stat[1])
        sample_stat[0]['Sample'] = sample[0].name
        sample_stat[0]['Sample::Analysis'] = '[any]'
        sample_stat[1]['Sample'] = sample[0].name
        sample_stat[1]['Sample::Analysis'] = '[all]'

        for stat in fieldnames[2::2]:
            sample_stat[0][stat + ' (%)'] = "{0:.2f}%".format(sample_stat[0][stat] / reference_length * 100)
            sample_stat[1][stat + ' (%)'] = "{0:.2f}%".format(sample_stat[1][stat] / reference_length * 100)

    # Calculate overall stat percentages
    for stat in fieldnames[2::2]:
        any[stat + ' (%)'] = "{0:.2f}%".format(any[stat] / reference_length * 100)
        all[stat + ' (%)'] = "{0:.2f}%".format(any[stat] / reference_length * 100)

    with open(filepath, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        handle.write('\n')
        # Write the any and all summaries genome.
        writer.writerow(any)
        writer.writerow(all)

        # Join sample names and analysis identifiers with the analysis stats before writing them to file.
        for sample_stat, sample in zip(sample_stats, sample_groups):
            handle.write('\n')
            writer.writerow(sample_stat[0])
            writer.writerow(sample_stat[1])
            # Join sample identifiers with sample data then write to file.
            for analysis_stats, identifier in zip(sample_stat[2:], sample):
                # Calculate per-analysis stat percentages
                for stat in fieldnames[2::2]:
                    analysis_stats[stat + ' (%)'] = "{0:.2f}%".format(analysis_stats[stat] / reference_length * 100)
                analysis_stats['Sample'] = identifier.name
                analysis_stats['Sample::Analysis'] = identifier.identifier
                writer.writerow(analysis_stats)


def _sum_contig_stats(summation, contig_stat):
    summation.update(contig_stat)
    return summation


def write_general_stats(filepath, contig_stats):
    """
    Args:
        filepath (str):
        contig_stats (tuple of Counter):

    Return:
        int: Total reference length.
    """
    whole_genome_stats = functools.reduce(_sum_contig_stats, contig_stats, Counter({'Contig': ''}))
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
                contig_stat[stat + ' (%)'] = "{0:.2f}%".format(contig_stat[stat] / reference_length * 100)
            writer.writerow(contig_stat)

    return whole_genome_stats['reference_length']


def write_master_matrix(filepath, contig_name, identifiers):
    """
    Args:
        filepath (str): Path to the output file.
        contig_name (str): Name
        identifiers (tuple of
    """
    with open(filepath, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=get_header('all_callable', identifiers), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1
            call_str = str(row.call_str, encoding='utf-8')
            call_was_made = str(row.CallWasMade, encoding='utf-8')
            passed_depth_filter = str(row.PassedDepthFilter, encoding='utf-8')
            passed_proportion_filter = str(row.PassedProportionFilter, encoding='utf-8')
            # num_samples is the number of analyses not including the reference.
            num_samples = len(call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - call_was_made.count('N'), num_samples),
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
                'CallWasMade': call_was_made,
                'PassedDepthFilter': passed_depth_filter,
                'PassedProportionFilter': passed_proportion_filter,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update({k: v for k, v in zip(identifiers, call_str[1:])})
            writer.writerow(line)


def write_missingdata_matrix(filepath, contig_name, identifiers):
    with open(filepath, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=get_header('missing_data', identifiers), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not (row.is_missing_matrix):
                continue

            # call_str = str(row.call_str, encoding='utf-8')
            call_was_made = str(row.CallWasMade, encoding='utf-8')
            passed_depth_filter = str(row.PassedDepthFilter, encoding='utf-8')
            passed_proportion_filter = str(row.PassedProportionFilter, encoding='utf-8')
            # num_samples is the number of analyses not including the reference.
            num_samples = len(row.masked_call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': row.masked_call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - call_was_made.count('N'), num_samples),
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
                'CallWasMade': call_was_made,
                'PassedDepthFilter': passed_depth_filter,
                'PassedProportionFilter': passed_proportion_filter,
                'Pattern': "".join(row.Pattern)
            }
            # Match each base call with its sample analysis column.
            line.update({k: v for k, v in zip(identifiers, row.masked_call_str[1:])})
            writer.writerow(line)


# def write_missingdata_vcf(filepath, identifiers, contigs, version):
#     with open(filepath, 'w') as handle:
#         handle.write(get_vcf_metadata(version, identifiers, contigs))
#         writer = csv.DictWriter(handle, fieldnames=get_header('vcf', identifiers), delimiter='\t')
#         writer.writeheader()
#         position = 0
#         while True:
#             row = yield
#             position += 1
#             writer.writerow({
#                 '#CHROM': '',
#                 'POS': position,
#                 'ID': '',
#                 'REF': '',
#                 'ALT': '',
#                 'QUAL': '',
#                 'FILTER': '',
#                 'INFO': '',
#                 'FORMAT': ''
#             })


def write_bestsnp_matrix(filepath, contig_name, sample_groups):
    sample_names = tuple(sample[0].name for sample in sample_groups)

    # first_analysis_index is a list of the index of the first analysis for each sample in the call string
    # Unlike the
    first_analysis_index = []
    # num_analyses starts at 1 to skip the index call
    num_analyses = 1
    for sample in sample_groups:
        first_analysis_index.append(num_analyses)
        num_analyses += len(sample)

    with open(filepath, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=get_header('best_snp', sample_names), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        position = 0
        while True:
            row = yield
            position += 1

            if not row.is_best_snp:
                continue

            call_str = str(row.call_str, encoding='utf-8')
            call_was_made = str(row.CallWasMade, encoding='utf-8')
            # num_samples is the number of analyses not including the reference.
            num_samples = len(call_str) - 1

            line = {
                'LocusID': "{0}::{1}".format(contig_name, position),
                'Reference': call_str[0],
                '#SNPcall': row.called_snp,
                # TODO: replace with n/a
                '#Indelcall': '0',
                '#Refcall': row.called_reference,
                '#CallWasMade': "{0:d}/{1:d}".format(num_samples - call_was_made.count('N'), num_samples),
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
                {sample_name: call_str[index] for sample_name, index in zip(sample_names, first_analysis_index)})

            writer.writerow(line)


            # def write_fasta(filepath, contig_name):
            #     with open(filepath, 'w') as handle:
            #         handle.write('>' + contig_name + '\n')
            #         num_chars = 0
            #         while True:
            #             row = yield
            #             handle.write(row)
            #             num_chars += 1
            #             if num_chars == 80:
            #                 handle.write('\n')
            #                 num_chars = 0

            # for matrix_format in all_vcfs:
            # if len(encountered_calls) > 0:
            #             matrix_format['linetowrite'] += ",".join(encountered_calls)
            #         else:
            #             matrix_format['linetowrite'] += "."
            #         matrix_format['linetowrite'] += "\t.\tPASS\tAN={0};NS={1}\tGT:FT".format(len(encountered_calls)+1, str(call_data['snpcall']+call_data['indelcall']+call_data['refcall']))
            #         for vcf_current_data in vcf_pending_data:
            #             matrix_format['linetowrite'] += "\t{0}:".format(str(vcf_current_data['GT']))
            #             if not vcf_current_data['was_called']:
            #                 matrix_format['linetowrite'] += "NoCall"
            #             elif not vcf_current_data['passed_coverage']:
            #                 matrix_format['linetowrite'] += "CovFail"
            #             elif not vcf_current_data['passed_proportion']:
            #                 matrix_format['linetowrite'] += "PropFail"
            #             else:
            #                 matrix_format['linetowrite'] += "PASS"
            #         matrix_format['linetowrite'] += "\n"


            # def send_to_matrix_handles(self, matrix_formats):
            #     """
            #     Writes headers and handles per-matrix logic.  Calls _write_matrix_line
            #     to handle the per-line computation and analysis.
            #     """
            #     for matrix_format in matrix_formats:
            #         if matrix_format['dataformat'] == 'matrix':
            #             matrix_format['handle'].write("LocusID\tReference\t")
            #             for genome in self._genomes:
            #                 matrix_format['handle'].write("{0}\t".format(genome.identifier()))
            #             for genome_path in self._failed_genomes:
            #                 matrix_format['handle'].write("{0}\t".format(genome_path))
            #             if matrix_format['filter'] == 'bestsnp' or matrix_format['filter'] == 'includeref':
            #                 matrix_format['handle'].write(
            #                     "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tPattern\tPattern#\n")
            #             else:
            #                 # must be all callable or missing data
            #                 matrix_format['handle'].write(
            #                     "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\t" +
            #                     "Contig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\tPattern#\n")
            #             matrix_format['linetowrite'] = ''
            #         elif matrix_format['dataformat'] == 'fasta':
            #             matrix_format['fastadata'] = GenomeStatus()
            #             for genome in self._genomes:
            #                 matrix_format['fastadata'].add_contig(genome.identifier())
            #         elif matrix_format['dataformat'] == 'vcf':
            #             matrix_format['handle'].write("##fileFormat=VCFv4.2\n##source=NASPv{0}\n".format(__version__))
            #             for current_contig in self._reference.get_contigs():
            #                 matrix_format['handle'].write("##contig=<ID=\"{0}\",length={1}>\n".format(current_contig,
            #                                                                                           self._reference.get_contig_length(
            #                                                                                               current_contig)))
            #             for genome in self._genomes:
            #                 matrix_format['handle'].write(
            #                     "##SAMPLE=<ID=\"{0}\",Genomes=\"{0}\",Mixture=1.0>\n".format(genome.identifier()))
            #             matrix_format['handle'].write(
            #                 "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n")
            #             matrix_format['handle'].write(
            #                 "##FILTER=<ID=NoCall,Description=\"No call for this sample at this position\">\n")
            #             matrix_format['handle'].write(
            #                 "##FILTER=<ID=CovFail,Description=\"Insufficient depth of coverage for this sample at this position\">\n")
            #             matrix_format['handle'].write(
            #                 "##FILTER=<ID=PropFail,Description=\"Insufficient proportion of reads were variant for this sample at this position\">\n")
            #             matrix_format['handle'].write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            #             matrix_format['handle'].write(
            #                 "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filters that failed for this sample at this position\">\n")
            #             matrix_format['handle'].write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            #             for genome in self._genomes:
            #                 matrix_format['handle'].write("\t" + genome.identifier())
            #             matrix_format['handle'].write("\n")
            #     # Key None stores next unused
            #     pattern_data = {None: 1}
            #     for current_contig in self.get_contigs():
            #         for current_pos in range(1, self._reference.get_contig_length(current_contig) + 1):
            #             self._format_matrix_line(current_contig, current_pos, matrix_formats, pattern_data)
            #             for matrix_format in matrix_formats:
            #                 if matrix_format['dataformat'] == 'matrix' and matrix_format['linetowrite'] is not None:
            #                     matrix_format['handle'].write(matrix_format['linetowrite'])
            #                     matrix_format['linetowrite'] = ''
            #                 if matrix_format['dataformat'] == 'vcf' and matrix_format['linetowrite'] is not None:
            #                     matrix_format['handle'].write(matrix_format['linetowrite'])
            #                     matrix_format['linetowrite'] = ''
            #     for matrix_format in matrix_formats:
            #         if matrix_format['dataformat'] == 'fasta':
            #             matrix_format['fastadata'].send_to_fasta_handle(matrix_format['handle'])