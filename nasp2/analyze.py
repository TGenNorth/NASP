__author__ = 'jtravis'

from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple, Counter
from struct import Struct
import itertools

# PositionInfo is all the data collected for a single position across all the SampleAnalyses.
PositionInfo = namedtuple('PositionInfo', [
    # True if all samples called A/C/G/T
    'is_all_called',
    'is_reference_clean',
    # True if reference call is in a duplicated region.
    'is_reference_duplicated',
    'is_all_passed_coverage',
    'is_all_passed_proportion',
    # True if all sample calls match.
    'is_all_passed_consensus',
    'is_quality_breadth',
    'is_any_snp',
    # Count of sample calls that differ from the reference.
    'any_snps',
    # True if all sample calls differ from the reference and passed all the filters (coverage, proportion, consensus).
    'is_best_snp',
    #
    'call_str',
    # Count of sample calls that match the reference.
    'called_reference',
    'num_A',
    'num_C',
    'num_G',
    'num_T',
    'num_N',
    'CallWasMade',
    'PassedDepthFilter',
    'PassedProportionFilter',
    'Pattern'
])


def analyze_position(reference_position, dups_position, samples):
    """
    Args:
        reference_position (SampleInfo): A single position from the reference genome.
        dups_position (SampleInfo): A single position from the duplicates file which indicates if the reference is in a duplicate region.
        samples (list of SampleInfo): A single position from all the samples.

    Returns:
        PositionInfo:
    """
    # TODO: Remove debugging threshold variables.
    coverage_threshold = 10.0
    proportion_threshold = .9
    # call_data = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'indel': 0, 'snpcall': 0, 'indelcall': 0, 'refcall': 0,
    # 'callstring': '', 'covstring': '', 'propstring': '', 'called': 0, 'passcov': 0, 'passprop': 0}
    count = Counter({
        # 'A': 0,
        # 'C': 0,
        # 'G': 0,
        # 'T': 0,
        # 'N': 0,
        # 'indel': 0,
        # 'snpcall': 0,
        # 'indelcall': 0,
        # 'refcall': 0,
        'CallWasMade': '',
        'PassedDepthFilter': '',
        'PassedProportionFilter': '',
        'Pattern': ''
    })
    # Assume true until proven otherwise.
    # General Stats
    is_reference_clean = reference_position.simple_call != 'N'
    is_reference_duplicated = dups_position.call == '1'
    is_all_called = True
    is_all_passed_coverage = True
    is_all_passed_proportion = True
    is_all_passed_consensus = True
    is_quality_breadth = True
    is_any_snps = True
    is_best_snp = True
    prev_sample_call = None

    call_str = bytearray(reference_position.call, encoding='utf-8')

    for sample in samples:
        is_pass_coverage = True
        is_pass_proportion = True

        call_str.append(ord(sample.call))

        if sample.call in ['X', 'N']:
            count['CallWasMade'] += 'N'
            is_all_called = False
        else:
            count['CallWasMade'] += 'Y'

        # FIXME: coverage/proportion may be "PASS" if the value was not deprecated in the VCF contig parser.
        # get_proportion/get_coverage

        if sample.coverage is '-':
            # Cannot determine proportion due to missing data.
            count['PassedDepthFilter'] += '-'
        elif sample.coverage == '?':
            # Missing VCF position.
            count['PassedDepthFilter'] += '?'
            is_pass_coverage = False
            is_all_passed_coverage = False
        elif sample.coverage >= coverage_threshold:
            count['PassedDepthFilter'] += 'Y'
        else:
            count['PassedDepthFilter'] += 'N'
            is_pass_coverage = False
            is_all_passed_coverage = False

        # Cannot determine proportion due to missing data.
        if sample.proportion is '-':
            count['PassedProportionFilter'] += '-'
        # Missing VCF position.
        elif sample.proportion == '?':
            count['PassedProportionFilter'] += '?'
            is_pass_proportion = False
            is_all_passed_proportion = False
        elif sample.proportion >= proportion_threshold:
            count['PassedProportionFilter'] += 'Y'
        else:
            count['PassedProportionFilter'] += 'N'
            is_pass_proportion = False
            is_all_passed_proportion = False

        # Samples have consensus if they all have the same call at the current position.
        if prev_sample_call is None:
            prev_sample_call = sample.simple_call
        elif prev_sample_call != sample.simple_call:
            # FIXME: all samples must also pass called, coverage, and proportion filters
            is_all_passed_consensus = False

        count[sample.simple_call] += 1

        if sample.simple_call not in ['X', 'N'] and is_pass_coverage and is_pass_proportion and\
                    reference_position.simple_call != 'N' and dups_position.call != '1':
            if sample.call == reference_position.call:
                count['called_reference'] += 1
            else:
                count['called_snp'] += 1

    # FIXME: all_passed_proportion, all_passed_consensus, is_quality_breadth, is_best_snp

    return PositionInfo(
        # General Stats
        is_all_called=is_all_called,
        is_reference_clean=is_reference_clean,
        is_reference_duplicated=is_reference_duplicated,
        is_all_passed_coverage=is_all_passed_coverage,
        is_all_passed_proportion=is_all_passed_proportion,
        is_all_passed_consensus=is_all_passed_consensus,
        is_quality_breadth=is_quality_breadth,
        is_any_snp=count['called_snp'] > 1,
        is_best_snp='',
        # Matrix
        call_str=call_str,
        called_reference=count['called_reference'],
        num_A=count['A'],
        num_C=count['C'],
        num_G=count['G'],
        num_T=count['T'],
        num_N=count['N'],
        any_snps=count['any_snps'],
        CallWasMade=bytearray(count['CallWasMade'], encoding='utf-8'),
        PassedDepthFilter=bytearray(count['PassedDepthFilter'], encoding='utf-8'),
        PassedProportionFilter=bytearray(count['PassedProportionFilter'], encoding='utf-8'),
        Pattern=bytearray(count['Pattern'], encoding='utf-8')
    )


def analyze_contig(reference_contig, dups_contig, sample_analyses, offset):
    """
    analyze_contig expects reference_contig, dups_contig, and sample_analyses to represent the same contig.
    It reads all of their positions and writes their analysis to the outfile.

    Args:
        outfile (str): Path of the binary file to write analysis data.
        offset (int): File position to start writing analysis data.
        reference_contig (FastaContig):
        dups_contig (FastaContig):
        sample_analyses (tuple of SampleAnalysis):

    Returns:
        (str, collections.Counter): The contig name and statistics gathered across all the contig positions.
    """
    contigs = (sample_analysis.get_contig(reference_contig.name).positions for sample_analysis in sample_analyses)

    string_pattern = str(len(sample_analyses)) + 's'
    s = Struct('HHHHHHH?' + string_pattern*5)

    # Open the file without truncating. The file must already exist.
    with open('partial.nasp', 'r+b') as handle:
        handle.seek(offset)
        contig_stats = Counter()
        for position in map(analyze_position, reference_contig.positions, dups_contig.positions, zip(*contigs)):
            contig_stats['reference_length'] += 1
            contig_stats['reference_clean'] += 1 if position.is_reference_clean else 0
            contig_stats['reference_duplicated'] += 1 if position.is_reference_duplicated else 0
            contig_stats['all_called'] += 1 if position.is_all_called else 0
            contig_stats['all_passed_coverage'] += 1 if position.is_all_passed_coverage else 0
            contig_stats['all_passed_proportion'] += 1 if position.is_all_passed_proportion else 0
            contig_stats['all_passed_consensus'] += 1 if position.is_all_passed_consensus else 0
            contig_stats['quality_breadth'] += 1 if position.is_quality_breadth else 0
            contig_stats['any_snps'] += 1 if position.is_any_snp else 0
            contig_stats['best_snps'] += 1 if position.is_best_snp else 0

            # # TODO: Remove print
            # print(contig_stats['reference_length'], position)
            # if contig_stats['reference_length'] == 392:
            #     exit(1)
            #
            # foo=position.PassedProportionFilter
            # foo.append(ord('\n'))
            # handle.write(foo)

            # handle.write(s.pack(position.called_reference, position.num_A, position.num_C, position.num_G,
            #                     position.num_T, position.num_N, position.any_snps, position.is_reference_duplicated,
            #                     position.call_str, position.CallWasMade,
            #                     position.PassedDepthFilter, position.PassedProportionFilter, position.Pattern))
    return reference_contig.name, contig_stats


def analyze_samples(reference_fasta, reference_dups, sample_analyses):
    string_pattern = str(len(sample_analyses)) + 's'
    s = Struct('HHHHHH?' + string_pattern*5)

    # for k, v in itertools.groupby(sample_analyses, lambda x: x.name):
    #     print(k, v)
    # exit(1)

    with open('partial.nasp', 'w') as handle:
        # handle.write('\t'.join((sample_analysis.identifier for sample_analysis in sample_analyses)) + '\n')
        # handle.write('\t'.join((contig.name for contig in matrix_parameters.reference_fasta.contigs)) + '\n')
        header_offset = handle.tell()

    partial_file_positions = itertools.chain((header_offset,), (len(contig)*s.size for contig in reference_fasta.contigs))



    # with ProcessPoolExecutor() as executor:
    #     for result in executor.map(analyze_contig, matrix_parameters.reference_fasta.contigs, matrix_parameters.reference_dups.contigs, itertools.repeat(sample_analyses), partial_file_positions):
    #         print(result)
    for result in map(analyze_contig, reference_fasta.contigs, reference_dups.contigs, itertools.repeat(sample_analyses), partial_file_positions):
        print(result)

    # with open('partial.nasp', 'br') as handle:
    #     sample_analyses = str(handle.readline(), encoding='utf-8').rstrip().split('\t')
    #     contigs = str(handle.readline(), encoding='utf-8').rstrip().split('\t')
    #     print(sample_analyses)
    #     print(contigs)
    #     for chunk in iter(lambda: handle.read(s.size), b''):
    #         print(s.unpack(chunk))