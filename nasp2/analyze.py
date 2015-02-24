"""
The analyze module takes a collection of SampleAnalyses
"""

__author__ = 'jtravis'

from collections import namedtuple, Counter
from concurrent.futures import ProcessPoolExecutor
import itertools

from nasp2.write_matrix import write_master_matrix


# TODO: Remove profile. It is a dummy wrapper to leave @profile decorators in place when not profiling.
def profile(fn):
    return fn

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
    'is_all_quality_breadth',
    'is_any_snp',
    # Count of sample calls that differ from the reference.
    'any_snps',
    # True if all sample calls differ from the reference and passed all the filters (coverage, proportion, consensus).
    'is_best_snp',
    #
    'call_str',
    # Count of sample calls that match the reference.
    'called_reference',
    'passed_coverage_filter',
    'passed_proportion_filter',
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


# @profile
def sample_positions(contig_name, sample_groups):
    """
    sample_positions is a generator that yields positions across all the SampleAnalyses
    for the given contig grouped by sample name. This is the magic sauce that allows all
    the files to be read in unison as if they were a single file. The element

    TODO: Add See Also section referencing SampleAnalysis.

    Args:
        contig_name (str): Name of a contig from the reference.
        sample_groups (tuple of SampleAnalysis tuples):

    Return:
        tuple of Position tuples generator: Yields tuples of Position grouped by sample name.

    Example:

        sample_groups is the SampleAnalyses grouped by the name of the sample they are analyzing.

        +------------------------------------------------------------------------------------+
        |              Sample Group                |             Sample Group                |
        +==========================================+=========================================+
        |  SampleAnalysis    | SampleAnalysis      |  SampleAnalysis   | SampleAnalysis      |
        +------------------------------------------+-----------------------------------------+

        contig_groups is a mirror of sample_groups with each SampleAnalysis replaced
        with a contig position generator corresponding to the reference contig.

        +------------------------------------------------------------------------------------+
        |              Sample Group                |              Sample Group               |
        +==========================================+=========================================+
        | <contig_positions> | <contig_positions>  | <contig_positions> | <contig_positions> |
        +------------------------------------------+-----------------------------------------+

        The sample positions generator yields a tuple of SampleInfo tuples:
        ( (SampleInfo, SampleInfo,), (SampleInfo, SampleInfo), )

        +------------------------------------------------------------------------------------+
        |              Sample Group                |              Sample Group               |
        +==========================================+=========================================+
        | Position           | Position            | Position           | Position           |
        +------------------------------------------+-----------------------------------------+
    """
    # sample_groups is converted to contig_groups by calling get_contig() on each SampleAnalysis.
    contig_groups = tuple(
        tuple(map(lambda sample_analysis: sample_analysis.get_contig(contig_name).positions, sample_group))\
        for sample_group in sample_groups
    )

    while True:
        # The inner tuple is all the sample analyses for a single sample.
        # The outer tuple is all the sample analyses.
        # ( (SampleInfo, SampleInfo,), (SampleInfo, SampleInfo), )
        yield tuple(tuple(map(next, contig_group)) for contig_group in contig_groups)


# @profile
def analyze_position(reference_position, dups_position, samples):
    """
    analyze_position compares all the SampleAnalyses at a single position.

    Args:
        reference_position (SampleInfo): A single position from the reference genome.
        dups_position (SampleInfo): A single position from the duplicates file which indicates if the reference is in a duplicate region.
        samples (tuple of SampleInfo tuples): A single position across all SampleAnalyses grouped by sample name.

    Returns:
        PositionInfo: All the information gathered
    """
    # TODO: Remove debugging threshold variables.
    coverage_threshold = 10.0
    proportion_threshold = .9


    # General Stats
    # Assume true until proven otherwise.
    is_reference_clean = reference_position.simple_call != 'N'
    is_reference_duplicated = dups_position.call == '1'
    is_all_called = True
    is_all_passed_coverage = True
    is_all_passed_proportion = True
    is_all_passed_consensus = True
    # General Stats quality breadth passes if all the following conditions for the position are true:
    # - Reference called and not duplicated.
    # - All analyses called and passed all filters.
    # - All samples have consensus within their aligner/snpcaller combinations.
    is_all_quality_breadth = True
    # General Stats quality breadth filter passes and all analyses are SNPs.
    is_best_snp = True

    # TODO: consolidate is_all_sample_consensus vs is_all_passed_consensus
    is_all_sample_consensus = True

    # position_stats is a counter across all SampleAnalyses at the current position. It is used in the Matrix files.
    position_stats = Counter({
        # TODO: verify 'N' represents NoCall, AnyCall, Deletion, etc.
        # Tally of A/C/G/T calls where everything else is called N representing
        # 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,

        # The character at each position of the following strings indicates whether a corresponding SampleAnalysis
        # passed the filter:
        'CallWasMade': '',
        'PassedDepthFilter': '',
        'PassedProportionFilter': '',
        'Pattern': ''
    })

    # call_str is the base call for each SampleAnalysis starting with the reference.
    call_str = bytearray(reference_position.call, encoding='utf-8')

    # all_sample_stats is a list of sample_stats. It represents all the SampleAnalysis stats grouped by sample name.
    all_sample_stats = []
    for sample in samples:
        # Except for the first entry, sample_stats is composed of analysis_stats dictionaries for each aligner/snpcaller
        # combination used on the sample. The first entry is a special consensus boolean that is true if all the
        # analyses have the same DNA base call at this position.
        # [is_consensus, analysis_stats, analysis_stats, analysis_stats, ...]
        sample_stats = [True]
        is_sample_consensus = True
        prev_analysis_call = None
        for analysis in sample:
            # analysis_stats are unique to each SampleAnalysis. It is used in the Sample Statistics file.
            analysis_stats = {
                'was_called': True,
                'passed_coverage_filter': True,
                'passed_proportion_filter': True,
                'quality_breadth': True,
                'called_reference': False,
                'called_snp': False,
                'called_dgen': False
            }

            call_str.append(ord(analysis.call))

            # TODO: Must also be called (not X/N) and pass coverage/proportion filters
            # A sample has position consensus if all of its SampleAnalysis combinations have the same call.
            if prev_analysis_call is None:
                prev_analysis_call = analysis.simple_call
            elif prev_analysis_call != analysis.simple_call:
                is_sample_consensus = False
                is_all_passed_consensus = False
                is_all_quality_breadth = False

            # X means no value, N means any value.
            if analysis.call in ['X', 'N']:
                position_stats['CallWasMade'] += 'N'
                analysis_stats['was_called'] = False
                is_all_called = False
                is_all_quality_breadth = False
            else:
                position_stats['was_called'] += 1
                position_stats['CallWasMade'] += 'Y'

            # FIXME: coverage/proportion may be "PASS" if the value was not deprecated in the VCF contig parser.
            # get_proportion/get_coverage

            # Cannot determine proportion due to missing data.
            if analysis.coverage is '-':
                position_stats['PassedDepthFilter'] += '-'
                position_stats['passed_coverage_filter'] += 1
            # Missing VCF position.
            elif analysis.coverage == '?':
                position_stats['PassedDepthFilter'] += '?'
                analysis_stats['passed_coverage_filter'] = False
                is_all_passed_coverage = False
            elif analysis.coverage >= coverage_threshold:
                position_stats['PassedDepthFilter'] += 'Y'
                position_stats['passed_coverage_filter'] += 1
            else:
                position_stats['PassedDepthFilter'] += 'N'
                analysis_stats['passed_proportion_filter'] = False
                is_all_passed_coverage = False

            # Cannot determine proportion due to missing data.
            if analysis.proportion is '-':
                position_stats['PassedProportionFilter'] += '-'
                position_stats['passed_proportion_filter'] += 1
            # Missing VCF position.
            elif analysis.proportion == '?':
                position_stats['PassedProportionFilter'] += '?'
                analysis_stats['passed_proportion_filter'] = False
                is_all_passed_proportion = False
            elif analysis.proportion >= proportion_threshold:
                position_stats['PassedProportionFilter'] += 'Y'
                position_stats['passed_proportion_filter'] += 1
            else:
                position_stats['PassedProportionFilter'] += 'N'
                analysis_stats['passed_proportion_filter'] = False
                is_all_passed_proportion = False

            # Count the number of A/C/G/T calls. Any other value is N. See SampleInfo.simple_call()
            position_stats[analysis.simple_call] += 1

            # TODO: Clarify the purpose of this if statement.
            # Only count significant measurements.
            if analysis_stats['passed_coverage_filter'] and analysis_stats['passed_proportion_filter'] and\
                    reference_position.simple_call != 'N' and dups_position.call != '1':

                analysis_stats['quality_breadth'] = True

                if analysis.call in ['X', 'N']:
                    position_stats['called_degen'] += 1
                    analysis_stats['called_degen'] = True
                elif analysis.call == reference_position.call:
                    position_stats['called_reference'] += 1
                    analysis_stats['called_reference'] = True
                else:
                    position_stats['called_snp'] += 1
                    analysis_stats['called_snp'] = True
            else:
                is_all_quality_breadth = False

            sample_stats.append(analysis_stats)

        # Consensus for the last sample is available after checking all its analyses.
        sample_stats[0] = is_sample_consensus

        all_sample_stats.append(sample_stats)

    # FIXME: all_passed_proportion, all_passed_consensus, is_quality_breadth, is_best_snp

    return PositionInfo(
        # General Stats
        is_all_called=is_all_called,
        is_reference_clean=is_reference_clean,
        is_reference_duplicated=is_reference_duplicated,
        is_all_passed_coverage=is_all_passed_coverage,
        is_all_passed_proportion=is_all_passed_proportion,
        is_all_passed_consensus=is_all_passed_consensus,
        is_all_quality_breadth=is_all_quality_breadth,
        is_any_snp=position_stats['called_snp'] > 0,
        is_best_snp='',
        # TODO: return all_sample_stats
        # NASP Master Matrix
        call_str=call_str,
        called_reference=position_stats['called_reference'],
        passed_coverage_filter=position_stats['passed_coverage_filter'],
        passed_proportion_filter=position_stats['passed_proportion_filter'],
        num_A=position_stats['A'],
        num_C=position_stats['C'],
        num_G=position_stats['G'],
        num_T=position_stats['T'],
        num_N=position_stats['N'],
        # TODO: Remove any_snp if unused in any output file.
        # is_any_snp is a similar looking stat in General Stats. If this stat is required, add documentation explaining
        # the difference.
        any_snps=position_stats['any_snps'],
        CallWasMade=bytearray(position_stats['CallWasMade'], encoding='utf-8'),
        PassedDepthFilter=bytearray(position_stats['PassedDepthFilter'], encoding='utf-8'),
        PassedProportionFilter=bytearray(position_stats['PassedProportionFilter'], encoding='utf-8'),
        Pattern=bytearray(position_stats['Pattern'], encoding='utf-8')
    )


# @profile
def analyze_contig(reference_contig, dups_contig, sample_groups):
    """
    analyze_contig expects reference_contig, dups_contig, and sample_analyses to represent the same contig from
    different sources. It reads all of their positions and writes an analysis to the outfile.

    Args:
        outfile (str): Path of the binary file to write analysis data.
        offset (int): File position to start writing analysis data.
        reference_contig (FastaContig):
        dups_contig (FastaContig):
        sample_groups (tuple of SampleAnalysis tuples): SampleAnalyses grouped by sample name.

    Returns:
        (str, collections.Counter): The contig name and statistics gathered across all the contig positions.
    """
    contig_stats = Counter()
    # FIXME: Must the identifiers be recalculated for every contig?
    identifiers = tuple(analysis.identifier for sample in sample_groups for analysis in sample)

    # Initialize outfile coroutine.
    outfile = write_master_matrix(reference_contig.name + '.tsv', reference_contig.name, identifiers)
    outfile.send(None)

    for position in map(analyze_position, reference_contig.positions, dups_contig.positions, sample_positions(reference_contig.name, sample_groups)):
        contig_stats['reference_length'] += 1
        contig_stats['reference_clean'] += 1 if position.is_reference_clean else 0
        contig_stats['reference_duplicated'] += 1 if position.is_reference_duplicated else 0
        contig_stats['all_called'] += 1 if position.is_all_called else 0
        contig_stats['all_passed_coverage'] += 1 if position.is_all_passed_coverage else 0
        contig_stats['all_passed_proportion'] += 1 if position.is_all_passed_proportion else 0
        contig_stats['all_passed_consensus'] += 1 if position.is_all_passed_consensus else 0
        contig_stats['quality_breadth'] += 1 if position.is_all_quality_breadth else 0
        contig_stats['any_snps'] += 1 if position.is_any_snp else 0
        contig_stats['best_snps'] += 1 if position.is_best_snp else 0

        outfile.send(position)

            #
            # foo=position.PassedProportionFilter
            # foo.append(ord('\n'))
            # handle.write(foo)

            # handle.write(s.pack(position.called_reference, position.num_A, position.num_C, position.num_G,
            # position.num_T, position.num_N, position.any_snps, position.is_reference_duplicated,
            #                     position.call_str, position.CallWasMade,
            #                     position.PassedDepthFilter, position.PassedProportionFilter, position.Pattern))
    return reference_contig.name, contig_stats


# @profile
def analyze_samples(reference_fasta, reference_dups, sample_analyses):
    """
    Args:
        reference_fasta (Fasta):
        reference_dups (Fasta):
        sample_analyses (list of SampleAnalysis): The sample analysis are grouped by Sample.
    """
    # Group the analyses by sample name in order to collect sample-level statistics.
    # The SampleAnalyses are sorted before grouping because groups are determined by when the key changes.
    sample_groups = tuple(tuple(v) for _, v in itertools.groupby(sorted(sample_analyses), lambda x: x.name))
    with ProcessPoolExecutor() as executor:
        # for result in executor.map(analyze_contig, reference_fasta.contigs, reference_dups.contigs, itertools.repeat(sample_analyses), partial_file_positions):
        for result in executor.map(analyze_contig, reference_fasta.contigs, reference_dups.contigs, itertools.repeat(sample_groups)):
            print(result)
    # whole_genome_stats = Counter()
    # cols = ['reference_length', 'reference_clean', 'reference_duplicated', 'all_called', 'all_passed_coverage', 'all_passed_proportion', 'all_passed_consensus', 'quality_breadth', 'any_snps', 'best_snps']
    # for contig in map(analyze_contig, reference_fasta.contigs, reference_dups.contigs, itertools.repeat(sample_groups),
    #                   partial_file_positions):
    #     print(contig[0], [(x, contig[1][x]) for x in cols])
    #     whole_genome_stats += contig[1]
    # print(whole_genome_stats)

        # with open('partial.nasp', 'br') as handle:
        #     sample_analyses = str(handle.readline(), encoding='utf-8').rstrip().split('\t')
        #     contigs = str(handle.readline(), encoding='utf-8').rstrip().split('\t')
        #     print(sample_analyses)
        #     print(contigs)
        #     for chunk in iter(lambda: handle.read(s.size), b''):
        #         print(s.unpack(chunk))