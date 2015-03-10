"""
The analyze module takes a collection of SampleAnalyses.
"""

__author__ = 'jtravis'

from collections import namedtuple, Counter
from concurrent.futures import ProcessPoolExecutor
import itertools
import os

from nasp2.write_matrix import write_master_matrix, write_bestsnp_matrix, write_missingdata_matrix, write_general_stats, write_sample_stats

# TODO: Remove constant.
OUTPUT_DIR = './'


# TODO: Remove profile. It is a dummy wrapper to leave @profile decorators in place when not profiling.
# def profile(fn):
# return fn

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
    # True if all sample calls differ from the reference and passed all the filters (coverage, proportion, consensus).
    'is_best_snp',

    'all_sample_stats',

    'is_missing_matrix',

    # Count of sample calls that match the reference.
    'called_reference',
    'called_snp',
    'passed_coverage_filter',
    'passed_proportion_filter',
    'num_A',
    'num_C',
    'num_G',
    'num_T',
    'num_N',

    #
    'call_str',
    'masked_call_str',
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
    the files to be read in unison as if they were a single file.

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

        The sample positions generator yields a tuple of Position tuples:
        ( (SampleInfo, SampleInfo,), (SampleInfo, SampleInfo), )

        +------------------------------------------------------------------------------------+
        |              Sample Group                |              Sample Group               |
        +==========================================+=========================================+
        | Position           | Position            | Position           | Position           |
        +------------------------------------------+-----------------------------------------+

        The order or Positions matching the order of SampleAnalyses is how we know which Position
        belongs to which SampleAnalysis.
    """
    # sample_groups is converted to contig_groups by calling get_contig() on each SampleAnalysis.
    contig_groups = tuple(
        tuple(map(lambda sample_analysis: sample_analysis.get_contig(contig_name).positions, sample_group)) \
        for sample_group in sample_groups
    )

    while True:
        # The inner tuple is all the sample analyses for a single sample.
        # The outer tuple is all the sample analyses.
        # ( (SampleInfo, SampleInfo,), (SampleInfo, SampleInfo), )
        yield tuple(tuple(map(next, contig_group)) for contig_group in contig_groups)


def _is_pass_filter(value, threshold):
    """
    _filter determines if a position passed/failed the coverage/proportion filter.

    Args:
        value (int or str): coverage/proportion value from the file parser.
        threshold (int): The minimum value that will pass.

    Returns:
        (bool, str) True if it passed and a character representing the reason it passed/failed.
    """
    if value is '-':
        # Cannot determine due to no value specified. Assume it passed.
        # Always the case for Fasta files, sometimes the case for VCF.
        return True, '-'
    if value is '?':
        # Missing VCF position.
        return False, '?'
    if value >= threshold:
        # PASS - value is greater than threshold
        return True, 'Y'
    # FAIL - value is less than threshold
    return False, 'N'


# @profile
def analyze_position(reference_position, dups_position, samples):
    """
    analyze_position compares all the SampleAnalyses at a single position.

    Args:
        reference_position (Position): A single position from the reference genome.
        dups_position (Position): A single position from the duplicates file which indicates if the reference is in a duplicate region.
        samples (tuple of tuples of Position): A single position across all SampleAnalyses grouped by sample name.

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
    is_all_sample_consensus = True
    # General Stats quality breadth passes if all the following conditions for the position are true:
    # - Reference called and not duplicated.
    # - All analyses called and passed all filters.
    # - All samples have consensus within their aligner/snpcaller combinations.
    is_all_quality_breadth = True

    # A position is included in the missing data matrix if at least one analysis passes quality breadth and is a SNP.
    is_missing_matrix = False

    # position_stats is a counter across all SampleAnalyses at the current position. It used in the Matrix files.
    position_stats = Counter({
        # TODO: verify 'N' represents NoCall, AnyCall, Deletion, etc.
        # Tally of SampleAnalysis A/C/G/T calls where anything else is called N.
        # 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
        #
        # Tally of SampleAnalyses with any base call. Excludes None (X) and Any (N) calls.
        # 'was_called': 0,
        #
        # Tally of SampleAnalyses with coverage/proportion at least as high as the threshold.
        # 'passed_coverage_filter': 0, 'passed_proportion_filter': 0
        #
        # Tally of SampleAnalyses where the call was a degenerate, matched, or differed from the reference call.
        # 'called_degen': 0, 'called_reference': 0, 'called_snp': 0
        #
        # The character at each position of the following filter strings indicates whether a corresponding SampleAnalysis
        # passed the filter:
        'CallWasMade': '',
        'PassedDepthFilter': '',
        'PassedProportionFilter': '',
        'Pattern': ''
    })

    # TODO: replace bytearray with a list.
    # call_str is the base call for each SampleAnalysis starting with the reference.
    call_str = bytearray(reference_position.call, encoding='utf-8')
    # masked_call_str is the same as the call string except calls that did not pass the coverage/proportion filter are
    # masked out with 'N' indicating they did not pass.
    masked_call_str = [reference_position.call]

    # all_sample_stats is a list of sample_stats. It represents all the SampleAnalysis stats grouped by sample name.
    all_sample_stats = []
    for sample in samples:
        # sample_stats is composed of analysis_stats dictionaries for each aligner/snpcaller used on the sample.
        # The first two entries are special in that they summarize the analysis_stats. They are using a tripwire approach
        # where the flag will toggle if the assumption fails.
        sample_stats = [
            # any - True if true in any of the analysis_stats.
            Counter({
                'was_called': 0,
                'passed_coverage_filter': 0,
                'passed_proportion_filter': 0,
                'quality_breadth': 0,
                'called_reference': 0,
                'called_snp': 0,
                'called_degen': 0
            }),
            # all - True if true in all of the analysis_stats.
            Counter({
                'was_called': 1,
                'passed_coverage_filter': 1,
                'passed_proportion_filter': 1,
                'quality_breadth': 1,
                'called_reference': 1,
                'called_snp': 1,
                'called_degen': 1
            })
        ]

        # prev_analysis_call is used to detect sample consensus.
        prev_analysis_call = None
        for analysis in sample:
            # analysis_stats are unique to each SampleAnalysis. It is used in the Sample Statistics file.
            analysis_stats = Counter({
                'was_called': 0,
                'passed_coverage_filter': 0,
                'passed_proportion_filter': 0,
                'quality_breadth': 0,
                'called_reference': 0,
                'called_snp': 0,
                'called_degen': 0
            })

            # TODO: Must also be called (not X/N) and pass coverage/proportion filters
            # A sample has position consensus if all of its SampleAnalysis combinations have the same call.
            if prev_analysis_call is None:
                prev_analysis_call = analysis.simple_call
            elif prev_analysis_call != analysis.simple_call:
                is_all_sample_consensus = False
                is_all_quality_breadth = False

            # X means no value, N means any value.
            if analysis.call in ['X', 'N']:
                position_stats['CallWasMade'] += 'N'
                is_all_called = False
                is_all_quality_breadth = False
            else:
                analysis_stats['was_called'] = 1
                position_stats['was_called'] += 1
                position_stats['CallWasMade'] += 'Y'

            # FIXME: coverage/proportion may be "PASS" if the value is not deprecated in the VCF contig parser.
            # get_proportion/get_coverage


            is_pass_coverage, pass_str = _is_pass_filter(analysis.coverage, coverage_threshold)
            position_stats['PassedDepthFilter'] += pass_str
            if is_pass_coverage:
                analysis_stats['passed_coverage_filter'] = 1
                position_stats['passed_coverage_filter'] += 1
            else:
                is_all_passed_coverage = False

            is_pass_proportion, pass_str = _is_pass_filter(analysis.proportion, proportion_threshold)
            position_stats['PassedProportionFilter'] += pass_str
            if is_pass_proportion:
                analysis_stats['passed_proportion_filter'] = 1
                position_stats['passed_proportion_filter'] += 1
            else:
                is_all_passed_proportion = False

            call_str.append(ord(analysis.call))

            #
            if not analysis_stats['was_called'] or is_pass_coverage and is_pass_proportion:
                masked_call_str.append(analysis.call)
            else:
                masked_call_str.append('N')

            # Count the number of A/C/G/T calls. Any other value is N. See Position.simple_call().
            position_stats[analysis.simple_call] += 1

            # TODO: Clarify the purpose of this if statement.
            # Only count significant measurements.
            if analysis_stats['passed_coverage_filter'] and analysis_stats['passed_proportion_filter'] and \
                            reference_position.simple_call != 'N' and dups_position.call != '1':

                analysis_stats['quality_breadth'] = 1

                # TODO: Should the first if be simplecall?
                if analysis.call in ['X', 'N']:
                    position_stats['called_degen'] += 1
                    analysis_stats['called_degen'] = 1
                elif analysis.call == reference_position.call:
                    position_stats['called_reference'] += 1
                    analysis_stats['called_reference'] = 1
                else:
                    position_stats['called_snp'] += 1
                    analysis_stats['called_snp'] = 1
                    is_missing_matrix = True
            else:
                is_all_quality_breadth = False

            sample_stats.append(analysis_stats)


        # NOTE: There could be an index error due to [2:] if there is not at least one sample analysis in sample_stats.
        # Summarize the statistics for each sample where any or all of its SampleAnalysis combinations passed.
        # It is using a tripwire approach where the if statement toggles the flag if the assumption fails.
        # - any assumes all the stats failed
        # - all assumes all the stats passed
        for analysis_stats in sample_stats[2:]:
            for key, value in analysis_stats.items():
                if value:
                    # any - The stat passed in at least one analysis_stat.
                    sample_stats[0][key] = 1
                else:
                    # all - The stat failed in at least one analysis_stat.
                    sample_stats[1][key] = 0

        all_sample_stats.append(sample_stats)

    return PositionInfo(
        # General Stats
        is_all_called=is_all_called,
        is_reference_clean=is_reference_clean,
        is_reference_duplicated=is_reference_duplicated,
        is_all_passed_coverage=is_all_passed_coverage,
        is_all_passed_proportion=is_all_passed_proportion,
        is_all_passed_consensus=is_all_sample_consensus,
        is_all_quality_breadth=is_all_quality_breadth,
        is_best_snp=is_all_quality_breadth and position_stats['called_snp'] > 0,

        # NOTE: Would it increase performance if this were a tuple of tuples and we avoided using append?
        # Sample Stats - A list of list of Counters representing the stats for each analysis file grouped by sample.
        all_sample_stats=all_sample_stats,

        # Missing Data Matrix condition - at least one SampleAnalysis passes quality_breadth and is a SNP.
        is_missing_matrix=is_missing_matrix,

        # NASP Master Matrix
        # Counters
        called_reference=position_stats['called_reference'],
        called_snp=position_stats['called_snp'],
        passed_coverage_filter=position_stats['passed_coverage_filter'],
        passed_proportion_filter=position_stats['passed_proportion_filter'],
        num_A=position_stats['A'],
        num_C=position_stats['C'],
        num_G=position_stats['G'],
        num_T=position_stats['T'],
        num_N=position_stats['N'],
        # Strings
        call_str=call_str,
        masked_call_str=masked_call_str,
        CallWasMade=bytearray(position_stats['CallWasMade'], encoding='utf-8'),
        PassedDepthFilter=bytearray(position_stats['PassedDepthFilter'], encoding='utf-8'),
        PassedProportionFilter=bytearray(position_stats['PassedProportionFilter'], encoding='utf-8'),
        Pattern=bytearray(position_stats['Pattern'], encoding='utf-8')
    )


# @profile
def analyze_contig(reference_contig, dups_contig, sample_groups):
    """
    analyze_contig expects reference_contig, dups_contig, and sample_analyses to represent the same contig from
    different sources.

    Args:
        reference_contig (FastaContig):
        dups_contig (FastaContig):
        sample_groups (tuple of SampleAnalysis tuples): SampleAnalyses grouped by sample name.

    Returns:
        (list of lists of dict, collections.Counter): SampleAnalysis stats grouped by sample name and contig stats.

    Note:
        Side-effect:
    """
    # sample stats is a list of dict lists representing the SampleAnalysis stats in the same order as sample_groups.
    sample_stats = None
    contig_stats = Counter({'Contig': reference_contig.name})

    # FIXME: Must the identifiers be recalculated for every contig?
    identifiers = tuple(analysis.identifier for sample in sample_groups for analysis in sample)

    # Initialize outfile coroutines. This will create the file, write the header, and wait for each line of data.
    master_matrix = write_master_matrix(os.path.join(OUTPUT_DIR, reference_contig.name + '_master.tsv'), reference_contig.name,
                                        identifiers)
    bestsnp_matrix = write_bestsnp_matrix(os.path.join(OUTPUT_DIR, reference_contig.name + '_best.tsv'), reference_contig.name,
                                        sample_groups)
    missing_matrix = write_missingdata_matrix(os.path.join(OUTPUT_DIR, reference_contig.name + '_missing.tsv'), reference_contig.name,
                                        identifiers)
    master_matrix.send(None)
    bestsnp_matrix.send(None)
    missing_matrix.send(None)

    for position in map(analyze_position, reference_contig.positions, dups_contig.positions,
                        sample_positions(reference_contig.name, sample_groups)):
        contig_stats['reference_length'] += 1
        contig_stats['reference_clean'] += 1 if position.is_reference_clean else 0
        contig_stats['reference_duplicated'] += 1 if position.is_reference_duplicated else 0
        contig_stats['all_called'] += 1 if position.is_all_called else 0
        contig_stats['all_passed_coverage'] += 1 if position.is_all_passed_coverage else 0
        contig_stats['all_passed_proportion'] += 1 if position.is_all_passed_proportion else 0
        contig_stats['all_passed_consensus'] += 1 if position.is_all_passed_consensus else 0
        contig_stats['quality_breadth'] += 1 if position.is_all_quality_breadth else 0
        contig_stats['any_snps'] += 1 if position.called_snp > 0 else 0
        contig_stats['best_snps'] += 1 if position.is_best_snp else 0

        # Sum the SampleAnalysis stats preserving sample_groups order.
        if sample_stats is None:
            sample_stats = position.all_sample_stats
        else:
            for i, sample in enumerate(position.all_sample_stats):
                for j, analysis in enumerate(sample):
                    sample_stats[i][j].update(analysis)

        # Write position to Contig Master Matrix.
        master_matrix.send(position)
        bestsnp_matrix.send(position)
        missing_matrix.send(position)

    return sample_stats, contig_stats


# TODO: Should analyze_samples return the contig_stats generator so that a unit test can check the return value and let
# another function actually write the stats file?
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
        contig_stats = []
        # NOTE: if the loop does not run, None will be passed to write_general_stats.
        sample_stats = None

        # for sample_stat, contig_stat in executor.map(analyze_contig, reference_fasta.contigs, reference_dups.contigs,
        for sample_stat, contig_stat in map(analyze_contig, reference_fasta.contigs, reference_dups.contigs,
                                                     itertools.repeat(sample_groups)):
            # TODO: use the contig name from contig stat to begin reassembling the matrices.
            contig_stats.append(contig_stat)

            # Sum the SampleAnalysis stats preserving the sample_groups order.
            if sample_stats is None:
                sample_stats = sample_stat
            else:
                for i, sample in enumerate(sample_stats):
                    for j, analysis in enumerate(sample):
                        sample_stats[i][j].update(analysis)

        write_sample_stats(os.path.join(OUTPUT_DIR, 'sample_stats.tsv'), sample_stats, sample_groups)

        write_general_stats(os.path.join(OUTPUT_DIR, 'general_stats.tsv'), contig_stats)