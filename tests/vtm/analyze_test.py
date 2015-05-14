__author__ = 'jtravis'

import unittest
from io import StringIO
from collections import Counter

from nasp.vtm.analyze import GenomeAnalysis, PositionInfo
from nasp.vtm.parse import Position


class SamplePositionsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        vcf = StringIO(initial_value='')
        fasta = StringIO(initial_value='')
        self.analysis = GenomeAnalysis(10, 0.9)

    def tearDown(self):
        pass

    @unittest.skip('Not Implemented')
    def test_foo(self):
        self.analysis.sample_positions('test_contig', self.sample_groups)


class AnalyzePositionDuplicate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analysis = GenomeAnalysis(10, 0.9)

    def test_duplicate_position(self):
        reference_position = Position(
            call='G',
            coverage='-',
            proportion='-'
        )

        dups_position = Position(
            call='1',
            coverage='-',
            proportion='-'
        )

        samples = (
            (
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                ),
            ),
        )
        self.analysis.analyze_position(reference_position, dups_position, samples)

    def test_unique_position(self):
        reference_position = Position(
            call='G',
            coverage='-',
            proportion='-'
        )

        dups_position = Position(
            # Unique position
            call='0',
            coverage='-',
            proportion='-'
        )

        samples = (
            (
                # Single Nucleotide Monomorphism
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                ),
                # Single Nucleotide Polymorphism
                Position(
                    call='A',
                    coverage='-',
                    proportion='-'
                ),
            ),
        )
        self.analysis.analyze_position(reference_position, dups_position, samples)

    def test_no_duplicate_information_position(self):
        """
        Scenario: The reference was not scanned for duplicate positions.
        As a result, it is assumed all positions passed.

        The following fields should be affected:
        - is_reference_duplicated is False
        - called_snp is incremented for snps, but not for monomorphisms
        - is_missing_matrix is True
        - is_all_quality_breadth is True
        """

        reference_position = Position(
            call='G',
            coverage='-',
            proportion='-'
        )

        dups_position = Position(
            # No duplicate information
            call='X',
            coverage='-',
            proportion='-'
        )

        samples = (
            (
                # Single Nucleotide Monomorphism
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                ),
                # Single Nucleotide Polymorphism
                Position(
                    call='A',
                    coverage='-',
                    proportion='-'
                ),
            ),
        )

        expected = PositionInfo(
            is_all_called=True,
            is_reference_clean=True,
            is_reference_duplicated=False,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=True,
            is_all_quality_breadth=True,
            is_best_snp=False,
            all_sample_stats=[
                (
                    Counter(
                        {'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1,
                         'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter(
                        {'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1,
                         'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})),
                [
                    Counter(
                        {'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1,
                         'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter(
                        {'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1,
                         'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter(
                        {'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1,
                         'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})
                ]
            ],
            is_missing_matrix=False,
            called_reference=1,

            called_snp=0,
            passed_coverage_filter=1,
            passed_proportion_filter=1,
            num_A=0,
            num_C=0,
            num_G=1,
            num_T=0,
            num_N=0,
            call_str=['G', 'G'],
            masked_call_str=['G', 'G'],
            CallWasMade='Y',
            PassedDepthFilter='-',
            PassedProportionFilter='-',
            Pattern=['1', '1']
        )

        expected = PositionInfo(
            is_all_called=True,
            is_reference_clean=True,
            is_reference_duplicated=False,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=False,
            is_all_quality_breadth=False,
            is_best_snp=False,
            all_sample_stats=[
                [
                    Counter({'quality_breadth': 1, 'called_reference': 1, 'called_snp': 1, 'passed_coverage_filter': 1,
                             'passed_proportion_filter': 1, 'was_called': 1, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'passed_coverage_filter': 1, 'passed_proportion_filter': 1,
                             'was_called': 1, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0})
                ], [
                    Counter({'quality_breadth': 1, 'called_reference': 1, 'called_snp': 1, 'passed_coverage_filter': 1,
                             'passed_proportion_filter': 1, 'was_called': 1, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'passed_coverage_filter': 1, 'passed_proportion_filter': 1,
                             'was_called': 1, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'called_reference': 1, 'passed_coverage_filter': 1,
                             'passed_proportion_filter': 1, 'was_called': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'called_snp': 1, 'passed_coverage_filter': 1,
                             'passed_proportion_filter': 1, 'was_called': 1, 'called_reference': 0, 'called_degen': 0})
                ]
            ],
            is_missing_matrix=True,
            called_reference=1,
            called_snp=1,
            passed_coverage_filter=2,
            passed_proportion_filter=2,
            num_A=1,
            num_C=0,
            num_G=1,
            num_T=0,
            num_N=0,
            call_str=['G', 'G', 'A'],
            masked_call_str=['G', 'G', 'A'],
            CallWasMade='YY',
            PassedDepthFilter='--',
            PassedProportionFilter='--',
            Pattern=['1', '1', '2']
        )
        self.assertEqual(expected, self.analysis.analyze_position(reference_position, dups_position, samples))

    def test_duplicate_information_position(self):
        """
        Scenario: The reference was scanned for duplicate positions.
        Positions marked as duplicates are not quality positions and should not increment the sample statistics.

        - is_reference_duplicated is True
        - is_missing_matrix is False
        - is_all_quality_breadth is False
        - is_all_passed_consensus is False
        - None of the sample stats are incremented
        """

        reference_position = Position(
            call='G',
            coverage='-',
            proportion='-'
        )

        dups_position = Position(
            # No duplicate information
            call='1',
            coverage='-',
            proportion='-'
        )

        samples = (
            (
                # Single Nucleotide Monomorphism
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                ),
                # Single Nucleotide Polymorphism
                Position(
                    call='A',
                    coverage='-',
                    proportion='-'
                ),
            ),
        )

        expected = PositionInfo(
            is_all_called=True,
            is_reference_clean=True,
            is_reference_duplicated=True,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=False,
            is_all_quality_breadth=False,
            is_best_snp=False,
            all_sample_stats=[
                [
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0})
                ], [
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'was_called': 1, 'passed_proportion_filter': 1,
                             'called_reference': 0, 'quality_breadth': 0, 'called_snp': 0, 'called_degen': 0})
                ]
            ],
            is_missing_matrix=False,
            called_reference=1,
            called_snp=0,
            passed_coverage_filter=2,
            passed_proportion_filter=2,
            num_A=1,
            num_C=0,
            num_G=1,
            num_T=0,
            num_N=0,
            call_str=['G', 'G', 'A'],
            masked_call_str=['G', 'G', 'A'],
            CallWasMade='YY',
            PassedDepthFilter='--',
            PassedProportionFilter='--',
            Pattern=['1', '1', '2']
        )
        self.assertEqual(expected, self.analysis.analyze_position(reference_position, dups_position, samples))


class AnalyzePositionTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.analysis = GenomeAnalysis(10, 0.9)

    def tearDown(self):
        pass

    @unittest.skip('Not Implemented')
    def test_it_filters_position_by_coverage(self):
        """
        Change the coverage to check for a hardcoded value.
        """
        self.analysis.analyze_position((), (), ())

    @unittest.skip('Not Implemented')
    def test_it_filters_position_by_proportion(self):
        # Change the proportion to check for a hardcoded value.
        self.analysis.analyze_position((), (), ())

    def test_fasta_position(self):
        reference_position = Position(
            call='G',
            coverage='-',
            proportion='-'
        )

        dups_position = Position(
            call='0',
            coverage='-',
            proportion='-'
        )

        samples = (
            (
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                ),
                Position(
                    call='G',
                    coverage='-',
                    proportion='-'
                )
            ),
        )

        expected = PositionInfo(
            is_all_called=True,
            is_reference_clean=True,
            is_reference_duplicated=False,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=True,
            is_all_quality_breadth=True,
            is_best_snp=False,
            all_sample_stats=[
                [
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0})
                ], [
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'called_reference': 1, 'passed_coverage_filter': 1, 'quality_breadth': 1, 'was_called': 1,
                             'passed_proportion_filter': 1, 'called_snp': 0, 'called_degen': 0})
                ]
            ],
            is_missing_matrix=False,
            called_reference=2,
            called_snp=0,
            passed_coverage_filter=2,
            passed_proportion_filter=2,
            num_A=0,
            num_C=0,
            num_G=2,
            num_T=0,
            num_N=0,
            call_str=['G', 'G', 'G'],
            masked_call_str=['G', 'G', 'G'],
            CallWasMade='YY',
            PassedDepthFilter='--',
            PassedProportionFilter='--',
            Pattern=['1', '1', '1']
        )

        self.assertEqual(expected, self.analysis.analyze_position(reference_position, dups_position, samples))

        # # Infinite loop
        # def test_sample_positions_empty_sample_analyses(self):
        # raise NotImplementedError
        # contig_name = 'foo'
        # sample_positions = tuple(tuple())
        # positions = analyze.sample_positions(contig_name, sample_positions)
        # for position in positions:
        # print(position)


from tests import testdata
from nasp.vtm.parse import Fasta


class AnalyzeContigTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        fasta = Fasta(testdata.PARSE_FASTA, 'test_data', None, False)
        ref_fasta = Fasta(testdata.PARSE_FASTA, 'test_data', None, True)

        self.analysis = GenomeAnalysis(10, .9)
        self.coroutine_fn = lambda *arg: (self.CoroutineMock(), self.CoroutineMock())
        self.dups_contig = fasta
        self.reference_contig = ref_fasta.get_contig('ContigWithFilePositionOffset')
        self.sample_groups = (
            (
                fasta,
                fasta
            ), (
                fasta,
            ), (
                fasta,
                fasta
            )
        )

    def tearDown(self):
        pass

    # def test_it_calls_the_coroutine_function_once(self):
    # coroutine_fn_call_count = 0
    #
    # def coroutine_fn():
    # return lambda *args: coroutine_fn_call_count += 1
    #
    # self.analysis.analyze_contig(coroutine_fn, sample_groups, dups_contig, reference_contig)

    class CoroutineMock:
        def send(self, arg):
            pass

    def test_it_collects_analysis_and_overall_contig_stats(self):
        expected_analysis_stats = [
            [
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}),
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0})
            ],
            [
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}), Counter(
                {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                 'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}), Counter(
                {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                 'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}), Counter(
                {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                 'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0})
            ], [
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}), Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}), Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0})
            ], [
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}),
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}),
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0}),
                Counter(
                    {'quality_breadth': 8, 'passed_proportion_filter': 8, 'called_reference': 8, 'was_called': 8,
                     'passed_coverage_filter': 8, 'called_snp': 0, 'called_degen': 0})
            ]
        ]
        expected_contig_stats = Counter(
            {'reference_clean': 8, 'all_passed_proportion': 8, 'quality_breadth': 8,
             'Contig': 'ContigWithFilePositionOffset', 'reference_length': 8, 'best_snps': 0,
             'all_passed_coverage': 8, 'reference_duplicated': 0, 'all_passed_consensus': 8, 'all_called': 8,
             'any_snps': 0}
        )
        observed_analysis_stats, observed_contig_stats = self.analysis.analyze_contig(self.coroutine_fn,
                                                                                      self.sample_groups,
                                                                                      self.dups_contig,
                                                                                      self.reference_contig)

        self.assertEqual(expected_analysis_stats, observed_analysis_stats)
        self.assertEqual(expected_contig_stats, observed_contig_stats)
