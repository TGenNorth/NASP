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
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})),
                [
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})
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
        self.assertEqual(expected, self.analysis.analyze_position(reference_position, dups_position, samples))

    def test_duplicate_information_position(self):
        """
        Scenario: The reference was scanned for duplicate positions.
        Positions marked as duplicates should not

        TODO:
        - is_reference_duplicated is True
        - is_missing_matrix is False
        - is_all_quality_breadth is False
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
            is_reference_duplicated=False,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=True,
            is_all_quality_breadth=True,
            is_best_snp=False,
            all_sample_stats=[
                (
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})),
                [
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'quality_breadth': 1, 'was_called': 1, 'called_reference': 1, 'passed_proportion_filter': 1, 'passed_coverage_filter': 1, 'called_snp': 0, 'called_degen': 0})
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
            is_reference_duplicated=True,
            is_all_passed_coverage=True,
            is_all_passed_proportion=True,
            is_all_passed_consensus=True,
            is_all_quality_breadth=False,
            is_best_snp=False,
            all_sample_stats=[
                [
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0})
                ],
                [
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0}),
                    Counter({'passed_coverage_filter': 1, 'passed_proportion_filter': 1, 'was_called': 1, 'quality_breadth': 0, 'called_reference': 0, 'called_snp': 0, 'called_degen': 0})
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
        #     positions = analyze.sample_positions(contig_name, sample_positions)
        #     for position in positions:
        #         print(position)


class AnalyzeContigTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

        # # Infinite loop
        # def test_sample_positions_empty_sample_analyses(self):
        # raise NotImplementedError
        # contig_name = 'foo'
        # sample_positions = tuple(tuple())
        #     positions = analyze.sample_positions(contig_name, sample_positions)
        #     for position in positions:
        #         print(position)
