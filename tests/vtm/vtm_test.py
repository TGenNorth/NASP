import itertools
from nasp.vtm.matrix_DTO import NaspFile

__author__ = 'jtravis'

import unittest
from nasp.vtm import vtm
from nasp.vtm.parse import Fasta
from tests import testdata


# TODO: What if fastas and vcfs are mixed in the dto?
# What if no samples are present?


# class ParseArgsTestCase(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         pass
#
#     def setUp(self):
#         pass
#
#     def tearDown(self):
#         pass


class IndexContigsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.frankenfastas = []
        self.vcfs = []
        for a, b, c in itertools.permutations(('a', 'b', 'c')):
            self.frankenfastas.append(
                NaspFile(
                    name=a,
                    path=testdata.PARSE_FASTA,
                    aligner=b,
                    snpcaller=c
                )
            )

            self.vcfs.append(
                NaspFile(
                    name=a,
                    path=testdata.GATK_VCF,
                    aligner=b,
                    snpcaller=c
                )
            )

    def tearDown(self):
        pass

    def test__index_contig(self):
        (reference, dups, sample_groups) = vtm._index_contigs(
            testdata.REFERENCE_FASTA,
            testdata.REFERENCE_DUPS,
            self.frankenfastas,
            self.vcfs
        )

        self.assertIsInstance(reference, Fasta)
        self.assertIsInstance(dups, Fasta)

        # The individual sample analyses should be grouped by sample name; in this case a, b, and c
        # In this test, each sample should have 4 analyses
        for expected_sample_name, sample in zip(('a', 'b', 'c'), sample_groups):
            self.assertEqual(4, len(sample))
            for analysis in sample:
                self.assertEqual(expected_sample_name, analysis.name)


    # def test__index_contig_no_samples(self):
    #     (reference, dups, sample_groups) = vtm._index_contigs(data.REFERENCE_FASTA, data.REFERENCE_DUPS, (), ())
    #
    #     self.assertIsInstance(reference, Fasta)
    #     self.assertIsInstance(dups, Fasta)

    def test__index_contig_sorts_samples_by_identifier(self):
        expected = (
            'a::b',
            'a::b,c',
            'a::c',
            'a::c,b',
            'b::a',
            'b::a,c',
            'b::c',
            'b::c,a',
            'c::a',
            'c::a,b',
            'c::b',
            'c::b,a'
        )

        (reference, dups, sample_groups) = vtm._index_contigs(
            testdata.REFERENCE_FASTA,
            testdata.REFERENCE_DUPS,
            self.frankenfastas,
            self.vcfs
        )

        for expect, observe in itertools.zip_longest(expected, itertools.chain.from_iterable(sample_groups)):
            self.assertEqual(expect, observe.identifier)
