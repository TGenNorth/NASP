# coding: utf-8

__author__ = 'jtravis'

import unittest
from nasp.nasp_objects import Genome


class GenomeTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = Genome()

    @unittest.skip("Covered by GenomeStatus.set_value()")
    def test_set_call(self):
        pass

    @unittest.skip("Covered by GenomeStatus.get_value()")
    def test_get_call(self):
        pass

    # FIXME: it should throw an exception or the prefix should be optional
    # If the prefix is missing, the contig will be appended to the previous contig
    def test__import_fasta_line_missing_prefix(self):
        expected = 'SEQUENCE'
        prefix = 'prefix'
        identifier = '>' + expected
        self.genome._import_fasta_line(identifier, prefix)
        self.assertListEqual([expected], self.genome.get_contigs())

    # FIXME: assertRaises a specific Exception
    def test__import_fasta_line_missing_identifier(self):
        sequence = 'ABCDGHMNRSTUVWXY'
        with self.assertRaises(Exception):
            self.genome._import_fasta_line(sequence)

    # FIXME: assertRaises a specific Exception
    def test__import_fasta_line_missing_contig(self):
        with self.assertRaises(Exception):
            self.genome._import_fasta_line('SEQUENCE1')
            self.genome._import_fasta_line('SEQUENCE2')

    def test__import_fasta_line_identifier_contains_spaces(self):
        identifier = '>prefixName Description'
        prefix = 'prefix'
        expected = 'Name'
        self.genome._import_fasta_line(identifier, prefix)
        self.assertListEqual([expected], self.genome.get_contigs())

    @unittest.skip("Covered by _import_fasta_file tests")
    def test_import_fasta_file(self):
        pass

    def test_reverse_complement(self):
        dna_string = 'ABCDGHMNRSTUVWXYabcdghmnrstuvwxy'
        expected = 'rxwbaasynkdchgvtRXWBAASYNKDCHGVT'
        dna_string2 = 'ABCDGHKNRSTTVWXYabcdghknrsttvwxy'
        self.assertEqual(expected, self.genome.reverse_complement(dna_string))
        self.assertEqual(dna_string2, self.genome.reverse_complement(expected))

    def test_simple_call(self):
        expected = ['A', 'C', 'G', 'T']
        for expect in expected:
            self.assertEqual(expect, self.genome.simple_call(expect.lower()))
        # It should check the base at position one
        self.assertEqual('A', self.genome.simple_call('agctn'))
        # It should replace uracil with thymine
        self.assertEqual('T', self.genome.simple_call('u'))
        # It should replace X with N if not allowed
        self.assertEqual('N', self.genome.simple_call('X', allow_x=False))
        self.assertEqual('X', self.genome.simple_call('X', allow_x=True))
        # It should replace . with N if deletions are not allowed
        self.assertEqual('N', self.genome.simple_call('.', allow_del=False))
        self.assertEqual('.', self.genome.simple_call('.', allow_del=True))
        # It should replace degeneracies with N
        self.assertEqual('N', self.genome.simple_call('d'))

    def test_simple_call_with_empty(self):
        self.assertEqual('N', self.genome.simple_call('', allow_del=False))
        self.assertEqual('.', self.genome.simple_call('', allow_del=True))

    def test_simple_call_with_none(self):
        self.assertEqual('N', self.genome.simple_call(None, allow_del=False))
        self.assertEqual('.', self.genome.simple_call(None, allow_del=True))
