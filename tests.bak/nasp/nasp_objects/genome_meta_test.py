# coding: utf-8

__author__ = 'jtravis'

import unittest
from nasp.nasp_objects import GenomeMeta


class GenomeMetaTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = GenomeMeta()

    # TODO: Should the side effect of setting an undefined nickname be tested?
    def test_set_file_path(self):
        expected = 'foo'
        self.genome.set_file_path(expected)
        self.assertEqual(expected, self.genome.file_path())

    def test_set_file_type(self):
        expected = 'foo'
        self.genome.set_file_type(expected)
        self.assertEqual(expected, self.genome.file_type())

    def test_set_nickname(self):
        expected = 'foo'
        self.genome.set_nickname(expected)
        self.assertEqual(expected, self.genome.nickname())

    @unittest.skip("Not Implemented")
    def test_add_generators(self):
        pass

    def test_file_path(self):
        expected = 'foo'
        self.assertIsNone(self.genome.file_path())
        self.genome.set_file_path(expected)
        self.assertEqual(expected, self.genome.file_path())

    def test_file_type(self):
        expected = 'foo'
        self.assertIsNone(self.genome.file_type())
        self.genome.set_file_type(expected)
        self.assertEqual(expected, self.genome.file_type())

    def test_nickname(self):
        expected = 'foo'
        self.assertIsNone(self.genome.nickname())
        self.genome.set_nickname(expected)
        self.assertEqual(expected, self.genome.nickname())
        pass

    def test_identifier(self):
        generators = ['bar', 'baz', 'quox']
        nickname = 'foo'
        self.assertIsNone(self.genome.identifier())
        self.genome.set_nickname(nickname)
        self.assertEqual(nickname, self.genome.identifier())
        self.genome.add_generators(generators)
        self.assertEqual(nickname + '::' + ','.join(generators), self.genome.identifier())

    def test_generate_nickname_from_filename(self):
        expected = 'foo'
        paths = [
            '',
            '/',
            '/bar/baz/quox/'
            'bar/baz/quox/'
        ]
        extensions = [
            '.fRaNkEnFaStA',
            '.FaStA',
            '.fAs',
            '.VcF'
        ]
        for extension in extensions:
            filename = expected + extension
            self.assertEqual(expected, self.genome.generate_nickname_from_filename(filename))

        for path in paths:
            filename = path + expected + '.vcf'
            self.assertEqual(expected, self.genome.generate_nickname_from_filename(filename))

    @unittest.skip("Duplicate of Genome.reverse_complement")
    def test_reverse_complement(dna_string):
        pass