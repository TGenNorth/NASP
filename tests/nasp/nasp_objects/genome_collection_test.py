# coding: utf-8

__author__ = 'jtravis'

import tempfile
import unittest
from nasp.nasp_objects import ReferenceGenome, FastaGenome, GenomeCollection

from tests import testdata

class GenomeCollectionTestCase(unittest.TestCase):
    def setUp(self):
        reference_path = testdata.REFERENCE_FASTA
        dups_path = testdata.REFERENCE_DUPS
        reference = ReferenceGenome()
        reference.import_fasta_file(reference_path)
        reference.import_dups_file(dups_path)
        self.genome = GenomeCollection()
        fasta = FastaGenome()
        self.genome.add_genome(fasta)
        self.genome.set_reference(reference)
        self.matrix_formats = [
            {
                'dataformat': 'matrix',
                'handle': tempfile.TemporaryFile(mode='w'),
                'filter': ''
            }
        ]

    def tearDown(self):
        for matrix in self.matrix_formats:
            matrix['handle'].close()

    def test_send_to_matrix_handles(self):
        self.genome.send_to_matrix_handles(self.matrix_formats)
