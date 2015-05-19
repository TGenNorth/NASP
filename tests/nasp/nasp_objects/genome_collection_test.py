# coding: utf-8

__author__ = 'jtravis'

import os
import tempfile
import unittest
from nasp.nasp_objects import ReferenceGenome, FastaGenome, GenomeCollection

from tests import testdata

class GenomeCollectionStatsTestCase(unittest.TestCase):
    def setUp(self):
        reference_path = testdata.REFERENCE_FASTA
        dups_path = testdata.REFERENCE_DUPS
        reference = ReferenceGenome()
        reference.import_fasta_file(reference_path)
        reference.import_dups_file(dups_path)
        self.genome = GenomeCollection()
        self.genome.set_reference(reference)
        fasta = FastaGenome()
        fasta.import_fasta_file(reference_path)
        fasta.set_file_path(reference_path)
        fasta.set_file_type('fasta')
        fasta.set_nickname('test_nickname')
        fasta.add_generators(['genA', 'genB'])
        self.genome.add_genome(fasta)

        self.tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

        # self.matrix_formats = [
        #     {
        #         'dataformat': 'matrix',
        #         'handle': self.tmpfile,
        #         'filter': ''
        #     }
        # ]

    def tearDown(self):
        os.remove(self.tmpfile.name)

    def test_send_to_matrix_handles(self):
        from vcf_to_matrix import write_stats_data
        expected = ''

        with tempfile.TemporaryDirectory() as tmpdir:
            write_stats_data(self.genome, tmpdir)

            os.listdir(tmpdir)

            with open(tmpdir + '/general_stats.tsv') as handle:
                for line in handle:
                    print(line, end='')

class GenomeCollectionMasterMatrixFastaTestCase(unittest.TestCase):
    def setUp(self):
        reference_path = testdata.REFERENCE_FASTA
        dups_path = testdata.REFERENCE_DUPS
        reference = ReferenceGenome()
        reference.import_fasta_file(reference_path)
        reference.import_dups_file(dups_path)
        self.genome = GenomeCollection()
        fasta = FastaGenome()
        fasta.import_fasta_file(reference_path)
        self.genome.add_genome(fasta)
        self.genome.set_reference(reference)

        self.tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

        self.matrix_formats = [
            {
                'dataformat': 'matrix',
                'handle': self.tmpfile,
                'filter': ''
            }
        ]

    def tearDown(self):
        os.remove(self.tmpfile.name)

    def test_send_to_matrix_handles(self):
        expected = """LocusID	Reference	None	#SNPcall	#Indelcall	#Refcall	#CallWasMade	#PassedDepthFilter	#PassedProportionFilter	#A	#C	#G	#T	#Indel	#NXdegen	Contig	Position	InDupRegion	SampleConsensus	CallWasMade	PassedDepthFilter	PassedProportionFilter	Pattern	Pattern#
500WT1_test::1	A	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	1	False	False	'1N'	1
500WT1_test::2	G	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	2	False	False	'1N'	1
500WT1_test::3	C	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	3	False	False	'1N'	1
500WT1_test::4	G	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	4	False	False	'1N'	1
500WT1_test::5	C	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	5	False	False	'1N'	1
500WT1_test::6	C	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	6	False	False	'1N'	1
500WT1_test::7	C	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	7	False	False	'1N'	1
500WT1_test::8	A	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	8	False	False	'1N'	1
500WT1_test::9	A	0	0	0	0/1	1/1	1/1	0	0	0	0	0	1	500WT1_test	9	False	False	'1N'	1"""

        self.genome.send_to_matrix_handles(self.matrix_formats)

        with open(self.tmpfile.name) as handle:
            for expect, observe in zip(expected.split('\n'), handle):
                self.assertEqual(expect + '\n', observe)

class GenomeCollectionBestSnpVcfTestCase(unittest.TestCase):
    def setUp(self):
        reference = ReferenceGenome()
        with tempfile.NamedTemporaryFile(mode='w+') as tmpfile:
            # Position 5 and 6 can be anything, the others must match the VCF REF column or an exception will be thrown.
            tmpfile.write('>500WT1_test\nCCTGGGGA')
            tmpfile.seek(0)
            reference.import_fasta_file(tmpfile.name)
        self.genome = GenomeCollection()
        from nasp.vcf_to_matrix import read_vcf_file
        for genome in read_vcf_file(reference, 10, .9, testdata.GATK_VCF):
            self.genome.add_genome(genome)
        self.genome.set_reference(reference)

        self.tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

        self.matrix_formats = [
            {
                'dataformat': 'vcf',
                'handle': self.tmpfile,
                'filter': ''
            }
        ]

    def tearDown(self):
        os.remove(self.tmpfile.name)

    def test_send_to_matrix_handles(self):
        expected = ''

        self.genome.send_to_matrix_handles(self.matrix_formats)

        with open(self.tmpfile.name) as handle:
            for line in handle:
                print(line)
                # Why is it empty?
                self.assertEqual('fdsfdsfds', line)
