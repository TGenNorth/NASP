# coding: utf-8

__author__ = 'jtravis'

import os
import tempfile
import unittest
from nasp.nasp_objects import ReferenceGenome, FastaGenome, GenomeCollection

from tests import testdata

class GenomeCollectionStatsTestCase(unittest.TestCase):
    maxDiff = None

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
        # Statistics are gathered when the matrices are created
        self.genome.write_to_matrices({})

        self.tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

    def tearDown(self):
        os.remove(self.tmpfile.name)

    def test_send_to_matrix_handles(self):
        from vcf_to_matrix import write_stats_data
        expected = (
            "Contig	reference_length	reference_clean	reference_clean (%)	reference_duplicated	reference_duplicated (%)	all_called	all_called (%)	all_passed_coverage	all_passed_coverage (%)	all_passed_proportion	all_passed_proportion (%)	all_passed_consensus	all_passed_consensus (%)	quality_breadth	quality_breadth (%)	any_snps	any_snps (%)	best_snps	best_snps (%)	\n"
	        "	stat descriptions go here\n"
            "Whole Genome	3977	3977	100.00%	0	0.00%	3977	100.00%	3977	100.00%	3977	100.00%	3977	100.00%	3977	100.00%	0	0.00%	0	0.00%	\n"
            "500WT1_test	3977	3977	100.00%	0	0.00%	3977	100.00%	3977	100.00%	3977	100.00%	3977	100.00%	3977	100.00%	0	0.00%	0	0.00%	\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            write_stats_data(self.genome, tmpdir)

            os.listdir(tmpdir)

            with open(tmpdir + '/general_stats.tsv') as handle:
                self.assertEqual(expected, ''.join(handle.readlines()))

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
500WT1_test::1	A	0	0	1	1/1	1/1	1/1	1	0	0	0	0	0	500WT1_test	1	False	True	'11'	1
500WT1_test::2	G	0	0	1	1/1	1/1	1/1	0	0	1	0	0	0	500WT1_test	2	False	True	'11'	1
500WT1_test::3	C	0	0	1	1/1	1/1	1/1	0	1	0	0	0	0	500WT1_test	3	False	True	'11'	1
500WT1_test::4	G	0	0	1	1/1	1/1	1/1	0	0	1	0	0	0	500WT1_test	4	False	True	'11'	1
500WT1_test::5	C	0	0	1	1/1	1/1	1/1	0	1	0	0	0	0	500WT1_test	5	False	True	'11'	1
500WT1_test::6	C	0	0	1	1/1	1/1	1/1	0	1	0	0	0	0	500WT1_test	6	False	True	'11'	1
500WT1_test::7	C	0	0	1	1/1	1/1	1/1	0	1	0	0	0	0	500WT1_test	7	False	True	'11'	1
500WT1_test::8	A	0	0	1	1/1	1/1	1/1	1	0	0	0	0	0	500WT1_test	8	False	True	'11'	1
500WT1_test::9	A	0	0	1	1/1	1/1	1/1	1	0	0	0	0	0	500WT1_test	9	False	True	'11'	1
500WT1_test::10	T	0	0	1	1/1	1/1	1/1	0	0	0	1	0	0	500WT1_test	10	False	True	'11'	1"""
        # ...
        # There might be more data, but 10 samples is enough

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
            self.assertEqual(expected, handle.readlines())
