'''
Created on Nov 19, 2013

@author: dlemmer
'''

import logging
import nasp.convert_external_genome as convert_external_genome
import unittest
import os


class ConvertExternalGenomeTestCase(unittest.TestCase):
    maxDiff = None

    @classmethod
    def setUpClass(cls):
        import re
        import subprocess
        cls.nucmer_path = ""
        cls.deltafilter_path = ""
        try:
            cls.nucmer_path = subprocess.check_output(["which", "nucmer"])
        except subprocess.CalledProcessError:
            try:
                cls.nucmer_path = subprocess.check_output("find ~ -name nucmer", shell=True)
            except subprocess.CalledProcessError as cpe:
                cls.nucmer_path = cpe.output
        cls.nucmer_path = str(cls.nucmer_path, encoding='utf-8')
        match = re.search('^(.*)\n.*$', cls.nucmer_path)
        if match:
            cls.nucmer_path = match.group(1)
        print("found nucmer at: %s" % cls.nucmer_path)
        try:
            cls.deltafilter_path = subprocess.check_output(["which", "delta-filter"])
        except subprocess.CalledProcessError:
            try:
                cls.deltafilter_path = subprocess.check_output("find ~ -name delta-filter", shell=True)
            except subprocess.CalledProcessError as cpe:
                cls.deltafilter_path = cpe.output
        cls.deltafilter_path = str(cls.deltafilter_path, encoding='utf-8')
        match = re.search('^(.*)\n.*$', cls.deltafilter_path)
        if match:
            cls.deltafilter_path = match.group(1)
        print("found delta-filter at: %s" % cls.deltafilter_path)

    def setUp(self):
        self.reference = open("convert_external_genome_reference.fasta", "w")
        self.reference.write(">reference_sequence_1\n")
        self.reference.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatc\n")
        self.reference.write(">reference_sequence_2\n")
        self.reference.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgagctagtcgatcgatttgcatgcagcgcgcgatcattatgcccccat\n")
        self.reference.write(">reference_sequence_3\n")
        self.reference.write("gagcgatcgggatcgagctagctagcacgttcgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatg\n")
        self.reference.close()

        self.external = open("convert_external_genome_external.fasta", "w")
        self.external.write(">external_sequence_1\n")
        self.external.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatc\n")
        self.external.write(">external_sequence_2\n")
        self.external.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgagctagtcgatcgatttgcatgcagcgcgcgatcattatgcccccat\n")
        self.external.write(">external_sequence_3\n")
        self.external.write("gagcgatcgggatcgagctagctagcacgttcgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatg\n")
        self.external.close()

        self.filtered_delta = "external.filtered.delta"

    def tearDown(self):
        os.remove(self.reference.name)
        os.remove(self.external.name)
        if os.path.exists(self.filtered_delta) : os.remove(self.filtered_delta)

    def test_parse_args(self):
        pass #not worthwhile to test

    def test_generate_delta_file(self):
        # Checking if the PermissionError is handled confuses users who expect nucmer is installed correctly.
        # try:
        #     convert_external_genome.generate_delta_file("", "", self.deltafilter_path, "external", self.reference.name, self.external.name)
        # except PermissionError:
        #     self.fail("Unhandled PermissionError. Triggered when nucmer is not installed: self.nucmer_path=\"\"")
        # convert_external_genome.generate_delta_file(self.nucmer_path, "", self.deltafilter_path, "external", self.reference.name, self.external.name)
        convert_external_genome.generate_delta_file(self.nucmer_path, "", self.deltafilter_path, "external", self.reference.name, self.external.name, 'convert_external_genome_external.fasta')

    @unittest.skip("Throws TypeError: missing required arguments instead of asserted OSError. What is this testing?")
    def test_update_genome_from_delta_data(self):
        self.assertRaises(OSError, convert_external_genome._update_genome_from_delta_data, "", "")

    def test_parse_delta_file(self):
        from tests import testdata

        from tempfile import NamedTemporaryFile
        from nasp.nasp_objects import Genome
        franken_genome = Genome()
        external_genome = Genome()
        external_genome.import_fasta_file(testdata.REFERENCE_FASTA)
        convert_external_genome.parse_delta_file(testdata.REFERENCE_DELTA, franken_genome, external_genome)
        with NamedTemporaryFile() as tmpfile:
            franken_genome.write_to_fasta_file(tmpfile.name)

            with open(testdata.REFERENCE_FASTA) as expected, open(tmpfile.name) as actual:
                self.assertEqual(expected.readlines(), actual.readlines())

    def test_parse_delta_line(self):
        pass #tested by parse_delta_file tests

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
