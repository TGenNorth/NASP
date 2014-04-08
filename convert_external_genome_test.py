'''
Created on Nov 19, 2013

@author: dlemmer
'''

import convert_external_genome
import unittest
import os

class ConvertExternalGenomeTestCase(unittest.TestCase):

    def setUp(self):
        import re
        import subprocess
        self.nucmer_path = ""
        self.deltafilter_path = ""
        try:
            self.nucmer_path = subprocess.check_output(["which", "nucmer"])
        except subprocess.CalledProcessError:
            try:
                self.nucmer_path = subprocess.check_output("find ~ -name nucmer", shell=True)
            except subprocess.CalledProcessError as cpe:
                self.nucmer_path = cpe.output
        match = re.search('^(.*)\n.*$', self.nucmer_path)
        if match: self.nucmer_path = match.group(1)
        print("found nucmer at: %s" % self.nucmer_path)
        try:
            self.deltafilter_path = subprocess.check_output(["which", "delta-filter"])
        except subprocess.CalledProcessError:
            try:
                self.deltafilter_path = subprocess.check_output("find ~ -name delta-filter", shell=True)
            except subprocess.CalledProcessError as cpe:
                self.deltafilter_path = cpe.output
        match = re.search('^(.*)\n.*$', self.deltafilter_path)
        if match: self.deltafilter_path = match.group(1)
        print("found delta-filter at: %s" % self.deltafilter_path)
        
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
        convert_external_genome.generate_delta_file(self.nucmer_path, "", self.deltafilter_path, "external", self.reference.name, self.external.name)

    def test_update_genome_from_delta_data(self):
        self.assertRaises(OSError, convert_external_genome._update_genome_from_delta_data, "", "")

    def test_parse_delta_file(self):
        from nasp_objects import GenomeStatus
        convert_external_genome.generate_delta_file(self.nucmer_path, "", self.deltafilter_path, "external", self.reference.name, self.external.name)
        dups_data = GenomeStatus()
        convert_external_genome.parse_delta_file(self.delta, dups_data)
        print(dups_data)
                
    def test_parse_delta_line(self):  
        pass #tested by parse_delta_file tests

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()