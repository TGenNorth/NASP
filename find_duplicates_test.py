#!/usr/bin/env python3

'''
Created on Oct 18, 2013

@author: dlemmer
'''
import find_duplicates
import unittest
import os

class FindDuplicatesTestCase(unittest.TestCase):

    def setUp(self):
        import re
        import subprocess
        self.nucmer_path = ""
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
        
        self.fasta = open("find_dups_test.fasta", "w")
        self.fasta.write(">test_sequence_1\n")
        self.fasta.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatc\n")
        self.fasta.write(">test_sequence_2\n")
        self.fasta.write("ctgcgatcgggatcgagctagctagcacgtacgtacgtacgtacgtacgtacgagctagtcgatcgatttgcatgcagcgcgcgatcattatgcccccat\n")
        self.fasta.write(">test_sequence_3\n")
        self.fasta.write("gagcgatcgggatcgagctagctagcacgttcgtacgtacgtacgtacgtacgtacgtacgatcagctagcagcatcgatcgatcgacgatacgatcatg\n")
        self.fasta.close()
        
        self.delta = "reference.delta"

    def tearDown(self):
        os.remove(self.fasta.name)
        if os.path.exists(self.delta) : os.remove(self.delta)

    def test_parse_args(self):
        pass #not worthwhile to test

    def test_run_nucmer_on_reference(self):
        find_duplicates.run_nucmer_on_reference(self.nucmer_path, self.fasta.name)

    def test_run_nucmer_on_reference_bad_arguments(self):
        self.assertRaises(OSError, find_duplicates.run_nucmer_on_reference, "", "")

    def test_parse_delta_file(self):
        from nasp_objects import GenomeStatus
        find_duplicates.run_nucmer_on_reference(self.nucmer_path, self.fasta.name)
        dups_data = GenomeStatus()
        find_duplicates.parse_delta_file(self.delta, dups_data)
        print(dups_data)
                
    def test_parse_delta_line(self):
        pass #tested by parse_delta_file tests

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()