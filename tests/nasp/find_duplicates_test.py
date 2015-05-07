#!/usr/bin/env python3

'''
Created on Oct 18, 2013

@author: dlemmer
'''
import logging
from nasp import find_duplicates
import unittest
import os


class FindDuplicatesTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import re
        import subprocess
        cls.nucmer_path = ""
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

    def setUp(self):
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
        if os.path.exists(self.delta):
            os.remove(self.delta)

    def test_parse_args(self):
        pass #not worthwhile to test

    def test_run_nucmer_on_reference(self):
        try:
            find_duplicates.run_nucmer_on_reference("", self.fasta.name)
        except PermissionError:
            self.fail("Unhandled PermissionError. Triggered when nucmer is not installed: self.nucmer_path=\"\"")
        find_duplicates.run_nucmer_on_reference(self.nucmer_path, self.fasta.name)

    def test_run_nucmer_on_reference_bad_arguments(self):
        self.assertRaises(OSError, find_duplicates.run_nucmer_on_reference, "", "")

    def test_parse_delta_file(self):
        from nasp.nasp_objects import GenomeStatus
        try:
            find_duplicates.run_nucmer_on_reference("", self.fasta.name)
        except PermissionError:
            self.fail("Unhandled PermissionError. Triggered when nucmer is not installed: self.nucmer_path=\"\"")
        find_duplicates.run_nucmer_on_reference(self.nucmer_path, self.fasta.name)
        dups_data = GenomeStatus()
        find_duplicates.parse_delta_file(self.delta, dups_data)
        print(dups_data)
                
    def test_parse_delta_line(self):
        pass #tested by parse_delta_file tests

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
