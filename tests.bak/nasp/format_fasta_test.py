__author__ = 'jtravis'

import unittest

class FormatFastaTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_format_fasta_partitions_into_80_character_lines(self):
        from tempfile import NamedTemporaryFile
        from tests import testdata
        from nasp.format_fasta import format_fasta

        expected_max_line_length = 80 + 1 # +1 for the newline character


        with NamedTemporaryFile() as tmpfile:
            format_fasta(testdata.PARSE_FASTA, tmpfile.name)

            with open(tmpfile.name) as observed:
                for line in observed:
                    self.assertLessEqual(len(line), expected_max_line_length)
