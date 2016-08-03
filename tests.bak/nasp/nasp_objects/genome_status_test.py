# coding: utf-8

__author__ = 'jtravis'

import unittest
from nasp.nasp_objects import InvalidContigName, GenomeStatus
from io import StringIO


class GenomeStatusTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = GenomeStatus()

    def test_add_contig(self):
        expected = ['foo']
        self.genome.add_contig(expected[0])
        self.assertEqual(expected, self.genome.get_contigs())

    def test_add_contig_with_contig_name_none_raises_invalidcontigname(self):
        with self.assertRaises(InvalidContigName):
            self.genome.add_contig(None)

    def test_add_contig_with_contig_name_empty_raises_invalidcontigname(self):
        with self.assertRaises(InvalidContigName):
            self.genome.add_contig("")

    def test_add_contig_with_duplicate_contig_name_raises_invalidcontigname(self):
        with self.assertRaises(InvalidContigName):
            self.genome.add_contig('duplicate')
            self.genome.add_contig('duplicate')

    def test_set_current_contig_with_contig_name_none(self):
        self.assertIsNone(self.genome.set_current_contig(None))

    # Special case: getter behavior when argument is None
    def test_set_current_contig_with_contig_name_none_returns_current_config(self):
        expected = 'foo'
        self.genome.add_contig(expected)
        self.assertEqual(expected, self.genome.set_current_contig(None))

    def test_set_current_contig_with_contig_name_invalid_and_create_contig_false_raises_invalidcontigname(self):
        with self.assertRaises(InvalidContigName):
            self.genome.set_current_contig('invalid_contig_name', False)

    def test_set_current_contig_with_contig_name_valid_and_create_contig_false(self):
        expected = 'foo'
        self.genome.add_contig(expected)
        self.genome.add_contig('bar')
        self.assertEqual(expected, self.genome.set_current_contig(expected, False))

    def test_get_contigs_empty(self):
        self.assertListEqual([], self.genome.get_contigs())

    def test_get_contigs(self):
        expected = ['a', 'b', 'c']
        # Reverse the list to verify the returned list is sorted
        for contig in reversed(expected):
            self.genome.add_contig(contig)
        self.assertListEqual(expected, self.genome.get_contigs())

    def test_append_contig_with_single_nucleotide(self):
        expected = ['A']
        contig_name = 'foo'
        self.genome.append_contig(expected[0], contig_name)
        self.assertListEqual(expected, self.genome.get_value(1, -1, contig_name))

    def test_append_contig_with_list(self):
        expected = ['A', 'C', 'G', 'T', 'U']
        contig_name = 'foo'
        self.genome.append_contig(expected, contig_name)
        self.assertListEqual(expected, self.genome.get_value(1, -1, contig_name))

    # TODO: make placeholder an optional final parameter
    def test_extend_contig(self):
        data = ['A', 'C', 'G', 'T', 'U']
        placeholder = '@'
        contig_name = 'foo'
        expected_length = len(data) + 5
        expected_content = list(data)
        expected_content.extend([placeholder] * 5)
        self.genome.append_contig(data, contig_name)
        self.genome.extend_contig(expected_length, placeholder, contig_name)
        # It allocated the additional space
        self.assertEqual(expected_length, self.genome.get_contig_length(contig_name))
        # It filled the allocated space after the data with the placeholder
        self.assertListEqual(expected_content, self.genome.get_value(1, -1, contig_name))

    def test_extend_contig_when_empty(self):
        expected = 42
        placeholder = '@'
        contig_name = 'foo'
        self.genome.add_contig(contig_name)
        self.genome.extend_contig(expected, placeholder, contig_name)
        # It allocated the additional space
        self.assertEqual(expected, self.genome.get_contig_length(contig_name))
        # It filled the undefined areas with the placeholder
        self.assertListEqual([placeholder] * expected, self.genome.get_value(1, -1, contig_name))

    def test_extend_contig_does_not_shrink_contig(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        expected = len(data)
        placeholder = '@'
        self.genome.append_contig(data, contig_name)
        self.genome.extend_contig(0, placeholder, contig_name)
        # It did not shrink the contig
        self.assertEqual(expected, self.genome.get_contig_length(contig_name))
        # It did not overwrite the existing data with a placeholder
        self.assertListEqual(data, self.genome.get_value(1, -1, contig_name))

    def test_set_value(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.set_value(data, 1, contig_name=contig_name)
        # Given a list, it should set the range of elements beginning with the start position
        self.assertListEqual(data, self.genome.get_value(1, -1, contig_name))
        self.genome.set_value('A', 3, contig_name=contig_name)
        data[2] = 'A'
        # Given a character, it should set the element at the start position
        self.assertListEqual(data, self.genome.get_value(1, -1, contig_name))

    # TODO: Should it instead raise an exception?
    def test_set_value_out_of_bounds(self):
        contig_name = 'foo'
        placeholder = '!'
        expected = [placeholder] * 4
        expected.append('A')
        # It should fill the interim elements with the placeholder when setting an out of bounds element
        self.genome.set_value('A', 5, missing_range_filler="!", contig_name=contig_name)
        self.assertEqual(expected, self.genome.get_value(1, -1, contig_name))

    # FIXME: get_value returns 3 data types - list, string, None
    def test_get_value_with_first_position(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.append_contig(data, contig_name)
        self.assertEqual('A', self.genome.get_value(1, contig_name=contig_name))

    def test_get_value_with_last_position(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.append_contig(data, contig_name)
        self.assertEqual(data, self.genome.get_value(1, -1, contig_name))

    def test_get_value_empty(self):
        contig_name = 'foo'
        self.genome.add_contig(contig_name)
        self.assertEqual([], self.genome.get_value(1, contig_name=contig_name))

    def test_get_value_with_contig_name_invalid(self):
        self.assertEqual([], self.genome.get_value(1, contig_name='invalid'))

    def test_get_value_with_invalid_position_raises_indexerror(self):
        with self.assertRaises(IndexError):
            contig_name = 'foo'
            data = ['A', 'C', 'G', 'T', 'U']
            self.genome.append_contig(data, contig_name)
            self.genome.get_value(0, contig_name=contig_name)

    def test_get_value_range_of_one(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.append_contig(data, contig_name)
        self.assertListEqual(['C'], self.genome.get_value(2, 2, contig_name))

    def test_get_value_reversed_position_values(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.append_contig(data, contig_name)
        self.assertEqual([], self.genome.get_value(2, 1, contig_name=contig_name))

    # TODO: Could this hide index errors?
    def test_get_value_outside_range(self):
        contig_name = 'foo'
        expected = ['A', 'C', 'G', 'T', 'U']
        self.genome.append_contig(expected, contig_name)
        self.assertListEqual(expected, self.genome.get_value(1, 10, contig_name=contig_name))

    # TODO: Is there a conflict if there are different placeholders?
    # Example: contig_extend fills allocated region with placeholder A,
    # but get_value retrieves an even larger region filling it with placeholder B
    # the result might be: [d,a,t,a,A,A,A,A,B,B,B,B]
    def test_get_value_outside_range_with_filler(self):
        contig_name = 'foo'
        placeholder = '@'
        data = ['A', 'C', 'G', 'T', 'U']
        expected = list(data)
        expected.extend([placeholder] * 5)
        self.genome.append_contig(data, contig_name)
        # It should fill out of bounds region with the placeholder
        self.assertListEqual(expected, self.genome.get_value(1, 10, contig_name=contig_name, placeholder=placeholder))

    def test_get_contig_length(self):
        contig_name = 'foo'
        data = ['A', 'C', 'G', 'T', 'U']
        self.genome.add_contig(contig_name)
        self.genome.append_contig(data, contig_name)
        self.assertEqual(len(data), self.genome.get_contig_length(contig_name))

    def test_get_contig_length_when_empty(self):
        contig_name = 'foo'
        self.genome.add_contig(contig_name)
        self.assertEqual(0, self.genome.get_contig_length(contig_name))

    def test_get_contig_length_with_contig_name_invalid(self):
        self.assertEqual(0, self.genome.get_contig_length('invalid'))

    def test_get_contig_length_with_contig_name_none(self):
        # FIXME: test fails without self.genome.add_contig(None)
        self.assertEqual(0, self.genome.get_contig_length(None))

    def test_send_to_fasta_handle(self):
        # contigs are sorted by name, not the order they were added
        expected = ">bar\nUTGCA\n>baz\nAUCTG\n>foo\nACGTU\n"
        data = [
            ('foo', ['A', 'C', 'G', 'T', 'U']),
            ('bar', ['U', 'T', 'G', 'C', 'A']),
            ('baz', ['A', 'U', 'C', 'T', 'G'])
        ]
        for contig in data:
            self.genome.append_contig(contig[1], contig[0])
        # Raises exception if contigs include None
        # TypeError: unorderable types: NoneType() < str()
        #self.genome.add_contig(None)
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle)
            self.assertEqual(expected, handle.getvalue())

    def test_send_to_fasta_handle_with_prefix(self):
        # contigs are sorted by name, not the order they were added
        expected = '>prefixbar\nUTGCA\n>prefixbaz\nAUCTG\n>prefixfoo\nACGTU\n'
        data = [
            ('foo', ['A', 'C', 'G', 'T', 'U']),
            ('bar', ['U', 'T', 'G', 'C', 'A']),
            ('baz', ['A', 'U', 'C', 'T', 'G'])
        ]
        for contig in data:
            self.genome.append_contig(contig[1], contig[0])
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle, contig_prefix='prefix')
            self.assertEqual(expected, handle.getvalue())

    def test_send_to_fasta_handle_when_empty(self):
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle)
            self.assertEqual("", handle.getvalue())

    def test_send_to_fasta_handle_with_max_chars_per_line(self):
        expected = '>foobarbaz\nACG\nTUA\nCGT\nUAC\nGTU\n'
        self.genome.append_contig(['A', 'C', 'G', 'T', 'U'] * 3, 'foobarbaz')
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle, max_chars_per_line=3)
            self.assertEqual(expected, handle.getvalue())

    def test_send_to_fasta_handle_without_max_chars_per_line(self):
        expected = '>foobarbaz\nACGTUACGTUACGTU\n'
        self.genome.append_contig(['A', 'C', 'G', 'T', 'U'] * 3, 'foobarbaz')
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle, max_chars_per_line=0)
            self.assertEqual(expected, handle.getvalue())

    def test_send_to_fasta_handle_appends_to_existing_files(self):
        expected = '>foo\nACGTU\n'
        self.genome.append_contig(['A', 'C', 'G', 'T', 'U'], 'foo')
        with StringIO() as handle:
            self.genome.send_to_fasta_handle(handle)
            self.assertEqual(expected, handle.getvalue())
            self.genome.send_to_fasta_handle(handle)
            self.assertEqual(expected + expected, handle.getvalue())

    def test_write_to_fasta_file(self):
        pass

if __name__ == '__main__':
    unittest.main()
