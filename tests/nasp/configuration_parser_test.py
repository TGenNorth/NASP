#!/usr/bin/env python3
# coding=utf-8

__author__ = "jtravis"

from nasp import configuration_parser
import os
import unittest
import tempfile
from xml.etree import ElementTree

# If parse_config accepted a file object it could be passed a tempfile.SpooledTemporaryFile or StringIO
# instead of a temporary file
#
# The NaspInputSchema.xsd does not match the NaspInputExample*.xml
# Is it valid for ReadFolder to contain both Read and ReadPair elements?
#
# Remove global variables
#
# configuration_parser is used by dispatcher and nasp to write and parse config
#
# configuration: written to by _parse_{applications, files, options} and returned by parse_config
# read_list is written to by _{get, find}_reads and read by _parse_files. _find_reads is called by _get_reads only
# fasta_list ditto
# bam_list written to by _{find_files, get_bams} and read by _parse_files
# vcf_list written to by _{find_files, get_vcfs} and read by _parse_files
# aligner_list written and read by _parse_applications
# snpcaller_list ditto
#
# _getReads returns num_reads which is never read and writes global read_list that is read


class ConfigurationParserTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # @unittest.skip("Not Implemented")
    # def test_parse_args(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_parse_options(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_find_reads(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_find_files(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_get_reads(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_get_fastas(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_get_bams(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_get_vcfs(self):
    #     raise NotImplementedError

    def test_parse_files_empty_files_node(self):
        expected = {
            'reads': [],
            'assemblies': [],
            'alignments': [],
            'vcfs': []
        }
        files_node = ElementTree.fromstring('<Files></Files>')
        configuration_parser._parse_files(files_node)
        self.assertTrue(all(configuration_parser.configuration[key] == expected[key] for key in expected.keys()),
                        "The configuration must contain all the expected key value pairs")

    # TODO: should have the same outcome as test_parse_files_missing_read_folder_path_attribute
    def test_parse_files_empty_read_folder_node(self):
        expected = {
            'reads': [],
            'assemblies': [],
            'alignments': [],
            'vcfs': []
        }
        files_node = ElementTree.fromstring('<Files>'
                                            '   <ReadFolder></ReadFolder>'
                                            '</Files>')
        configuration_parser._parse_files(files_node)
        self.assertTrue(all(configuration_parser.configuration[key] == expected[key] for key in expected.keys()),
                        "The configuration must contain all the expected key value pairs")

    # AttributeError: 'NoneType' object has no attribute 'startswith'
    def test_parse_files_missing_read_folder_path_attribute(self):
        expected = {
            'reads': [],
            'assemblies': [],
            'alignments': [],
            'vcfs': []
        }
        files_node = ElementTree.fromstring('<Files>'
                                            '   <ReadFolder>'
                                            '       <Read></Read>'
                                            '   </ReadFolder>'
                                            '</Files>')
        self.assertRaises(AttributeError, configuration_parser._parse_files, files_node)
        # configuration_parser._parse_files(files_node)
        # self.assertTrue(all(configuration_parser.configuration[key] == expected[key] for key in expected.keys()),
        #                "The configuration must contain all the expected key value pairs")

    # @unittest.skip("Not Implemented")
    # def test_get_application(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_parse_applications(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_write_reads(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_write_files(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_write_application(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_write_config_node(self):
    #     raise NotImplementedError
    #
    # @unittest.skip("Not Implemented")
    # def test_write_config(self):
    #     raise NotImplementedError

    # FileNotFound
    def test_parse_config_when_file_does_not_exist(self):
        self.assertRaises(FileNotFoundError, configuration_parser.parse_config, '')

    # xml.etree.ElementTree.ParseError: no element found
    def test_parse_config_empty_config_raises_parse_error(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            self.assertRaises(ElementTree.ParseError, configuration_parser.parse_config, config.name)

    # xml.etree.ElementTree.ParseError: mismatched tag
    def test_parse_config_malformed_xml_raises_parse_error(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('')
            config.seek(0)
            self.assertRaises(ElementTree.ParseError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'findtext'
    def test_parse_config_root_element_only_raises_attribute_error(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData></NaspInputData>')
            config.seek(0)
            self.assertRaises(AttributeError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'get'
    def test_parse_config_missing_reference_node(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData>'
                         '  <Options></Options>'
                         '  <Files></Files>'
                         '  <ExternalApplications></ExternalApplications>'
                         '</NaspInputData>')
            config.seek(0)
            self.assertRaises(AttributeError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'find'
    def test_parse_config_missing_filters_node(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData>'
                         '  <Options>'
                         '    <Reference></Reference>'
                         '  </Options>'
                         '  <Files></Files>'
                         '  <ExternalApplications></ExternalApplications>'
                         '</NaspInputData>')
            config.seek(0)
            self.assertRaises(AttributeError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'get'
    def test_parse_config_missing_samtools_node(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData>'
                         '  <Options>'
                         '    <Reference></Reference>'
                         '    <Filters></Filters>'
                         '  </Options>'
                         '  <Files></Files>'
                         '  <ExternalApplications></ExternalApplications>'
                         '</NaspInputData>')
            config.seek(0)
            self.assertRaises(AttributeError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'get'
    def test_parse_config_missing_samtools_node(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData>'
                         '  <Options>'
                         '    <Reference></Reference>'
                         '    <Filters></Filters>'
                         '  </Options>'
                         '  <Files></Files>'
                         '  <ExternalApplications>'
                         '      <Samtools></Samtools>'
                         '  </ExternalApplications>'
                         '</NaspInputData>')
            config.seek(0)
            self.assertRaises(AttributeError, configuration_parser.parse_config, config.name)

    # AttributeError: 'NoneType' object has no attribute 'find'
    def test_parse_config_with_minimum_elements(self):
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write('<?xml version="1.0" encoding="UTF-8"?>'
                         '<NaspInputData>'
                         '  <Options>'
                         '    <Reference></Reference>'
                         '    <Filters></Filters>'
                         '  </Options>'
                         '  <Files></Files>'
                         '  <ExternalApplications>'
                         '      <Samtools></Samtools>'
                         '      <MatrixGenerator></MatrixGenerator>'
                         '  </ExternalApplications>'
                         '</NaspInputData>')
            config.seek(0)
            expected = {
                'run_name': None,
                'output_folder': None,
                'reads': [],
                'vcfs': [],
                'matrix_generator': ('MatrixGenerator', None, '', {}),
                'aligners': [],
                'reference': (None, None),
                'snpcallers': [],
                'job_submitter': None,
                'alignments': [],
                'samtools': ('Samtools', None, '', {}),
                'assemblies': [],
                'find_dups': None
            }
            self.assertDictEqual(expected, configuration_parser.parse_config(config.name))


from tempfile import TemporaryDirectory


class ConfigurationParserWriteReadConfigTestCase(unittest.TestCase):
    maxDiff = None

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.expected = {
            'run_name': 'test_run_name',
            # 'output_folder': tmpdir,
            'reference': ('test_reference', 'test_reference.fasta'),
            # TODO: It should initialize to False if unspecified
            'find_dups': 'False',
            'job_submitter': 'None',
            'reads': [],
            'assemblies': [],
            'alignments': [],
            'vcfs': [],
            'index': ('Index', 'path/to/bin', '-args', {
                'name': 'index_job',
                'mem_requested': '1',
                'num_cpus': '2',
                'walltime': '3',
                'queue': 'four',
                'args': '-five'
            }),
            'matrix_generator': ('MatrixGenerator', 'path/to/bin/vtm', '-args', {
                'name': 'matrix_job',
                'mem_requested': '1',
                'num_cpus': '2',
                'walltime': '3',
                'queue': 'four',
                'args': '-five'
            }),
            'samtools': ('Samtools', 'path/to/bin/vtm', '-args', {}),
            'aligners': [
                ('BWA-mem', 'path/to/bwa', '-args', {
                    'name': 'bwa_align_job',
                    'mem_requested': '1',
                    'num_cpus': '2',
                    'walltime': '3',
                    'queue': 'four',
                    'args': '-five'
                }),
                ('Bowtie2', 'path/to/bt2', '-args', {
                    'name': 'bowtie_align_job',
                    'mem_requested': '1',
                    'num_cpus': '2',
                    'walltime': '3',
                    'queue': 'four',
                    'args': '-five'
                }),
            ],
            'snpcallers': [
                ('GATK', 'path/to/gatk', '-args', {
                    'name': 'gatk_snpcall_job',
                    'mem_requested': '1',
                    'num_cpus': '2',
                    'walltime': '3',
                    'queue': 'four',
                    'args': '-five'
                }),
                ('VarScan', 'path/to/varscan', '-args', {
                    'name': 'varscan_snpcall_job',
                    'mem_requested': '1',
                    'num_cpus': '2',
                    'walltime': '3',
                    'queue': 'four',
                    'args': '-five'
                }),
            ]
        }

    def tearDown(self):
        # FIXME: The configuration_parser global variables prevent the unit tests from running in isolation.
        configuration_parser.configuration = {}
        configuration_parser.read_list = []
        configuration_parser.fasta_list = []
        configuration_parser.bam_list = []
        configuration_parser.vcf_list = []
        configuration_parser.aligner_list = []
        configuration_parser.snpcaller_list = []


def test_it_reads_the_same_data_it_writes(self):
    with TemporaryDirectory(suffix='fjdksnfsdk') as tmpdir:
        self.expected['output_folder'] = tmpdir

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)


@unittest.skip("Not implemented")
def test_it_should_warn_when_the_read_element_does_not_contain_read_files(self):
    with TemporaryDirectory() as tmpdir:
        self.expected['output_folder'] = tmpdir

        self.expected['reads'] = [
            ('not_a_valid_read_file', '', '')
        ]

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)

def test_it_stores_reads(self):
    with TemporaryDirectory() as tmpdir:
        self.expected['output_folder'] = tmpdir

        self.expected['reads'] = [
            ('foo', 'bar', 'baz')
        ]

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)


def test_it_stores_assemblies(self):
    with TemporaryDirectory() as tmpdir:
        self.expected['output_folder'] = tmpdir

        self.expected['assemblies'] = [
            ('foo', 'bar')
        ]

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)


def test_it_stores_alignments(self):
    with TemporaryDirectory() as tmpdir:
        self.expected['output_folder'] = tmpdir

        self.expected['alignments'] = [
            ('foo', 'bar')
        ]

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)


def test_it_stores_vcfs(self):
    with TemporaryDirectory() as tmpdir:
        self.expected['output_folder'] = tmpdir

        self.expected['vcfs'] = [
            ('foo', 'bar')
        ]

        configuration_parser.write_config(self.expected)

        observed = configuration_parser.parse_config(os.path.join(tmpdir, 'test_run_name-config.xml'))

        self.assertDictEqual(self.expected, observed)


if __name__ == "__main__":
    unittest.main()
