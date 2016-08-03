__author__ = 'jtravis'

import unittest
import tempfile
from nasp.vtm.matrix_DTO import NaspFile, parse_dto, write_dto

# TODO: vtm with empty dto?
# TODO: Should parameters return an empty dict instead of None?
# dto = {
# 'parameters': None,
# 'vcf': (),
# 'frankenfasta': ()
# }

# FIXME: dto parsing is crude. The following is an example of the output where empty strings are assigned the value
# None and the parameters, instead of being a dict of its elements it is a string with a newline and pretty print
# indentation.
# + {'filter_matrix_format': None,
# +  'frankenfasta': (NaspFile(path='file1.frankenfasta', name='name1', aligner='aligner1', snpcaller=None),
# +                   NaspFile(path='file2.frankenfasta', name='name2', aligner='aligner2', snpcaller=None),
# +                   NaspFile(path='file3.frankenfasta', name='name3', aligner='aligner3', snpcaller=None)),
# +  'matrix_folder': None,
# +  'minimum_coverage': None,
# +  'minimum_proportion': None,
# +  'parameters': '\n          ',
# +  'reference_dups': None,
# +  'reference_fasta': None,
# +  'stats_folder': None,
# +  'vcf': (NaspFile(path='file1.vcf', name='name1', aligner='aligner1', snpcaller='snpcaller1'),
# +          NaspFile(path='file2.vcf', name='name2', aligner='aligner2', snpcaller='snpcaller2'),
# +          NaspFile(path='file3.vcf', name='name3', aligner='aligner3', snpcaller='snpcaller3'))}


# TODO: Parse empty dto? extra elements? typo? missing element?
class MatrixDtoTestCase(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.dto = {
            'parameters': {
                'stats_folder': '',
                'filter_matrix_format': '',
                'reference_fasta': '',
                'minimum_proportion': '',
                'reference_dups': '',
                'matrix_folder': '',
                'minimum_coverage': ''
            },
            'vcf': (
                ('name1', 'aligner1', 'snpcaller1', 'file1.vcf'),
                ('name2', 'aligner2', 'snpcaller2', 'file2.vcf'),
                ('name3', 'aligner3', 'snpcaller3', 'file3.vcf')
            ),
            'frankenfasta': (
                ('name1', 'aligner1', 'file1.frankenfasta'),
                ('name2', 'aligner2', 'file2.frankenfasta'),
                ('name3', 'aligner3', 'file3.frankenfasta')
            )

        }

    def test_write_matrix_dto(self):
        expected = {
            'minimum_proportion': None,
            'matrix_folder': None,
            'reference_dups': None,
            'reference_fasta': None,
            'vcf': (
                NaspFile(path='file1.vcf', name='name1', aligner='aligner1', snpcaller='snpcaller1'),
                NaspFile(path='file2.vcf', name='name2', aligner='aligner2', snpcaller='snpcaller2'),
                NaspFile(path='file3.vcf', name='name3', aligner='aligner3', snpcaller='snpcaller3')
            ),
            'minimum_coverage': None,
            'filter_matrix_format': None,
            'parameters': '\n          ',
            'stats_folder': None,
            'frankenfasta': (
                NaspFile(path='file1.frankenfasta', name='name1', aligner='aligner1', snpcaller=None),
                NaspFile(path='file2.frankenfasta', name='name2', aligner='aligner2', snpcaller=None),
                NaspFile(path='file3.frankenfasta', name='name3', aligner='aligner3', snpcaller=None)
            )
        }

        with tempfile.NamedTemporaryFile() as tmp:
            write_dto(self.dto['parameters'], self.dto['frankenfasta'], self.dto['vcf'], tmp.name)
            observed = parse_dto(tmp.name)

            self.assertDictEqual(expected, observed)
