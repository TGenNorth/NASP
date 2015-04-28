__author__ = 'jtravis'

import unittest
import tempfile
import types
import os

from nasp.vtm.parse import Contig, EmptyContig, FastaContig, VcfContig, Position


# TODO: It should handle \n and \r\n line endings


class PositionTestCase(unittest.TestCase):

    def test_simple_call_normalizes_to_uppercase_and_masks_degeneracies_with_N(self):
        observed = ['g', 'a', 't', 'c', 'd', 'x', '.']
        expected = ['G', 'A', 'T', 'C', 'N', 'N', 'N']

        for observe, expect in zip(observed, expected):
            pos = Position(call=observe, coverage='-', proportion='-')
            self.assertEqual(expect, pos.simple_call)

    def test_call_normalizes_to_uppercase(self):
        observed = ['g', 'a', 't', 'c', 'd', 'x', '.']
        expected = ['G', 'A', 'T', 'C', 'D', 'X', '.']

        for observe, expect in zip(observed, expected):
            pos = Position(call=observe, coverage='-', proportion='-')
            self.assertEqual(expect, pos.call)

    def test_multi_base_call_raises_exception(self):
        # TODO: use appropriate exception
        with self.assertRaises(Exception):
            Position(call='gatc', coverage='-', proportion='-')


class EmptyContigTestCase(unittest.TestCase):
    """
    EmptyContig is a placeholder for contigs that exist in the reference,
    but not in the sample. It should yield infinite empty positions.
    """

    def test_yields_empty_fasta_positions(self):
        contig_name = 'test_fasta'
        is_fasta = True
        contig = EmptyContig(contig_name, is_fasta=is_fasta)
        self.assertIsInstance(contig.positions, types.GeneratorType)
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(contig.positions))
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(contig.positions))
        self.assertEqual(contig_name, contig.name)
        self.assertEqual("EmptyContig(name='{contig_name}', is_fasta={is_fasta})".format(
            contig_name=contig_name,
            is_fasta=is_fasta
        ), repr(contig))

    def test_yields_empty_vcf_positions(self):
        contig_name = 'test_vcf'
        is_fasta = False
        contig = EmptyContig(contig_name, is_fasta=is_fasta)
        self.assertIsInstance(contig.positions, types.GeneratorType)
        self.assertEqual(Contig.VCF_EMPTY_POSITION, next(contig.positions))
        self.assertEqual(Contig.VCF_EMPTY_POSITION, next(contig.positions))
        self.assertEqual(contig_name, contig.name)
        self.assertEqual("EmptyContig(name='{contig_name}', is_fasta={is_fasta})".format(
            contig_name=contig_name,
            is_fasta=is_fasta
        ), repr(contig))


class FastaContigTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # Source fasta values
        contigs = (
            # No gap
            (
                ">contig0\n"
                "GATC\n"
                "GGAA\n"
            ),
            # Gap after contig
            (
                ">contig1\n"
                "GATC\n"
                "GGAA\n"
                "\n"
            # ),
            # # No linebreak
            # (
            #     ">contig2\n"
            #     "GATCGGAA"
            )
            # TODO: >80 characters contig?
        )

        # Expected values
        self.contigs_expected = (
            {
                'name': 'contig0',
                'file_position': 0 + len('>contig0\n'),
                'positions': (
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-'),
                    Position(call='T', coverage='-', proportion='-'),
                    Position(call='C', coverage='-', proportion='-'),
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-')
                )
            },
            {
                'name': 'contig1',
                'file_position': len(contigs[0]) + len('>contig1\n'),
                'positions': (
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-'),
                    Position(call='T', coverage='-', proportion='-'),
                    Position(call='C', coverage='-', proportion='-'),
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='G', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-'),
                    Position(call='A', coverage='-', proportion='-')
                )
            }
        )

        fasta_content = "".join(contigs)

        # Create a mock fasta
        self.fasta_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        self.fasta_file.write(fasta_content)
        self.fasta_file.seek(0)

        # Instantiate test contigs
        self.contig0 = FastaContig(
            self.contigs_expected[0]['name'],
            len(self.contigs_expected[0]['positions']),
            self.contigs_expected[0]['file_position'],
            self.fasta_file.name,
            is_reference=False
        )
        self.contig1 = FastaContig(
            self.contigs_expected[1]['name'],
            len(self.contigs_expected[1]['positions']),
            self.contigs_expected[1]['file_position'],
            self.fasta_file.name,
            is_reference=False
        )

        self.ref_contig0 = FastaContig(
            self.contigs_expected[0]['name'],
            len(self.contigs_expected[0]['positions']),
            self.contigs_expected[0]['file_position'],
            self.fasta_file.name,
            is_reference=True
        )
        self.ref_contig1 = FastaContig(
            self.contigs_expected[1]['name'],
            len(self.contigs_expected[1]['positions']),
            self.contigs_expected[1]['file_position'],
            self.fasta_file.name,
            is_reference=True
        )

    def tearDown(self):
        os.remove(self.fasta_file.name)

    def test_repr(self):
        contig = self.contigs_expected[0]

        expected = "FastaContig(name='{0}', length={1}, file_position={2}, file_path='{3}', is_reference={4})".format(
            contig['name'],
            len(contig['positions']),
            contig['file_position'],
            self.fasta_file.name,
            False
        )

        self.assertEqual(expected, repr(self.contig0))

    def test_name(self):
        self.assertEqual(self.contigs_expected[0]['name'], self.contig0.name)

    def test_length(self):
        self.assertEqual(len(self.contigs_expected[0]['positions']), len(self.contig0))

    def test_reference_fastas_return_finite_positions_zero_offset(self):
        """
        The reference fasta drives the analysis. When the reference has data, all the samples
        should return a corresponding value for the position. It is only file that should stop returning values
        when there is no data.

        It should read positions starting from the file_position.
        """
        expected = self.contigs_expected[0]['positions']
        positions = self.ref_contig0.positions

        # It should yield all the positions in the contig.
        for expected, observed in zip(expected, positions):
            self.assertEqual(expected, observed)

        # It should stop yielding positions after the contig is exhausted.
        with self.assertRaises(StopIteration):
            next(positions)

    def test_reference_fastas_return_finite_positions_nonzero_offset(self):
        """
        The reference fasta drives the analysis. When the reference has data, all the samples
        should return a corresponding value for the position. It is only file that should stop returning values
        when there is no data.

        It should read positions starting from the file_position.
        """
        expected = self.contigs_expected[1]['positions']
        positions = self.ref_contig1.positions

        # It should yield all the positions in the contig.
        for expected, observed in zip(expected, positions):
            self.assertEqual(expected, observed)

        # It should stop yielding positions after the contig is exhausted.
        with self.assertRaises(StopIteration):
            next(positions)

    def test_sample_fastas_return_infinite_positions(self):
        """
        A sample fasta must yield a corresponding position for each reference position.
        If there is no data, it must yield a empty position as a placeholder.
        """
        expected = self.contigs_expected[1]['positions']
        observed = self.contig1.positions

        # It should yield all the positions in the contig.
        for expect, observe in zip(expected, observed):
            self.assertEqual(expect, observe)

        # It should yield empty positions after the contig is exhausted.
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(observed))
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(observed))


class VcfContigTestCase(unittest.TestCase):

    @classmethod
    def setUp(self):
        self.file_path = './test_data/gatk.vcf'

        # # Expected values
        # self.contigs_expected = (
        #     {
        #         'name': 'contig0',
        #         'file_position': 0 + len('>contig0\n'),
        #         'positions': (
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='T', coverage='-', proportion='-'),
        #             Position(call='C', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-')
        #         )
        #     },
        #     {
        #         'name': 'contig1',
        #         'file_position': len(contigs[0]) + len('>contig1\n'),
        #         'positions': (
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='T', coverage='-', proportion='-'),
        #             Position(call='C', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-')
        #         )
        #     },
        #     {
        #         'name': 'contig2',
        #         'file_position': len(contigs[0]) + len(contigs[1]) + len('>contig2\n'),
        #         'positions': (
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='T', coverage='-', proportion='-'),
        #             Position(call='C', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='G', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-'),
        #             Position(call='A', coverage='-', proportion='-')
        #         ),
        #     }
        # )
        #
        # # Instantiate test contigs
        # self.contig0 = VcfContig(
        #     self.contigs_expected[0]['name'],
        #     len(self.contigs_expected[0]['positions']),
        #     self.contigs_expected[0]['file_position'],
        #     self.fasta_file.name,
        #     is_reference=False
        # )
        # self.contig1 = VcfContig(
        #     self.contigs_expected[1]['name'],
        #     len(self.contigs_expected[1]['positions']),
        #     self.contigs_expected[1]['file_position'],
        #     self.fasta_file.name,
        #     is_reference=False
        # )
        #
        # self.ref_contig0 = VcfContig(
        #     self.contigs_expected[0]['name'],
        #     len(self.contigs_expected[0]['positions']),
        #     self.contigs_expected[0]['file_position'],
        #     self.fasta_file.name,
        #     is_reference=True
        # )
        # self.ref_contig1 = VcfContig(
        #     self.contigs_expected[1]['name'],
        #     len(self.contigs_expected[1]['positions']),
        #     self.contigs_expected[1]['file_position'],
        #     self.fasta_file.name,
        #     is_reference=True
        # )

    def test_foo(self):
        from nasp.vtm.parse import Vcf
        vcf = Vcf(self.file_path, '', '', '')
        contig = vcf.get_contig('500WT1_test')

        positions = contig.positions

        for n, p in zip(range(2), positions):
            if n == 315:
                a = 'a'
                # break
            # self.assertEquals(contig.VCF_EMPTY_POSITION, p)

        print(next(positions))

        print(next(positions))

        for _, p in zip(range(6), contig.positions):
            print(p)
            self.assertNotEqual(contig.VCF_EMPTY_POSITION, p)

    def test_sample_vcfs_return_infinite_positions(self):
        contig_name = 'sample_contig'
        length = len(self.contig1_positions)
        file_position = 0
        filepath = self.fasta_file.name
        is_reference = False

        contig1 = FastaContig(contig_name, length, file_position, filepath, is_reference)
        positions = fasta.positions

        self.assertEqual(length, len(contig1))

        for expected, observed in zip(self.contig1_positions, positions):
            self.assertEqual(expected, observed)

        # It should yield empty positions after the contig is exhausted.
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(positions))
        self.assertEqual(Contig.FASTA_EMPTY_POSITION, next(positions))


        # class VcfTestCase(unittest.TestCase):
        #
        #     @classmethod
        #     def setUpClass(cls):
        #         pass
        #
        #     def setUp(self):
        #         # Create a mock fasta
        #         self.vcf_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        #         self.vcf_file.write(header)
        #         self.vcf_file.seek(0)
        #
        #     def tearDown(self):
        #         os.remove(self.vcf_file)
        #
        #     # VCF Syntax tests
        #
        #     # TODO: empty file?
        #
        #     # TODO: invalid rows?
        #
        #     def test_vcf_missing_header_raises_MalformedInputFileException(self):
        #         with self.assertRaises(MalformedInputFileException) as tempfile.NamedTemporaryFile(mode="w+") as handle:
        #             vcf = Vcf(handle.name, 'name', 'aligner', 'snpcaller')
        #
        #     def test_vcf_repr(self):
        #         with tempfile.NamedTemporaryFile(mode="w+") as handle:
        #             expected = "Vcf('{0}', 'name', 'aligner', 'snpcaller')".format(handle.name)
        #
        #             handle.write(header)
        #             handle.seek(0)
        #             vcf = Vcf(handle.name, 'name', 'aligner', 'snpcaller')
        #
        #         self.assertEqual(expected, repr(vcf))
        #
        #     def test_vcf_identifier(self):
        #         with tempfile.NamedTemporaryFile(mode="w+") as handle:
        #             expected = "name::aligner,snpcaller"
        #
        #             handle.write(header)
        #             handle.seek(0)
        #             vcf = Vcf(handle.name, 'name', 'aligner', 'snpcaller')
        #
        #         self.assertEqual(expected, vcf.identifier)
        #
        #     def test_vcf_contig_random_access(self):
        #         with tempfile.NamedTemporaryFile(mode='w+') as handle:
        #             handle.write(header)
        #             handle.seek(0)
        #
        #             vcf.get_contig('')
        #
        #             vcf = Vcf(handle.name, 'TestVcf', 'pre-aligned', 'pre-called')
        #
        #     def test_vcf_missing_contig(self):
        #
        #         with tempfile.NamedTemporaryFile(mode='w+') as handle:
        #             handle.write(header)
        #             handle.seek(0)
        #
        #             vcf = Vcf(handle.name, 'TestVcf', 'pre-aligned', 'pre-called')
        #             contig = vcf.get_contig('DoesNotExist')
        #             self.assertIsInstance(contig, EmptyContig)
