"""
The following tests verify the I/O functions create the expected files with the expected content and format.
They are not expected to validate the data passed to them.
"""

__author__ = 'jtravis'

import unittest
import os
from tempfile import TemporaryDirectory
from copy import deepcopy
from collections import Counter

from nasp.vtm.analyze import PositionInfo
from nasp.vtm import write_matrix
from nasp.vtm.parse import FastaContig


contig_name = 'TestContig'

identifiers = (
    'sample1::aligner,snpcaller',
    'sample2::aligner1,snpcaller',
    'sample2::aligner2,snpcaller',
    'sample3::aligner,snpcaller'
)

# TODO: It should handle \n and \r\n line endings
# TODO: vcf column assertions

# A non-missingdata position
position1 = PositionInfo(
    # General Stats
    is_all_called=True,
    is_reference_clean=True,
    is_reference_duplicated=True,
    is_all_passed_coverage=True,
    is_all_passed_proportion=True,
    is_all_passed_consensus=True,
    is_all_quality_breadth=True,
    is_best_snp=True,

    # NOTE: Would it increase performance if this were a tuple of tuples and we avoided using append?
    # Sample Stats - A list of list of Counters representing the stats for each analysis file grouped by sample.
    all_sample_stats=(
        Counter({
            'was_called': 0,
            'passed_coverage_filter': 0,
            'passed_proportion_filter': 0,
            'quality_breadth': 0,
            'called_reference': 0,
            'called_snp': 0,
            'called_degen': 0
        }),

        # all - True if true in all of the analysis_stats for the same sample.
        Counter({
            'was_called': 1,
            'passed_coverage_filter': 1,
            'passed_proportion_filter': 1,
            'quality_breadth': 1,
            'called_reference': 1,
            'called_snp': 1,
            'called_degen': 1
        })),

    # Missing Data Matrix condition - at least one SampleAnalysis passes quality_breadth and is a SNP.
    is_missing_matrix=True,

    # NASP Master Matrix
    # Counters
    called_reference=1,
    called_snp=2,
    passed_coverage_filter=3,
    passed_proportion_filter=4,
    num_A=5,
    num_C=6,
    num_G=7,
    num_T=8,
    num_N=9,
    # Strings
    call_str='ACGTRYKMSWBDHVN.acgtrykmswbdhvn',
    masked_call_str='ACGTRYKMSWBDHVN.acgtrykmswbdhvn',
    CallWasMade='YYYYNNNNNNNNNNNNYYYYNNNNNNNNNNN',
    PassedDepthFilter='',
    PassedProportionFilter='',
    Pattern=''
)

# A missingdata position
position2 = deepcopy(position1)

#
position3 = deepcopy(position1)


class VcfMetadataTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_vcf_metadata(self):
        identifiers = ('Sample{0}'.format(x) for x in range(5))
        contigs = tuple(FastaContig('contig{0}'.format(x), x, 0, '', True) for x in range(5))
        expected = '\n'.join((
            '##fileFormat=VCFv4.2',
            '##source=NASPvtest_version',
            '##contig=<ID="contig0",length=0>',
            '##contig=<ID="contig1",length=1>',
            '##contig=<ID="contig2",length=2>',
            '##contig=<ID="contig3",length=3>',
            '##contig=<ID="contig4",length=4>##SAMPLE=<ID="Sample0",Genomes="Sample0",Mixture=1.0>',
            '##SAMPLE=<ID="Sample1",Genomes="Sample1",Mixture=1.0>',
            '##SAMPLE=<ID="Sample2",Genomes="Sample2",Mixture=1.0>',
            '##SAMPLE=<ID="Sample3",Genomes="Sample3",Mixture=1.0>',
            '##SAMPLE=<ID="Sample4",Genomes="Sample4",Mixture=1.0>##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
            '##FILTER=<ID=NoCall,Description="No call for this sample at this position">',
            '##FILTER=<ID=CovFail,Description="Insufficient depth of coverage for this sample at this position">',
            '##FILTER=<ID=PropFail,Description="Insufficient proportion of reads were variant for this sample at this position">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=FT,Number=1,Type=String,Description="Filters that failed for this sample at this position">\n',
        ))

        observed = write_matrix.get_vcf_metadata('test_version', identifiers, contigs)

        self.assertEqual(expected, observed)


class GetHeaderTestCase(unittest.TestCase):
    def test_vcf(self):
        expected = (
            '#CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO',
            'FORMAT',
            'sample1::aligner,snpcaller',
            'sample2::aligner1,snpcaller',
            'sample2::aligner2,snpcaller',
            'sample3::aligner,snpcaller'
        )
        self.assertEqual(expected, write_matrix.get_header('vcf', identifiers))

    def test_master(self):
        expected = (
            'LocusID',
            'Reference',
            'sample1::aligner,snpcaller',
            'sample2::aligner1,snpcaller',
            'sample2::aligner2,snpcaller',
            'sample3::aligner,snpcaller',
            '#SNPcall',
            '#Indelcall',
            '#Refcall',
            '#CallWasMade',
            '#PassedDepthFilter',
            '#PassedProportionFilter',
            '#A',
            '#C',
            '#G',
            '#T',
            '#Indel',
            '#NXdegen',
            'Contig',
            'Position',
            'InDupRegion',
            'SampleConsensus',
            'CallWasMade',
            'PassedDepthFilter',
            'PassedProportionFilter',
            'Pattern',
            'Pattern#'
        )
        self.assertEqual(expected, write_matrix.get_header('master', identifiers))

    def test_bestsnp(self):
        expected = (
            'LocusID',
            'Reference',
            'sample1::aligner,snpcaller',
            'sample2::aligner1,snpcaller',
            'sample2::aligner2,snpcaller',
            'sample3::aligner,snpcaller',
            '#SNPcall',
            '#Indelcall',
            '#Refcall',
            '#CallWasMade',
            '#PassedDepthFilter',
            '#PassedProportionFilter',
            '#A',
            '#C',
            '#G',
            '#T',
            '#Indel',
            '#NXdegen',
            'Contig',
            'Position',
            'InDupRegion',
            'SampleConsensus',
            'Pattern',
            'Pattern#'
        )
        self.assertEqual(expected, write_matrix.get_header('bestsnp', identifiers))

    def test_missingdata(self):
        expected = (
            'LocusID',
            'Reference',
            'sample1::aligner,snpcaller',
            'sample2::aligner1,snpcaller',
            'sample2::aligner2,snpcaller',
            'sample3::aligner,snpcaller',
            '#SNPcall',
            '#Indelcall',
            '#Refcall',
            '#CallWasMade',
            '#PassedDepthFilter',
            '#PassedProportionFilter',
            '#A',
            '#C',
            '#G',
            '#T',
            '#Indel',
            '#NXdegen',
            'Contig',
            'Position',
            'InDupRegion',
            'SampleConsensus',
            'CallWasMade',
            'PassedDepthFilter',
            'PassedProportionFilter',
            'Pattern',
            'Pattern#'
        )
        self.assertEqual(expected, write_matrix.get_header('missingdata', identifiers))

    def test_withallrefpos(self):
        expected = (
            'LocusID',
            'Reference',
            'sample1::aligner,snpcaller',
            'sample2::aligner1,snpcaller',
            'sample2::aligner2,snpcaller',
            'sample3::aligner,snpcaller',
            '#SNPcall',
            '#Indelcall',
            '#Refcall',
            '#CallWasMade',
            '#PassedDepthFilter',
            '#PassedProportionFilter',
            '#A',
            '#C',
            '#G',
            '#T',
            '#Indel',
            '#NXdegen',
            'Contig',
            'Position',
            'InDupRegion',
            'SampleConsensus',
            'Pattern',
            'Pattern#'
        )
        self.assertEqual(expected, write_matrix.get_header('withallrefpos', identifiers))

    def test_undefined_header(self):
        # It should warn the developer if they request an undefined header.
        with self.assertRaises(ValueError):
            write_matrix.get_header('undefined', identifiers)


class WriteMissingMatrixTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata = 'fake metadata\n'
        cls.positions = (position1, position2, position3)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_it_includes_missingdata_positions(self):
        """
        See Also:
            nasp.vtm.analyze.analyze_position
        """
        pass

    def test_it_excludes_non_missingdata_positions(self):
        """
        See Also:
            nasp.vtm.analyze.analyze_position
        """
        pass

    def test_write_missingdata_vcf(self):
        expected_files = ['TestContig_missingdata.vcf']
        expected_lines = (
            self.metadata,
            '\t'.join(write_matrix.get_header('vcf', identifiers)) + '\n',
            # TODO: test thresholds aren't hardcoded
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_missingdata_vcf(tmpdir, contig_name, identifiers, self.metadata)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

    def test_write_missingdata_matrix(self):
        expected_files = ['TestContig_missingdata.tsv']
        expected_lines = (
            '\t'.join(write_matrix.get_header('missing_data', identifiers)) + '\n',
            '',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_missingdata_matrix(tmpdir, contig_name, identifiers)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

    def test_write_master_matrix(self):
        expected_files = ['TestContig_master.tsv']
        expected_lines = (
            '\t'.join(write_matrix.get_header('master', identifiers)) + '\n',
            '',
        )

        print(expected_lines)

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_master_matrix(tmpdir, contig_name, identifiers)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

    def test_write_bestsnp_vcf(self):
        expected_files = ['TestContig_bestsnp.vcf']
        expected_lines = (
            self.metadata,
            '\t'.join(write_matrix.get_header('vcf', identifiers)) + '\n',
            'TestContig	1	.	A	g,d,B,D,m,Y,M,T,r,R,k,v,K,W,y,H,V,N,a,S,t,.,G,n,s,C,b,h,w,c	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	k,D,W,r,g,m,T,K,w,c,a,y,H,h,v,s,N,M,C,n,t,b,.,V,G,B,R,d,Y,S	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	v,V,H,g,.,W,c,m,k,G,K,d,a,M,y,r,Y,t,B,T,D,S,b,n,s,h,N,C,w,R	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	S,V,w,b,M,C,K,k,H,T,R,W,g,n,d,a,B,Y,c,G,m,N,s,h,y,D,.,t,v,r	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	K,H,h,w,g,b,T,B,v,n,.,W,D,t,V,y,S,d,s,a,k,G,N,m,r,Y,M,R,C,c	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	c,R,t,n,C,K,d,H,N,a,r,k,Y,m,g,.,b,B,G,h,v,W,T,y,s,V,D,M,w,S	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	Y,H,h,N,v,T,m,G,S,B,w,s,R,c,W,n,d,r,a,k,K,.,V,y,b,g,C,M,D,t	.	PASS	AN=31;NS=3	GT:FT',
            'TestContig	1	.	A	T,B,k,G,Y,D,H,w,m,c,a,.,t,d,g,M,n,K,b,s,C,N,W,S,v,V,R,r,y,h	.	PASS	AN=31;NS=3	GT:FT				',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_bestsnp_vcf(tmpdir, contig_name, identifiers, self.metadata)
            writer.send(None)

            for position in self.positions:
                print("POSITION", position)
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    print(line)
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


    def test_write_sample_stats(self):
        expected_files = ['TestContig_sample_stats.tsv']
        expected_lines = (
            self.metadata,
            '\t'.join(write_matrix.get_header('vcf', identifiers)) + '\n',
            '',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_sample_stats(tmpdir, contig_name, identifiers, self.metadata)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


    def test_write_general_stats(self):
        expected_files = ['TestContig_general_stats.tsv']
        expected_lines = (
            self.metadata,
            '\t'.join(write_matrix.get_header('vcf', identifiers)) + '\n',
            '',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_general_stats(tmpdir, sample_stats)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


    def test_write_bestsnp_matrix(self):
        expected_files = ['TestContig_bestsnp_matrix.tsv']
        expected_lines = (
            self.metadata,
            '\t'.join(write_matrix.get_header('bestsnp', identifiers)) + '\n',
            '',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_bestsnp_matrix(tmpdir, contig_name, sample_groups)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


    def test_write_withallrefpos_matrix(self):
        expected_files = ['TestContig_withallrefpos.tsv']
        expected_lines = (
            '\t'.join(write_matrix.get_header('withallrefpos', identifiers)) + '\n',
            '\t'.join((
            'TestContig::1', 'A', 'C', 'G', 'T', 'R', '2', '0', '1', '7/30', '3/30', '4/30', '5', '6', '7', '8', '0',
            '9', 'TestContig', '1', 'True', 'True', '', '')) + '\n',
            '\t'.join((
            'TestContig::2', 'A', 'C', 'G', 'T', 'R', '2', '0', '1', '7/30', '3/30', '4/30', '5', '6', '7', '8', '0',
            '9', 'TestContig', '2', 'True', 'True', '', '')) + '\n',
            '\t'.join((
            'TestContig::3', 'A', 'C', 'G', 'T', 'R', '2', '0', '1', '7/30', '3/30', '4/30', '5', '6', '7', '8', '0',
            '9', 'TestContig', '3', 'True', 'True', '', '')) + '\n',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_withallrefpos_matrix(tmpdir, contig_name, identifiers)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


    def test_write_bestsnp_snpfasta(self):
        expected_files = [
            'TestContig_sample1::aligner,snpcaller_bestsnp.fasta',
            'TestContig_sample2::aligner1,snpcaller_bestsnp.fasta',
            'TestContig_sample2::aligner2,snpcaller_bestsnp.fasta',
            'TestContig_sample3::aligner,snpcaller_bestsnp.fasta'
        ]
        expected_lines = (
            'AAA',
        )

        with TemporaryDirectory() as tmpdir:
            writer = write_matrix.write_bestsnp_snpfasta(tmpdir, contig_name, identifiers)
            writer.send(None)

            for position in self.positions:
                writer.send(position)
            writer.close()

            # The file was created.
            self.assertListEqual(expected_files, os.listdir(tmpdir))

            with open(os.path.join(tmpdir, expected_files[0])) as handle:
                # The file contains all the expected rows.
                for expected_line, line in zip(expected_lines, handle):
                    self.assertEqual(expected_line, line)

                # The file does not contain any unexpected rows.
                self.assertEqual([], handle.readlines())

            # No other artifacts were created in the tmpdir.
            self.assertListEqual(expected_files, os.listdir(tmpdir))


class WriteMissingSnpfastaTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # cls.metadata = 'fake metadata\n'
        cls.positions = (position1, position2, position3)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_write_masked_base_calls_for_each_sample_analysis(self):
        expected_files = ['test_contig_a_missingdata.fasta', 'test_contig_b_missingdata.fasta',
                          'test_contig_c_missingdata.fasta']
        with TemporaryDirectory() as tmpdir:
            coroutine = write_matrix.write_missingdata_snpfasta(tmpdir, 'test_contig', ['a', 'b', 'c'])
            coroutine.send(None)

            # The files were created
            self.assertListEqual(os.listdir(tmpdir), expected_files)

            for position in self.positions:
                coroutine.send(position)

                # TODO: Assert each file has its corresponding calls


class WriteBestsnpSnpfastaTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # cls.metadata = 'fake metadata\n'
        cls.positions = (position1, position2, position3)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_write_masked_base_calls_for_each_sample_analysis(self):
        expected_files = ['test_contig_a_bestsnp.fasta', 'test_contig_b_bestsnp.fasta', 'test_contig_c_bestsnp.fasta']
        with TemporaryDirectory() as tmpdir:
            coroutine = write_matrix.write_bestsnp_snpfasta(tmpdir, 'test_contig', ['a', 'b', 'c'])
            coroutine.send(None)

            # The files were created
            self.assertListEqual(os.listdir(tmpdir), expected_files)

            for position in self.positions:
                coroutine.send(position)

                # TODO: Assert each file has its corresponding calls
