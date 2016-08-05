"""
The purpose of these tests is to ensure commands are properly shell escaped.

Characters that cause problems include pipes: '|', redirects: '>', newlines: '\n', etc
"""
__author__ = 'jtravis'

import unittest
from unittest.mock import Mock, call
from collections import namedtuple
#import os
import itertools
#from tempfile import TemporaryDirectory

from nasp import dispatcher

# TODO: @unittest.skipUnless(termios, 'tests require system with termios')

App = namedtuple('App', ['name', 'path', 'args', 'job_params'])
Sample = namedtuple('Sample', ['name', 'read1', 'read2'])


class DispatcherShellEscapeCommandsTestCase(unittest.TestCase):
    maxDiff = None

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # Monkey-patch mock _submit_job function
        # Each call to submit_job will return an incrementing "job_id"
        self.mock_submit_job = Mock(side_effect=(str(n) for n in itertools.count(start=1, step=1)))
        self.mock_release_hold = Mock()
        dispatcher._submit_job = self.mock_submit_job
        dispatcher._release_hold = self.mock_release_hold

        self.samples = {
            'paired_basic': Sample('NA10831_ATCACG_L002', 'NA10831_ATCACG_L002_R1_001.fastq.gz', 'NA10831_ATCACG_L002_R2_001.fastq.gz'),
            'paired_pipe': Sample('NA|10831_ATCACG_L002', 'NA|10831_ATCACG_L002_R1_001.fastq.gz', 'NA|10831_ATCACG_L002_R2_001.fastq.gz'),
            'single_basic': Sample('NA10831_ATCACG_L002', 'NA10831_ATCACG_L002.fastq.gz', ''),
            'single_pipe': Sample('NA|10831_ATCACG_L002', 'NA|10831_ATCACG_L002.fastq.gz', ''),
        }
        self.bwamem = App('bwamem', '/path/to/bwa', '-x "-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0"', {'num_cpus': 2})
        self.bwa = App('bwa', '/path/to/bwa', '-x "-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0"', {'num_cpus': 2})
        self.bowtie2 = App('bowtie2', '/path/to/bowtie2', '--very-sensitive-local --un pipe|in|name.fastq.gz --al "space in name.fastq.gz"', {'num_cpus': 2})
        self.novoalign = App('novoalign', '/path/to/novoalign', '-K mismatch:stats.txt -i MP 99-99 99,99', {'num_cpus': 2})
        self.snap = App('snap', '/path/to/snap', '--TODO', {'num_cpus': 2})
        self.samtools = App('samtools', '/path/to/samtools', '', {})
        self.job_submitter = 'pbs'
        self.reference = 'reference.fasta'
        self.output_folder = 'output_folder'
        self.index_job_id = ('jobid', 'action')


    def tearDown(self):
        pass

    def test_samtools_view_sort_index_pipe_command(self):
        tests = {
            'paired_basic': "/path/to/samtools view -S -b -h - | /path/to/samtools sort - NA10831_ATCACG_L002; /path/to/samtools index NA10831_ATCACG_L002.bam",

            'paired_pipe': "/path/to/samtools view -S -b -h - | /path/to/samtools sort - 'NA|10831_ATCACG_L002'; /path/to/samtools index 'NA|10831_ATCACG_L002.bam'",
        }

        for sample_type, expect in tests.items():
            result = dispatcher._samtools_view_sort_index_pipe_command(self.samtools.path, self.samples[sample_type].name)
            self.assertEqual(expect, result)


    def test_bwamem_command(self):
        tests = {
            'paired_basic': "/path/to/bwa mem -R '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0' -t 2 reference.fasta NA10831_ATCACG_L002_R1_001.fastq.gz NA10831_ATCACG_L002_R2_001.fastq.gz",

            'paired_pipe': "/path/to/bwa mem -R '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0' -t 2 reference.fasta 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz'",

            'single_basic': "/path/to/bwa mem -R '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0' -t 2 reference.fasta NA10831_ATCACG_L002.fastq.gz ",

            'single_pipe': "/path/to/bwa mem -R '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0' -t 2 reference.fasta 'NA|10831_ATCACG_L002.fastq.gz' "
        }

        for sample_type, expect in tests.items():
            result = dispatcher._bwamem_command(self.bwa.path, self.bwa.args, 2, self.reference, *self.samples[sample_type])
            self.assertEqual(expect, result)


    def test_bwa_command(self):
        tests = {
            'paired_basic': "/path/to/bwa aln  reference.fasta NA10831_ATCACG_L002_R1_001.fastq.gz -t 2 -f output_folder/bwa/NA10831_ATCACG_L002-R1.sai -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa aln  reference.fasta NA10831_ATCACG_L002_R2_001.fastq.gz -t 2 -f output_folder/bwa/NA10831_ATCACG_L002-R1.sai -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa sampe -r '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' reference.fasta output_folder/bwa/NA10831_ATCACG_L002-R1.sai output_folder/bwa/NA10831_ATCACG_L002-R2.sai NA10831_ATCACG_L002_R1_001.fastq.gz NA10831_ATCACG_L002_R2_001.fastq.gz -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'",

            'paired_pipe': "/path/to/bwa aln  reference.fasta 'NA|10831_ATCACG_L002_R1_001.fastq.gz' -t 2 -f 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa aln  reference.fasta 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -t 2 -f 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa sampe -r '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' reference.fasta 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' 'output_folder/bwa/NA|10831_ATCACG_L002-R2.sai' 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'",

            'single_basic': "/path/to/bwa aln  reference.fasta NA10831_ATCACG_L002.fastq.gz -t 2 -f output_folder/bwa/NA10831_ATCACG_L002.sai -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa samse -r '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' reference.fasta output_folder/bwa/NA10831_ATCACG_L002.sai NA10831_ATCACG_L002.fastq.gz -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'",

            'single_pipe': "/path/to/bwa aln  reference.fasta 'NA|10831_ATCACG_L002.fastq.gz' -t 2 -f 'output_folder/bwa/NA|10831_ATCACG_L002.sai' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa samse -r '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' reference.fasta 'output_folder/bwa/NA|10831_ATCACG_L002.sai' 'NA|10831_ATCACG_L002.fastq.gz' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'"
        }

        for sample_type, expect in tests.items():
            result = dispatcher._bwa_command(self.bwa.path, self.bwa.args, 2, self.reference, self.output_folder, *self.samples[sample_type])
            self.assertEqual(expect, result)


    def test_bowtie2_command(self):
        tests = {
            'paired_basic': "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg SM:NA10831_ATCACG_L002 --rg-id NA10831_ATCACG_L002 -x reference -1 NA10831_ATCACG_L002_R1_001.fastq.gz -2 NA10831_ATCACG_L002_R2_001.fastq.gz",

            'paired_pipe': "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg 'SM:NA|10831_ATCACG_L002' --rg-id 'NA|10831_ATCACG_L002' -x reference -1 'NA|10831_ATCACG_L002_R1_001.fastq.gz' -2 'NA|10831_ATCACG_L002_R2_001.fastq.gz'",

            'single_basic': "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg SM:NA10831_ATCACG_L002 --rg-id NA10831_ATCACG_L002 -x reference -U NA10831_ATCACG_L002.fastq.gz",

            'single_pipe': "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg 'SM:NA|10831_ATCACG_L002' --rg-id 'NA|10831_ATCACG_L002' -x reference -U 'NA|10831_ATCACG_L002.fastq.gz'"
        }

        for sample_type, expect in tests.items():
            result = dispatcher._bowtie2_command(self.bowtie2.path, self.bowtie2.args, 2, self.reference, *self.samples[sample_type])
            self.assertEqual(expect, result)


    def test_novoalign_command(self):
        tests = {
            'paired_basic': "/path/to/novoalign -d reference.fasta.idx -f NA10831_ATCACG_L002_R1_001.fastq.gz NA10831_ATCACG_L002_R2_001.fastq.gz -i PE 500,100 -c 2 -o SAM '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' -K mismatch:stats.txt -i MP 99-99 99,99",

            'paired_pipe': "/path/to/novoalign -d reference.fasta.idx -f 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -i PE 500,100 -c 2 -o SAM '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -K mismatch:stats.txt -i MP 99-99 99,99",

            'single_basic': "/path/to/novoalign -d reference.fasta.idx -f NA10831_ATCACG_L002.fastq.gz   -c 2 -o SAM '@RG\\tID:NA10831_ATCACG_L002\\tSM:NA10831_ATCACG_L002' -K mismatch:stats.txt -i MP 99-99 99,99",

            'single_pipe': "/path/to/novoalign -d reference.fasta.idx -f 'NA|10831_ATCACG_L002.fastq.gz'   -c 2 -o SAM '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -K mismatch:stats.txt -i MP 99-99 99,99",
        }

        for sample_type, expect in tests.items():
            result = dispatcher._novoalign_command(self.novoalign.path, self.novoalign.args, 2, self.reference, *self.samples[sample_type])
            self.assertEqual(expect, result)


    def test_snap_command(self):
        tests = {
            'paired_basic': "/path/to/snap paired output_folder/reference/snap NA10831_ATCACG_L002_R1_001.fastq.gz NA10831_ATCACG_L002_R2_001.fastq.gz -t 2 -b --TODO -o sam -",

            'paired_pipe': "/path/to/snap paired output_folder/reference/snap 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -t 2 -b --TODO -o sam -",

            'single_basic': "/path/to/snap single output_folder/reference/snap NA10831_ATCACG_L002.fastq.gz  -t 2 -b --TODO -o sam -",

            'single_pipe': "/path/to/snap single output_folder/reference/snap 'NA|10831_ATCACG_L002.fastq.gz'  -t 2 -b --TODO -o sam -",
        }

        for sample_type, expect in tests.items():
            result = dispatcher._snap_command(self.snap.path, self.snap.args, 2, self.reference, self.output_folder, *self.samples[sample_type])
            self.assertEqual(expect, result)


    def test_align_reads(self):
        configuration = {
            'job_submitter': self.job_submitter,
            'aligners': [
                self.bowtie2,
                self.bwamem,
                self.bwa,
                self.novoalign,
                self.snap
            ],
            'samtools': ['samtools', '/path/to/samtools'],
            'output_folder': self.output_folder
        }

        expected_submit_job_calls = [
            call(
                'pbs',
                "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg 'SM:NA|10831_ATCACG_L002' --rg-id 'NA|10831_ATCACG_L002' -x reference -1 'NA|10831_ATCACG_L002_R1_001.fastq.gz' -2 'NA|10831_ATCACG_L002_R2_001.fastq.gz'",
                {'work_dir': 'output_folder/bowtie2', 'num_cpus': 2, 'name': 'nasp_bowtie2_NA|10831_ATCACG_L002'},
                (('jobid', 'action'),)
            ),
            call(
                'pbs',
                "/path/to/bwa mem -R '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0' -t 2 reference.fasta 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz'",
                {'name': 'nasp_bwamem_NA|10831_ATCACG_L002', 'num_cpus': 2, 'work_dir': 'output_folder/bwamem'},
                (('jobid', 'action'),)
            ),
            call(
                'pbs',
                "/path/to/bwa aln  reference.fasta 'NA|10831_ATCACG_L002_R1_001.fastq.gz' -t 2 -f 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa aln  reference.fasta 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -t 2 -f 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'; /path/to/bwa sampe -r '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' reference.fasta 'output_folder/bwa/NA|10831_ATCACG_L002-R1.sai' 'output_folder/bwa/NA|10831_ATCACG_L002-R2.sai' 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -x '-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0'",
                {'num_cpus': 2, 'work_dir': 'output_folder/bwa', 'name': 'nasp_bwa_NA|10831_ATCACG_L002'},
                (('jobid', 'action'),)
            ),
            call(
                'pbs',
                "/path/to/novoalign -d reference.fasta.idx -f 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -i PE 500,100 -c 2 -o SAM '@RG\\tID:NA|10831_ATCACG_L002\\tSM:NA|10831_ATCACG_L002' -K mismatch:stats.txt -i MP 99-99 99,99",
                {'num_cpus': 2, 'work_dir': 'output_folder/novo', 'name': 'nasp_novo_NA|10831_ATCACG_L002'},
                (('jobid', 'action'),)
            ),
            call(
                'pbs',
                "/path/to/snap paired output_folder/reference/snap 'NA|10831_ATCACG_L002_R1_001.fastq.gz' 'NA|10831_ATCACG_L002_R2_001.fastq.gz' -t 2 -b --TODO -o sam -",
                {'work_dir': 'output_folder/snap', 'num_cpus': 2, 'name': 'nasp_snap_NA|10831_ATCACG_L002'},
                (('jobid', 'action'),)
            )
        ]

        dispatcher._align_reads(self.samples['paired_pipe'], configuration, self.index_job_id, self.reference)
        dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)

