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
        self.bowtie2 = App('bowtie2', '/path/to/bowtie2', '--very-sensitive-local --un pipe|in|name.fastq.gz --al "space in name.fastq.gz"', {'num_cpus': 2})
        self.samtools = App('samtools', '/path/to/samtools', '', {})
        self.job_submitter = 'pbs'
        self.reference = 'reference.fasta'
        self.output_folder = 'output_folder'
        self.index_job_id = ('jobid', 'action')

    def tearDown(self):
        pass

    def test_index_reference(self):
        pass

    def test_run_bwa(self):
        pass	

    def test_samtools_view_sort_index_pipe_command(self):
        # baseline test
        expect = "/path/to/samtools view -S -b -h - | /path/to/samtools sort - NA10831_ATCACG_L002; /path/to/samtools index NA10831_ATCACG_L002.bam"
        result = dispatcher._samtools_view_sort_index_command(self.samtools.path, self.samples['paired_basic'].name)
        self.assertEqual(expect, result)

        # special character test
        expect = "/path/to/samtools view -S -b -h - | /path/to/samtools sort - 'NA|10831_ATCACG_L002'; /path/to/samtools index 'NA|10831_ATCACG_L002.bam'"
        result = dispatcher._samtools_view_sort_index_command(self.samtools.path, self.samples['paired_pipe'].name)
        self.assertEqual(expect, result)

    def test_bowtie2_command(self):
        expect = "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg SM:NA10831_ATCACG_L002 --rg-id NA10831_ATCACG_L002 -x reference -1 NA10831_ATCACG_L002_R1_001.fastq.gz -2 NA10831_ATCACG_L002_R2_001.fastq.gz"
        result = dispatcher._bowtie2_command(self.bowtie2.path, self.bowtie2.args, 2, self.reference, *self.samples['paired_basic'])
        self.assertEqual(expect, result)

        expect = "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg 'SM:NA|10831_ATCACG_L002' --rg-id 'NA|10831_ATCACG_L002' -x reference -1 'NA|10831_ATCACG_L002_R1_001.fastq.gz' -2 'NA|10831_ATCACG_L002_R2_001.fastq.gz'"
        result = dispatcher._bowtie2_command(self.bowtie2.path, self.bowtie2.args, 2, self.reference, *self.samples['paired_pipe'])
        self.assertEqual(expect, result)

        expect = "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg SM:NA10831_ATCACG_L002 --rg-id NA10831_ATCACG_L002 -x reference -U NA10831_ATCACG_L002.fastq.gz"
        result = dispatcher._bowtie2_command(self.bowtie2.path, self.bowtie2.args, 2, self.reference, *self.samples['single_basic'])
        self.assertEqual(expect, result)

        expect = "/path/to/bowtie2 --very-sensitive-local --un 'pipe|in|name.fastq.gz' --al 'space in name.fastq.gz' --threads 2 --rg 'SM:NA|10831_ATCACG_L002' --rg-id 'NA|10831_ATCACG_L002' -x reference -U 'NA|10831_ATCACG_L002.fastq.gz'"
        result = dispatcher._bowtie2_command(self.bowtie2.path, self.bowtie2.args, 2, self.reference, *self.samples['single_pipe'])
        self.assertEqual(expect, result)

    def test_run_bowtie2(self):
        expected_submit_job_calls = [
            call(
                'pbs',
                "bowtie2  -p 42 --rg-id 'read1' --rg 'SM:read1' -x reference -U read2 | samtools view -S -b -h - | samtools sort - read1-bowtie2 \n a index read1-bowtie2.bam",
                {
                    'num_cpus': 42,
                    'name': 'nasp_bowtie2_read1',
                    'work_dir': 'output_folder/bowtie2'
                },
                (
                    ('jobid', 'action'),
                )
            ),
        ]

        bam_nickname, job_id, final_file = dispatcher._run_bowtie2(self.reads['paired_basic'], self.bowtie2, self.samtools, self.job_submitter, self.index_job_id, self.reference, self.output_folder)
        #bam_nickname, job_id, final_file = dispatcher._run_bowtie2(("foo|bar",""), ["", "bowtie2", "", {'num_cpus': 42}], ["", "samtools"], "pbs", ("jobid", "action"), "reference", "output_folder")
        print("FOOBAR", bam_nickname, job_id, final_file)
        print(dispatcher)

        dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)

