__author__ = 'jtravis'

import unittest
from unittest.mock import Mock, call
import os
import itertools
from tempfile import TemporaryDirectory

from nasp import dispatcher

# TODO: @unittest.skipUnless(termios, 'tests require system with termios')


class DispatcherShellEscapeCommandsTestCase(unittest.TestCase):
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

        # Minimal configuration for dispatcher.begin to complete
        self.configuration = {
                # 'output_folder': tmpdir,
                'index': ['', '', '', {}],
                'reference': ['', ''],
                'aligners': [],
                'snpcallers': [],
                'job_submitter': 'fake_job_submitter',
                'find_dups': 'False',
                'dup_finder': ('', '', '', {}),
                'assemblies': [],
                'alignments': [],
                'reads': [],
                'vcfs': [],
                # FIXME: It should crash or provide a reasonable default if path to matrix generator is unspecified.
                # ATM it returns a command with the correct arguments, but no executable/path.
                'matrix_generator': ['', '/path/to/vtm', '', {
                    # FIXME: It should default to at least one CPU if unspecified
                    'num_cpus': 0
                    }]
                }

        # Convenience references
        self.job_submitter = self.configuration['job_submitter']
        self.matrix_generator = self.configuration['matrix_generator'][1]

    def tearDown(self):
        pass

    def test_index_reference(self):
        pass

    def test_run_bwa(self):
        pass	

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
            )
        ]

        bam_nickname, job_id, final_file = dispatcher._run_bowtie2(("read1", "read2"), ["", "bowtie2", "", {'num_cpus': 42}], ["", "samtools"], "pbs", ("jobid", "action"), "reference", "output_folder")
        bam_nickname, job_id, final_file = dispatcher._run_bowtie2(("foo|bar",""), ["", "bowtie2", "", {'num_cpus': 42}], ["", "samtools"], "pbs", ("jobid", "action"), "reference", "output_folder")
        print("FOOBAR", bam_nickname, job_id, final_file)
        print(dispatcher)

        dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)

