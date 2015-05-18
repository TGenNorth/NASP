__author__ = 'jtravis'

import unittest
from unittest.mock import Mock, call
import os
import itertools
from tempfile import TemporaryDirectory

from nasp import dispatcher

# TODO: @unittest.skipUnless(termios, 'tests require system with termios')


class DispatcherTestCase(unittest.TestCase):
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

    def test_it_fails_early_if_submit_job_fails(self):
        dispatcher._submit_job = Mock(return_value=None)
        with self.assertRaises(SystemExit), TemporaryDirectory() as tmpdir:
            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)

        dispatcher._submit_job.assert_called_once_with(
            self.job_submitter,
            'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
            {'work_dir': '{0}/reference'.format(tmpdir)},
            hold=True
        )
        assert not dispatcher._release_hold.called

    def test_submit_jobs_when_configuration_empty(self):
        expected_output_folder_files = ['matrices', 'matrix_dto.xml', 'reference', 'statistics']

        with TemporaryDirectory() as tmpdir:
            self.configuration['output_folder'] = tmpdir
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            dispatcher.begin(self.configuration)
            self.assertListEqual(expected_output_folder_files, os.listdir(tmpdir))
        dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
        dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_removes_old_reference_fasta_files_from_the_output_folder(self):
        with TemporaryDirectory() as tmpdir:
            self.configuration['output_folder'] = tmpdir
            reference_dir = os.path.join(tmpdir, 'reference')
            reference_fasta = os.path.join(tmpdir, 'reference', 'reference.fasta')
            os.mkdir(reference_dir)
            with open(reference_fasta, 'w+'):
                pass
            self.assertListEqual(['reference.fasta'], os.listdir(reference_dir))
            dispatcher.begin(self.configuration)
            self.assertListEqual([], os.listdir(os.path.dirname(reference_fasta)))

    def test_it_submits_index_references_jobs_for_the_bwa_aligner(self):
        self.configuration['aligners'] = [
            ('bwa', 'path/to/bwa')
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/bwa index {0}/reference/reference.fasta'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_submits_index_reference_jobs_for_the_bowtie2_aligner(self):
        self.configuration['aligners'] = [
            ('bowtie2', 'path/to/bowtie2'),
            ('bt2', 'path/to/bowtie2')
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/bowtie2-build {0}/reference/reference.fasta reference\n'
                    'path/to/bowtie2-build {0}/reference/reference.fasta reference'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_submits_index_reference_jobs_for_the_novo_aligner(self):
        self.configuration['aligners'] = [
            ('novo', 'path/to/novo'),
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/novoindex {0}/reference/reference.fasta.idx {0}/reference/reference.fasta'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_submits_index_reference_jobs_for_the_snap_aligner(self):
        self.configuration['aligners'] = [
            ('snap', 'path/to/snap'),
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/snap index {0}/reference/reference.fasta {0}/reference/snap'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_ignores_unsupported_aligners(self):
        self.configuration['aligners'] = [
            ('fake_aligner', 'path/to/fake_aligner'),
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    # It provides reasonable defaults, so they're not really required
    def test_creates_a_Sequence_Dictionary_and_samtools_index_of_the_reference_when_using_GATK_with_default_values(
            self):
        # FIXME: It crashes if the configuration does not contain a picard or samtools key
        self.configuration['picard'] = ['', '', '', '']
        self.configuration['samtools'] = ['', '']

        self.configuration['snpcallers'] = [
            ('GATK', 'path/to/gatk'),
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'java -Xmx2G -jar CreateSequenceDictionary.jar R={0}/reference/reference.fasta O={0}/reference/reference.dict\n'
                    'samtools faidx {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_creates_a_Sequence_Dictionary_and_samtools_index_of_the_reference_when_using_GATK(self):
        # TODO: fill picard and samtool to check it does not use default values when values are specified
        self.configuration['picard'] = ['', '', '', '']
        self.configuration['samtools'] = ['', '']

        self.configuration['snpcallers'] = [
            ('GATK', 'path/to/gatk'),
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'java -Xmx2G -jar CreateSequenceDictionary.jar R={0}/reference/reference.fasta O={0}/reference/reference.dict\n'
                    'samtools faidx {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_checks_for_duplicate_positions(self):
        self.configuration["find_dups"] = "True"

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _find_dups
                call(
                    self.job_submitter,
                    'find_duplicates --nucmerpath  --reference {0}/reference/reference.fasta'.format(tmpdir),
                    {'name': 'nasp_', 'work_dir': '{0}/reference'.format(tmpdir)},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '{0} --mode xml --dto-file {1}/matrix_dto.xml --num-threads 0'.format(self.matrix_generator,
                                                                                          tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('2', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_normalizes_external_fastas(self):
        # FIXME: It crashes if the configuration does not contain an assembly_importer key
        self.configuration['assembly_importer'] = ['', '', '', {}]

        self.configuration['assemblies'] = [
            ['external', 'path/to/external.fasta']
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _convert_external_genome
                call(
                    self.job_submitter,
                    "format_fasta --inputfasta path/to/external.fasta --outputfasta {0}/external/external.fasta\n"
                    "convert_external_genome --nucmerpath  --nucmerargs '' --deltafilterpath  --deltafilterargs '' --reference {0}/reference/reference.fasta --external path/to/external.fasta --name external".format(
                        tmpdir),
                    {'work_dir': '{0}/external'.format(tmpdir), 'name': 'nasp__external'},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('2', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_indexes_pre_aligned_bam_files(self):
        # FIXME: It crashes if the configuration does not contain a bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\n'
                    '{1} index {0}/bams/name.bam'.format(tmpdir, self.configuration['samtools'][1]),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_snpcalls_pre_aligned_bam_files_with_gatk(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['picard'] = ['', '', '', {}]
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]
        self.configuration['snpcallers'] = [
            # FIXME: It should use reasonable defaults if snpcaller parameters are undefined
            ['gatk', 'path/to/gatk', '-args', {
                'num_cpus': 0,
                'mem_requested': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'java -Xmx2G -jar CreateSequenceDictionary.jar R={0}/reference/reference.fasta O={0}/reference/reference.dict\n'
                    'path/to/samtools faidx {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\npath/to/samtools index {0}/bams/name.bam'.format(tmpdir),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _call_snps
                call(
                    self.job_submitter,
                    'java -Xmx1G -jar path/to/gatk -T UnifiedGenotyper -dt NONE -glm BOTH -I {0}/bams/name.bam -R {0}/reference/reference.fasta -nt 0 -o name-gatk.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE -args'.format(
                        tmpdir),
                    {'num_cpus': 0, 'name': 'nasp_gatk_name', 'mem_requested': 1,
                     'work_dir': '{0}/gatk'.format(tmpdir)},
                    ('2',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'num_cpus': 0, 'work_dir': tmpdir},
                    ('3', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_snpcalls_pre_aligned_bam_files_with_solsnp(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['picard'] = ['', '', '', {}]
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]
        self.configuration['snpcallers'] = [
            # FIXME: It should use reasonable defaults if snpcaller parameters are undefined
            ['solsnp', 'path/to/solsnp', '-args', {
                'mem_requested': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\n'
                    'path/to/samtools index {0}/bams/name.bam'.format(tmpdir),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _call_snps
                call(
                    self.job_submitter,
                    'java -Xmx1G -jar path/to/solsnp INPUT={0}/solsnp/name REFERENCE_SEQUENCE={0}/reference/reference.fasta OUTPUT={0}/solsnp/name-solsnp.vcf SUMMARY=true CALCULATE_ALLELIC_BALANCE=true MINIMUM_COVERAGE=1 PLOIDY=Haploid STRAND_MODE=None OUTPUT_FORMAT=VCF OUTPUT_MODE=AllCallable -args'.format(
                        tmpdir),
                    {'mem_requested': 1,
                     'work_dir': '{0}/solsnp'.format(tmpdir),
                     'name': 'nasp_solsnp_name'},
                    ('2',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('3', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_snpcalls_pre_aligned_bam_files_with_varscan(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['picard'] = ['', '', '', {}]
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]
        self.configuration['snpcallers'] = [
            # FIXME: It should use reasonable defaults if snpcaller parameters are undefined
            ['varscan', 'path/to/varscan', '-args', {
                'mem_requested': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\npath/to/samtools index {0}/bams/name.bam'.format(tmpdir),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _call_snps
                call(
                    self.job_submitter,
                    'echo name > {0}/varscan/name.txt\npath/to/samtools mpileup -B -d 10000000 -f {0}/reference/reference.fasta {0}/bams/name.bam > {0}/bams/name.mpileup'
                    '\njava -Xmx1G -jar path/to/varscan mpileup2cns {0}/bams/name.mpileup --output-vcf 1 --vcf-sample-list {0}/varscan/name.txt > {0}/varscan/name-varscan.vcf -args'.format(
                        tmpdir),
                    {'mem_requested': 1, 'name': 'nasp_varscan_name',
                     'work_dir': '{0}/varscan'.format(tmpdir)},
                    ('2',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'num_cpus': 0, 'work_dir': '{0}'.format(tmpdir)},
                    ('3', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_snpcalls_with_samtools(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['picard'] = ['', '', '', {}]
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]
        self.configuration['snpcallers'] = [
            # FIXME: It should use reasonable defaults if snpcaller parameters are undefined
            ['samtools', 'path/to/samtools', '-args', {
                'mem_requested': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\n'
                    'path/to/samtools index {0}/bams/name.bam'.format(tmpdir),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _call_snps
                call(
                    self.job_submitter,
                    'path/to/samtools mpileup -uD -d 10000000 -f {0}/reference/reference.fasta {0}/bams/name.bam | path/to/samtools view -ceg -args - > {0}/samtools/name-samtools.vcf'.format(
                        tmpdir),
                    {'name': 'nasp_samtools_name', 'work_dir': '{0}/samtools'.format(tmpdir), 'mem_requested': 1},
                    ('2',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('3', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_ignores_unsupported_snpcallers(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['picard'] = ['', '', '', {}]
        self.configuration['samtools'] = ['', 'path/to/samtools']
        self.configuration['bam_index'] = ['', 'path/to/bam_index', '', {}]

        self.configuration['alignments'] = [
            ['name', 'path/to/bam']
        ]
        self.configuration['snpcallers'] = [
            # FIXME: It should use reasonable defaults if snpcaller parameters are undefined
            ['fake_snpcaller', 'path/to/fake_snpcaller', '-args', {
                'mem_requested': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _index_bams
                call(
                    self.job_submitter,
                    'ln -s -f path/to/bam {0}/bams/name.bam\n'
                    'path/to/samtools index {0}/bams/name.bam'.format(tmpdir),
                    {'work_dir': '{0}/bams'.format(tmpdir)},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'),
                    notify=True
                )
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_aligns_reads_with_bwa(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']

        self.configuration['reads'] = [
            ['name', 'path/to/fasta']
        ]
        self.configuration['aligners'] = [
            # FIXME: It should use reasonable defaults if job parameters are undefined
            ['bwa', 'path/to/bwa', '-args', {
                'num_cpus': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/bwa index {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True),
                # _align_reads
                call(
                    self.job_submitter,
                    "path/to/bwa aln  {0}/reference/reference.fasta path/to/fasta -t 1 -f {0}/bwa/name.sai -args\n"
                    "path/to/bwa samse -r '@RG\\tID:name\\tSM:name' {0}/reference/reference.fasta {0}/bwa/name.sai path/to/fasta -args | {1} view -S -b -h - | {1} sort - name-bwa \n"
                    " {1} index name-bwa.bam".format(tmpdir, self.configuration['samtools'][1]),
                    {'work_dir': '{0}/bwa'.format(tmpdir), 'num_cpus': 1,
                     'name': 'nasp_bwa_name'}, ('1',)),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'), notify=True)
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_aligns_reads_with_bowtie2(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']

        self.configuration['reads'] = [
            ['name', 'path/to/fasta']
        ]
        self.configuration['aligners'] = [
            # FIXME: It should use reasonable defaults if job parameters are undefined
            ['bt2', 'path/to/bt2', '-args', {
                'num_cpus': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/bowtie2-build {0}/reference/reference.fasta reference'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True),
                # _align_reads
                call(
                    self.job_submitter,
                    "path/to/bt2 -args -p 1 --rg-id 'name' --rg 'SM:name' -x {0}/reference/reference -U path/to/fasta | path/to/samtools view -S -b -h - | path/to/samtools sort - name-bowtie2"
                    " \n path/to/samtools index name-bowtie2.bam".format(tmpdir),
                    {'num_cpus': 1, 'work_dir': '{0}/bowtie2'.format(tmpdir),
                     'name': 'nasp_bowtie2_name'},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'), notify=True)
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_aligns_reads_with_novo(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']

        self.configuration['reads'] = [
            ['name', 'path/to/fasta']
        ]
        self.configuration['aligners'] = [
            # FIXME: It should use reasonable defaults if job parameters are undefined
            ['novo', 'path/to/novo', '-args', {
                'num_cpus': 1
            }],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/novoindex {0}/reference/reference.fasta.idx {0}/reference/reference.fasta'.format(tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True),
                # _align_reads
                call(
                    self.job_submitter,
                    "path/to/novo -f path/to/fasta   -c 1 -o SAM '@RG\\tID:name\\tSM:name' -d {0}/reference/reference.fasta.idx -args | path/to/samtools view -S -b -h - | path/to/samtools sort - name-novo"
                    " \n path/to/samtools index name-novo.bam".format(tmpdir),
                    {'name': 'nasp_novo_name', 'num_cpus': 1,
                     'work_dir': '{0}/novo'.format(tmpdir)},
                    ('1',)
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'), notify=True)
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_aligns_reads_with_snap(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']

        self.configuration['reads'] = [
            ['name', 'path/to/fasta']
        ]
        self.configuration['aligners'] = [
            # FIXME: It should use reasonable defaults if job parameters are undefined
            ['snap', 'path/to/snap', '-args', {}],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta\n'
                    'path/to/snap index {0}/reference/reference.fasta {0}/reference/snap'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _align_reads
                call(
                    self.job_submitter,
                    ' \n path/to/snap single {0}/reference/snap path/to/fasta  -o name-snap.sam -args \n'
                    ' path/to/samtools view -S -b -h name-snap.sam | path/to/samtools sort - name-snap \n'
                    ' path/to/samtools index name-snap.bam \n '.format(tmpdir),
                    {'work_dir': '{0}/snap'.format(tmpdir),
                     'name': 'nasp_snap_name'}, ('1',)),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'), notify=True)
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])

    def test_it_ignores_unknown_aligners(self):
        # FIXME: It crashes if the configuration does not contain a picard, samtools, or bam_index key
        self.configuration['samtools'] = ['', 'path/to/samtools']

        self.configuration['reads'] = [
            ['name', 'path/to/fasta']
        ]
        self.configuration['aligners'] = [
            # FIXME: It should use reasonable defaults if job parameters are undefined
            ['fake_aligner', 'path/to/fake_aligner', '-args', {}],
        ]

        with TemporaryDirectory() as tmpdir:
            expected_submit_job_calls = [
                # _index_reference
                call(
                    self.job_submitter,
                    'format_fasta --inputfasta  --outputfasta {0}/reference/reference.fasta'.format(
                        tmpdir),
                    {'work_dir': '{0}/reference'.format(tmpdir)},
                    hold=True
                ),
                # _create_matrices
                call(
                    self.job_submitter,
                    '/path/to/vtm --mode xml --dto-file {0}/matrix_dto.xml --num-threads 0'.format(tmpdir),
                    {'work_dir': tmpdir, 'num_cpus': 0},
                    ('', 'afterany'), notify=True)
            ]

            self.configuration['output_folder'] = tmpdir
            dispatcher.begin(self.configuration)
            dispatcher._submit_job.assert_has_calls(expected_submit_job_calls)
            dispatcher._release_hold.assert_has_calls([call(self.job_submitter, '1')])


class DispatcherNoJobManagerSubmitJobTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @unittest.mock.patch('subprocess.getoutput', return_value=42)
    @unittest.mock.patch('subprocess.Popen')
    @unittest.mock.patch('os.open', autospec=True)
    def test_it_manages_jobs_when_no_supported_job_manager_is_available(self, mock_getoutput, mock_subproc_popen, mock_open):
        self.job_submitter = 'fake_job_submitter'
        self.command = 'fake_command'
        self.job_parms = {
            'work_dir': 'job_parms',
            'mem_requested': -1,
            'name': 'fake_name'
        }
        self.waitfor_id = ('fake', 'id')

        expected = '43'

        process_mock = unittest.mock.Mock()
        attrs = {'pid': expected}
        process_mock.configure_mock(**attrs)
        mock_subproc_popen.return_value = process_mock


        mock_getoutput.return_value = expected

        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)
        self.assertEqual(expected, observed)
        self.assertTrue(mock_getoutput.called)
        self.assertTrue(mock_subproc_popen.called)
        self.assertTrue(mock_open.called)

class DispatcherPbsSubmitJobTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.job_submitter = 'PBS'
        self.command = 'fake_command'
        self.job_parms = {
            'queue': 'fake_queue',
            'args': '-args',
            'work_dir': 'my_work_dir',
            'num_cpus': -1,
            'mem_requested': -2,
            'walltime': -3,
            'name': 'fake_name'
        }

        self.waitfor_id = ('-4', 'afterok')
        self.expected_command = {
            'command': self.command,
            'work_dir': self.job_parms['work_dir'],
            'ncpus': self.job_parms['num_cpus'],
            'mem': self.job_parms['mem_requested'],
            'walltime': self.job_parms['walltime'],
            'name': self.job_parms['name'],
            'dependency_string': self.waitfor_id[1],
            'job_id': self.waitfor_id[0],
            'queue': self.job_parms['queue'],
            'args': self.job_parms['args']
        }

    def tearDown(self):
        pass


    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_submits_PBS_jobs(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -d \'{work_dir}\' -w \'{work_dir}\' -l ncpus={ncpus},mem={mem}gb,walltime={walltime}:00:00 -m a -N \'fake_name\' -W depend={dependency_string}:{job_id} -q {queue} {args} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='')
    def test_it_logs_an_error_if_PBS_does_not_return_a_job_id(self, mock_getoutput):

        expected = None
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -d \'{work_dir}\' -w \'{work_dir}\' -l ncpus={ncpus},mem={mem}gb,walltime={walltime}:00:00 -m a -N \'fake_name\' -W depend={dependency_string}:{job_id} -q {queue} {args} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_sets_hold_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, hold=True)

        self.assertEqual(expected, observed)

        self.expected_command['hold_flag'] = '-h'

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -d \'{work_dir}\' -w \'{work_dir}\' -l ncpus={ncpus},mem={mem}gb,walltime={walltime}:00:00 -m a -N \'fake_name\' -W depend={dependency_string}:{job_id} -q {queue} {args} {hold_flag} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_sets_notify_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, notify=True)

        self.assertEqual(expected, observed)

        self.expected_command['notify_flag'] = '-m e'

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -d \'{work_dir}\' -w \'{work_dir}\' -l ncpus={ncpus},mem={mem}gb,walltime={walltime}:00:00 -m a -N \'fake_name\' -W depend={dependency_string}:{job_id} -q {queue} {args} {notify_flag} - '.format(**self.expected_command))

class DispatcherSlurmSubmitJobTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.job_submitter = 'SLURM'
        self.command = 'fake_command'
        self.job_parms = {
            'queue': 'fake_queue',
            'args': '-args',
            'work_dir': 'my_work_dir',
            'num_cpus': -1,
            'mem_requested': -2,
            'walltime': -3,
            'name': 'fake_name'
        }
        self.waitfor_id = ('-4', 'afterok')

        self.expected_command = {
            'command': self.command,
            'work_dir': self.job_parms['work_dir'],
            'ncpus': self.job_parms['num_cpus'],
            'mem': self.job_parms['mem_requested'],
            'walltime': self.job_parms['walltime'],
            'name': self.job_parms['name'],
            'dependency_string': self.waitfor_id[1],
            'job_id': self.waitfor_id[0],
            'queue': self.job_parms['queue'],
            'args': self.job_parms['args']
        }

    def tearDown(self):
        pass


    @unittest.mock.patch('subprocess.getoutput', return_value='Submitted batch job 12345')
    def test_it_submits_jobs(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        # self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('sbatch -D \'{work_dir}\' -c{ncpus} --mem={mem}000 --time={walltime}:00:00 --mail-type=FAIL -J \'{name}\' -d {dependency_string}:{job_id} -p {queue} -args --wrap="{command}"'.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='')
    def test_it_logs_an_error_if_it_does_not_return_a_job_id(self, mock_getoutput):

        expected = None
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('sbatch -D \'{work_dir}\' -c{ncpus} --mem={mem}000 --time={walltime}:00:00 --mail-type=FAIL -J \'{name}\' -d {dependency_string}:{job_id} -p {queue} -args --wrap="{command}"'.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='Submitted batch job 12345')
    def test_it_sets_hold_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, hold=True)

        self.assertEqual(expected, observed)

        self.expected_command['hold_flag'] = '-H'

        mock_getoutput.assert_called_once_with('sbatch -D \'{work_dir}\' -c{ncpus} --mem={mem}000 --time={walltime}:00:00 --mail-type=FAIL -J \'{name}\' -d {dependency_string}:{job_id} -p {queue} -args {hold_flag} --wrap="{command}"'.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='Submitted batch job 12345')
    def test_it_set_notify_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, notify=True)

        self.assertEqual(expected, observed)

        self.expected_command['notify_flag'] = '--mail-type=END'

        mock_getoutput.assert_called_once_with('sbatch -D \'{work_dir}\' -c{ncpus} --mem={mem}000 --time={walltime}:00:00 --mail-type=FAIL -J \'{name}\' -d {dependency_string}:{job_id} -p {queue} -args {notify_flag} --wrap="{command}"'.format(**self.expected_command))

class DispatcherSGESubmitJobTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.job_submitter = 'SGE'
        self.command = 'fake_command'
        self.job_parms = {
            'queue': 'fake_queue',
            'args': '-args',
            'work_dir': 'my_work_dir',
            'num_cpus': -1,
            'mem_requested': -2,
            'walltime': -3,
            'name': 'fake_name'
        }
        self.waitfor_id = ('-4', 'afterok')

        self.expected_command = {
            'command': self.command,
            'work_dir': self.job_parms['work_dir'],
            'ncpus': self.job_parms['num_cpus'],
            'mem': self.job_parms['mem_requested'],
            'walltime': self.job_parms['walltime'],
            'name': self.job_parms['name'],
            'dependency_string': self.waitfor_id[1],
            'job_id': self.waitfor_id[0],
            'queue': self.job_parms['queue'],
            'args': self.job_parms['args']
        }

    def tearDown(self):
        pass


    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_submits_jobs(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        # self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -cwd \'{work_dir}\' -wd \'{work_dir}\' -l h_data=-2097152.0gb,h_rt=-3:00:00 -m a -N \'{name}\' -hold_jid -4 -q {queue} {args} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='')
    def test_it_logs_an_error_if_it_does_not_return_a_job_id(self, mock_getoutput):

        expected = None
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id)

        self.assertEqual(expected, observed)

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -cwd \'{work_dir}\' -wd \'{work_dir}\' -l h_data=-2097152.0gb,h_rt=-3:00:00 -m a -N \'{name}\' -hold_jid -4 -q {queue} {args} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_sets_hold_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, hold=True)

        self.assertEqual(expected, observed)

        self.expected_command['hold_flag'] = '-h'

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -cwd \'{work_dir}\' -wd \'{work_dir}\' -l h_data=-2097152.0gb,h_rt=-3:00:00 -m a -N \'{name}\' -hold_jid -4 -q {queue} {args} {hold_flag} - '.format(**self.expected_command))

    @unittest.mock.patch('subprocess.getoutput', return_value='12345.host')
    def test_it_set_notify_flags(self, mock_getoutput):

        expected = '12345'
        observed = dispatcher._submit_job(self.job_submitter, self.command, self.job_parms, self.waitfor_id, notify=True)

        self.assertEqual(expected, observed)

        self.expected_command['notify_flag'] = '-m e'

        mock_getoutput.assert_called_once_with('echo "{command}" | qsub -V -cwd \'{work_dir}\' -wd \'{work_dir}\' -l h_data=-2097152.0gb,h_rt=-3:00:00 -m a -N \'{name}\' -hold_jid -4 -q {queue} {args} {notify_flag} - '.format(**self.expected_command))
