#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "1.0.0"
__email__ = "dlemmer@tgen.org"

'''
Created on Mar 4, 2014

@author: dlemmer
'''

import logging
import shlex
from collections import namedtuple
import os

App = namedtuple('App', ['name', 'path', 'args', 'job_params'])
#Assembly = namedtuple('Assembly', ['name', 'read1', 'read2'])

def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--config", required=True, help="Path to the configuration xml file.")
    return parser.parse_args()


def _pbs_command(name, work_dir, mem_requested=1, num_cpus=1, walltime=1, queue='', args='', hold=False, notify=False, waitfor_id=None):

    job_resources = 'ncpus={ncpu},mem={mem}gb,walltime={hours}:00:00'.format(**{
        'ncpu': num_cpus,
        'mem': mem_requested,
        'hours': walltime
    })

    pbs_cmd = shlex.split("qsub -V -d {work_dir} -w {work_dir} -l {job_resources} -m a -N {job_name}".format(**{
        'work_dir': shlex.quote(work_dir),
        'job_resources': shlex.quote(job_resources),
        'job_name': shlex.quote(name),
    }))

    if waitfor_id:
        pbs_cmd.extend(['-W', 'depend={directive}:{job_ids}'.format(**{
            'directive': waitfor_id[1] if len(waitfor_id) > 1 else 'afterok',
            'job_ids': waitfor_id[0],
        })])
    if queue:
        pbs_cmd.extend(['-q', queue])
    if hold:
        pbs_cmd.append('-h')
    if notify:
        pbs_cmd.extend(['-m', 'e'])

    pbs_cmd.extend(shlex.split(args))

    return ' '.join(map(shlex.quote, pbs_cmd))

    logging.debug("submit_command = %s", submit_command)
    output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
    logging.debug("output = %s" % output)
    job_match = re.search('^(\d+)\..*$', output)
    if job_match:
        jobid = job_match.group(1)
    else:
        logging.warning("Job not submitted!!")
        print("WARNING: Job not submitted: %s" % output)


def _slurm_command(name, work_dir, mem_requested=1, num_cpus=1, walltime=1, queue='', args='', hold=False, notify=False, waitfor_id=None):
    slurm_cmd = shlex.split('sbatch -D {work_dir} -c {ncpu} --mem={mem_gb} --time={hours} --mail-type=FAIL -J {job_name}'.format(**{
        'work_dir': shlex.quote(work_dir),
        'ncpu': shlex.quote(str(num_cpus)),
        'mem_gb': shlex.quote(str(mem_requested) + '000'),
        'hours': shlex.quote(str(walltime) + ':00:00'),
        'job_name': shlex.quote(name),
    }))

    if waitfor_id:
        slurm_cmd.extend([
            '-d', '{directive}:{job_ids}'.format(**{
                'directive': waitfor_id[1] if len(waitfor_id) > 1 else 'afterok',
                'job_ids': waitfor_id[0],
            })
        ])

    if queue:
        slurm_cmd.extend(['-p', queue])
    if hold:
        slurm_cmd.append('-H') 
    if notify:
        slurm_cmd.append('--mail-type=END')

    slurm_cmd.extend(shlex.split(args))

    return ' '.join(map(shlex.quote, slurm_cmd))

    logging.debug("submit_command = %s" % submit_command)
    output = subprocess.getoutput("{0} --wrap={1}".format(submit_command, shlex.quote(command)))
    logging.debug("output = {0}".format(output))
    job_match = re.search('^Submitted batch job (\d+)$', output)
    if job_match:
        jobid = job_match.group(1)
    else:
        logging.warning("Job not submitted!!")
        print("WARNING: Job not submitted: %s" % output)


def _submit_job(job_submitter, command, job_parms, waitfor_id=None, hold=False, notify=False):
    import subprocess
    import re
    import os

    jobid = None
    logging.info("command = %s" % command)
    if job_submitter == "PBS":
        submit_command = _pbs_command(**job_params, hold=hold, notify=notify, waitfor_id=waitfor_id)
        logging.debug("submit_command = {0}".format(submit_command))
        output = subprocess.getoutput("echo {0} | {1} - ".format(shlex.quote(command), submit_command))
        logging.debug("output = {0}".format(output))
        job_match = re.search('^(\d+)\..*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SLURM":
        submit_command = _slurm_command(**job_params, hold=hold, notify=notify, waitfor_id=waitfor_id)
        logging.debug("submit_command = %s" % submit_command)
        output = subprocess.getoutput("%s --wrap=\"%s\"" % (submit_command, command))
        logging.debug("output = %s" % output)
        job_match = re.search('^Submitted batch job (\d+)$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SGE":
        command = re.sub('\n', '; ', command)
        waitfor = ""
        if waitfor_id:
            waitfor = "-hold_jid %s" % (re.sub(":", ",", waitfor_id[0]))
        queue = ""
        if job_parms["queue"]:
            queue = "-q %s" % job_parms["queue"]
        args = job_parms["args"]
        if hold:
            args = "-h " + args
        if notify:
            args = "-m e " + args
        mem_needed = float(job_parms['mem_requested']) * 1024 * 1024
        # Apparently the number of processors a job uses is controlled by the queue it is running on in SGE, so there is no way to request a specific number of CPUs??
        submit_command = "qsub -V -wd \'%s\' -l h_data=%s,h_rt=%s:00:00 -m a -N \'%s\' -b y %s %s %s" % (
            job_parms["work_dir"], mem_needed, job_parms['walltime'], job_parms['name'], waitfor,
            queue, args)
        logging.debug("submit_command = %s", submit_command)
        #output = subprocess.getoutput("echo \"%s\" | %s" % (command, submit_command))
        output = subprocess.getoutput("%s \"%s\"" % (submit_command, command))
        logging.debug("output = %s" % output)
        job_match = re.search('^\D*(\d+)\s.*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    else:
        command = re.sub('\n', '; ', command)
        work_dir = job_parms['work_dir']
        dependency_check = ""
        if waitfor_id:
            if re.search(":", waitfor_id[0]):
                pid_filename = os.path.join(work_dir, "%s_dependent_pids" % job_parms['name'])
                pid_file = open(pid_filename, 'w')
                pids = waitfor_id[0].split(":")
                pid_file.write("\n".join(pids))
                pid_file.close()
                dependency_check = "while [ -s %s ]; do sleep 600; for pid in `cat %s`; do kill -0 \"$pid\" 2>/dev/null || sed -i \"/^$pid$/d\" %s; done; done; rm %s; " % (
                    pid_filename, pid_filename, pid_filename, pid_filename)
            else:
                pid = waitfor_id[0]
                dependency_check = "while kill -0 %s; do sleep 300; done; " % pid
        # cpu_check = "grep 'cpu ' /proc/stat | awk '{usage=($2+$4)/($2+$4+$5)} END {print 1-usage}'"
        total_mem = subprocess.getoutput("free -g | grep Mem: | awk '{ print $2 }'")
        mem_requested = job_parms['mem_requested']
        if mem_requested > total_mem:
            mem_requested = total_mem
        print("total_mem = %s" % total_mem)
        mem_needed = float(mem_requested) * 750
        mem_check = "while [ `free -m | grep cache: | awk '{ print $4 }'` -lt %d ]; do sleep 300; done; " % mem_needed
        submit_command = "%s%s%s" % (dependency_check, mem_check, command)
        logging.debug("submit_command = %s" % submit_command)
        output_log = os.open(os.path.join(work_dir, "%s.out" % job_parms['name']), os.O_WRONLY | os.O_CREAT)
        proc = subprocess.Popen(submit_command, stderr=subprocess.STDOUT, stdout=output_log, shell=True, cwd=work_dir)
        jobid = str(proc.pid)
    logging.info("jobid = %s" % jobid)
    return jobid


def _release_hold(job_submitter, job_id):
    import subprocess

    if job_submitter == "PBS" or job_submitter == "SGE":
        command = "qrls %s" % job_id
    elif job_submitter == "SLURM":
        command = "scontrol release %s" % job_id
    else:
        return
    logging.info("command = %s", command)
    output = subprocess.getoutput(command)
    logging.debug("output = %s", output)


def _index_reference(configuration):
    import os
    import re

    output_folder = configuration["output_folder"]
    job_parms = configuration["index"][3]
    ref_path = configuration["reference"][1]
    ref_folder = os.path.join(output_folder, "reference")
    if not os.path.exists(ref_folder):
        os.makedirs(ref_folder)
    # Copy the reference as $output_folder/reference/reference.fasta, verifying its format first. Replace it if it already exists.
    reference = os.path.join(ref_folder, "reference.fasta")
    if os.path.exists(reference):
        os.remove(reference)
    index_commands = ["format_fasta --inputfasta %s --outputfasta %s" % (ref_path, reference)]

    # Gather all of the index commands that need to be run
    bwa_done = False
    for aligner in configuration["aligners"]:
        (name, path) = aligner[0:2]
        if re.search('bwa', name, re.IGNORECASE):
            if not bwa_done:
                index_commands.append("%s index %s" % (path, reference))
                bwa_done = True
        elif re.search('b(ow)?t(ie)?2', name, re.IGNORECASE):
            bt2path = os.path.split(path)[0]
            bt2_build_path = os.path.join(bt2path, "bowtie2-build")
            index_commands.append("%s %s reference" % (bt2_build_path, reference))
        elif re.search('novo', name, re.IGNORECASE):
            novopath = os.path.split(path)[0]
            novoindex_path = os.path.join(novopath, "novoindex")
            index_commands.append("%s %s.idx %s" % (novoindex_path, reference, reference))
        elif re.search('snap', name, re.IGNORECASE):
            index_commands.append("%s index %s %s" % (path, reference, os.path.join(ref_folder, "snap")))
        else:
            print("Unknown aligner \'%s\' found, don't know how to index the reference for it. Skipping..." % name)

    #if we are using GATK, we also need to create a Sequence Dictionary and samtools index of the reference
    if next((v for i, v in enumerate(configuration["snpcallers"]) if re.search('gatk', v[0], re.IGNORECASE)), None):
        #picard_path = configuration["picard"][1] or ""
        picard_path = configuration["picard"][1]
        picard_memory = 2
        if configuration["picard"][3]:
            picard_memory = configuration["picard"][3]['mem_requested'] or 2
        samtools_path = configuration["samtools"][1] or "samtools"
        out_file = os.path.join(ref_folder, "reference.dict")
        index_commands.append("java -Xmx%sG -jar %s CreateSequenceDictionary R=%s O=%s" % (picard_memory, picard_path, reference, out_file))
        index_commands.append("%s faidx %s" % (samtools_path, reference))

    command = "\n".join(index_commands)
    job_parms['work_dir'] = ref_folder
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, hold=True)
    return job_id, reference


def _samtools_view_sort_index_pipe_command(samtools_path, output_bam):
    (bam_prefix, _) = os.path.splitext(output_bam)
    return '{samtools} view -S -b -h - | {samtools} sort - {bam_prefix}; {samtools} index {bam_filename}'.format(**{
        'samtools': samtools_path,
        'bam_prefix': shlex.quote(bam_prefix),
        'bam_filename': shlex.quote(output_bam)
    })


def _bwamem_command(path, args, ncpu, reference, output_folder, sample_name, read1, read2=None):
    """
    Args:
        path (str): path to aligner executable
        args (str): raw arguments to be passed to the aligner
        ncpu: number of alignment threads to launch
        reference: (str): reference filename
        sample_name (str): 
        read1 (str): absolute path to read1 fastq[.gz|.bz2]
        read2 (str): absolute path to read2 fastq[.gz|.bz2]

    Returns:
        string: command to execute bowtie2 aligner
    """
    return '{bwa} mem -R {bam_string} {bam_args} -t {ncpu} {reference} {read1} {read2}'.format(**{
        'bwa': path,
        'bam_string': shlex.quote('@RG\\tID:{sample_name}\\tSM:{sample_name}'.format(sample_name=sample_name)),
        'bam_args': ' '.join(map(shlex.quote, shlex.split(args))),
        'ncpu': shlex.quote(str(ncpu)),
        'reference': shlex.quote(reference),
        'read1': shlex.quote(read1),
        'read2': shlex.quote(read2) if read2 else ''
    })


def _bwa_command(path, args, ncpu, reference, output_folder, sample_name, read1, read2=None):
    """
    Args:
        path (str): path to aligner executable
        args (str): raw arguments to be passed to the aligner
        ncpu: number of alignment threads to launch
        reference: (str): reference filename
        output_folder (str): directory for aligner to write output
        sample_name (str): 
        read1 (str): absolute path to read1 fastq[.gz|.bz2]
        read2 (str): absolute path to read2 fastq[.gz|.bz2]

    Returns:
        string: command to execute bowtie2 aligner
    """
    import re
    import os

    bam_string = '@RG\\tID:{sample_name}\\tSM:{sample_name}'.format(sample_name=sample_name)
    # Parse read file basename
    is_illumina_fastq = "-I" if re.search('(?:.*\/)?[^\/]+?_[12]_sequence\.txt(?:\.gz)?$', read1, re.IGNORECASE) else ""
    quoted_bwa_args = ' '.join(map(shlex.quote, shlex.split(args)))
    work_dir = os.path.join(output_folder, 'bwa')
    if read2:
        outfile1 = os.path.join(work_dir, "{0}-R1.sai".format(sample_name))
        outfile2 = os.path.join(work_dir, "{0}-R2.sai".format(sample_name))
        align_read1 = '{bwa} aln {is_illumina_fastq} {reference} {read1} -t {ncpu} -f {outfile1} {bwa_args}'.format(**{
            'bwa': path,
            'is_illumina_fastq': is_illumina_fastq,
            'reference': shlex.quote(reference),
            'read1': shlex.quote(read1),
            'ncpu': shlex.quote(str(ncpu)),
            'outfile1': shlex.quote(outfile1),
            'bwa_args': quoted_bwa_args
        })
        align_read2 = '{bwa} aln {is_illumina_fastq} {reference} {read2} -t {ncpu} -f {outfile1} {bwa_args}'.format(**{
            'bwa': path,
            'is_illumina_fastq': is_illumina_fastq,
            'reference': shlex.quote(reference),
            'read2': shlex.quote(read2),
            'ncpu': shlex.quote(str(ncpu)),
            'outfile1': shlex.quote(outfile1),
            'bwa_args': quoted_bwa_args
        })
        sampe_command = '{bwa} sampe -r {bam_string} {reference} {outfile1} {outfile2} {read1} {read2} {bwa_args}'.format(**{
            'bwa': path,
            'bam_string': shlex.quote(bam_string),
            'reference': shlex.quote(reference),
            'outfile1': shlex.quote(outfile1),
            'outfile2': shlex.quote(outfile2),
            'read1': shlex.quote(read1),
            'read2': shlex.quote(read2),
            'bwa_args': quoted_bwa_args
        })
        aligner_command = '; '.join([align_read1, align_read2, sampe_command])
    else:
        outfile = os.path.join(work_dir, "%s.sai" % sample_name)
        align_read = '{bwa} aln {is_illumina_fastq} {reference} {read1} -t {ncpu} -f {outfile} {bwa_args}'.format(**{
            'bwa': path,
            'is_illumina_fastq': is_illumina_fastq,
            'reference': shlex.quote(reference),
            'read1': shlex.quote(read1),
            'ncpu': shlex.quote(str(ncpu)),
            'outfile': shlex.quote(outfile),
            'bwa_args': quoted_bwa_args
        })
        samse_command = '{bwa} samse -r {bam_string} {reference} {outfile} {read1} {bwa_args}'.format(**{
            'bwa': path,
            'bam_string': shlex.quote(bam_string),
            'reference': shlex.quote(reference),
            'outfile': shlex.quote(outfile),
            'read1': shlex.quote(read1),
            'bwa_args': quoted_bwa_args
        })
        aligner_command = '; '.join([align_read, samse_command])

    return aligner_command


def _bowtie2_command(path, args, ncpu, reference, sample_name, read1, read2=None):
    """
    Args:
        path (str): path to aligner executable
        args (str): raw arguments to be passed to the aligner
        ncpu: number of alignment threads to launch
        reference: (str): reference filename
        sample_name (str): 
        read1 (str): absolute path to read1 fastq[.gz|.bz2]
        read2 (str): absolute path to read2 fastq[.gz|.bz2]

    Returns:
        string: command to execute bowtie2 aligner
    """
    import os

    reference_basename = os.path.splitext(reference)[0]

    quoted_args = ' '.join(map(shlex.quote, shlex.split(args)))
    quoted_read_args = "-1 {read1} -2 {read2}".format(read1=shlex.quote(read1), read2=shlex.quote(read2)) if read2 else "-U {read1}".format(read1=shlex.quote(read1))

    aligner_command = '{bowtie2} {bowtie2_args} --threads {ncpu} --rg {read_group} --rg-id {read_group_id} -x {bt2_index_prefix} {read_args}'.format(**{
        'bowtie2': path,
        'bowtie2_args': quoted_args,
        'ncpu': shlex.quote(str(ncpu)),
        'read_group': shlex.quote('SM:' + sample_name),
        'read_group_id': shlex.quote(sample_name),
        'bt2_index_prefix': shlex.quote(reference_basename),
        'read_args': quoted_read_args
    })
    return aligner_command


def _novoalign_command(path, args, ncpu, reference, sample_name, read1, read2=None):
    """
    Args:
        path (str): path to aligner executable
        args (str): raw arguments to be passed to the aligner
        ncpu: number of alignment threads to launch
        reference: (str): reference filename
        sample_name (str): 
        read1 (str): absolute path to read1 fastq[.gz|.bz2]
        read2 (str): absolute path to read2 fastq[.gz|.bz2]

    Returns:
        string: command to execute bowtie2 aligner
    """

    import os

    aligner_command = '{novoalign} -d {dbname} -f {read1} {read2} {paired_string} -c {ncpu} -o SAM {bam_string} {novoalign_args}'.format(**{
        'novoalign': path,
        'dbname': shlex.quote(reference + '.idx'),
        'read1': shlex.quote(read1),
        'read2': shlex.quote(read2) if read2 else '',
        'paired_string': '-i PE 500,100' if read2 else '',
        'ncpu': shlex.quote(str(ncpu)),
        'bam_string': shlex.quote('@RG\\tID:{sample_name}\\tSM:{sample_name}'.format(sample_name=sample_name)),
        'novoalign_args': ' '.join(map(shlex.quote, shlex.split(args)))
    })
    return aligner_command
    

def _snap_command(path, args, ncpu, reference, output_folder, sample_name, read1, read2=None):
    """
    Args:
        path (str): path to aligner executable
        args (str): raw arguments to be passed to the aligner
        ncpu: number of alignment threads to launch
        reference: (str): reference filename
        output_folder (str): directory for aligner output
        sample_name (str): 
        read1 (str): absolute path to read1 fastq[.gz|.bz2]
        read2 (str): absolute path to read2 fastq[.gz|.bz2]

    Returns:
        string: command to execute bowtie2 aligner
    """
    import os

    aligner_command = '{snap} {single_or_paired} {ref_dir} {read1} {read2} -t {ncpu} -b {snap_args} -o sam -'.format(**{
        'snap': path,
        'single_or_paired': 'paired' if read2 else 'single',
        'ref_dir': shlex.quote(os.path.join(output_folder, 'reference', 'snap')),
        'read1': shlex.quote(read1),
        'read2': shlex.quote(read2) if read2 else '',
        'ncpu': shlex.quote(str(ncpu)),
        'snap_args': ' '.join(map(shlex.quote, shlex.split(args)))
    })
    return aligner_command


def _gatk_command(path, args, ncpu, mem, output_folder, reference, bam):
    (bam_root, _)= os.path.splitext(bam)
    vcf = os.path.join(output_folder, 'gatk', '{0}-gatk.vcf'.format(bam_root))

    return 'java -Xmx{mem}G -jar {gatk} -T UnifiedGenotyper -dt NONE -glm BOTH -I {bam} -R {reference} -nt {ncpu} -o {vcf} -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE {gatk_args}'.format(**{
        'mem': shlex.quote(str(mem)),
        'gatk': shlex.quote(path),
        'bam': shlex.quote(bam),
        'reference': shlex.quote(reference),
        'ncpu': shlex.quote(str(ncpu)),
        'vcf': shlex.quote(vcf),
        'gatk_args': ' '.join(map(shlex.quote, shlex.split(args)))
    })


def _solsnp_command(path, args, ncpu, mem, output_folder, reference, bam):
    (bam_root, _)= os.path.splitext(bam)
    vcf = os.path.join(output_folder, 'solsnp', '{0}-solsnp.vcf'.format(bam_root))

    return 'java -Xmx{mem}G -jar {solsnp} INPUT={bam} REFERENCE_SEQUENCE={reference} OUTPUT={vcf} SUMMARY=true CALCULATE_ALLELIC_BALANCE=true MINIMUM_COVERAGE=1 PLOIDY=Haploid STRAND_MODE=None OUTPUT_FORMAT=VCF OUTPUT_MODE=AllCallable {solsnp_args}'.format(**{
        'mem': shlex.quote(str(mem)),
        'solsnp': path,
        'bam': shlex.quote(bam),
        'reference': shlex.quote(reference),
        'vcf': shlex.quote(vcf),
        'solsnp_args': ' '.join(map(shlex.quote, shlex.split(args)))
    })
    

def _varscan_command(varscan_path, varscan_args, ncpu, mem, output_folder, reference, bam, sample_name, samtools_path='samtools'):
    import os
    import re

    # TODO: assert bam_filename is not empty
    bam_filename = os.path.basename(bam)
    (bam_root, _)= os.path.splitext(bam_filename)

    vcf = os.path.join(output_folder, '{0}-varscan.vcf'.format(bam_root))
    pileup_file = os.path.join(os.path.dirname(bam), "{0}.mpileup".format(sample_name))
    sample_list = os.path.join(output_folder, "{0}.txt".format(sample_name))

    snpcall_command = '; '.join([
        "echo {sample_name} > {sample_list}".format(**{
            'sample_name': shlex.quote(sample_name),
            'sample_list': shlex.quote(sample_list),
        }),
         "{samtools} mpileup -B -d 10000000 -f {reference} {bam} > {pileup_file}".format(**{
             'samtools': samtools_path,
             'reference': shlex.quote(reference),
             'bam': shlex.quote(bam),
             'pileup_file': shlex.quote(pileup_file),
        }),
       "java -Xmx{mem_gb} -jar {varscan} mpileup2cns {pileup_file} --output-vcf 1 --vcf-sample-list {sample_list} {varscan_args} > {vcf}".format(**{
            'mem_gb': shlex.quote(str(mem) + 'G'),
            'varscan': varscan_path,
            'pileup_file': shlex.quote(pileup_file),
            'sample_list': shlex.quote(sample_list),
            'varscan_args': ' '.join(map(shlex.quote, shlex.split(varscan_args))),
            'vcf': shlex.quote(vcf),
        })
    ])

    return snpcall_command


#def _samtools_command(nickname, bam_file, snpcaller, samtools, job_submitter, aligner_job_id, reference, output_folder):
def _samtools_command(path, args, ncpu, mem, output_folder, reference, bam):
    import os
    (bam_root, _)= os.path.splitext(bam)
    vcf = os.path.join(output_folder, 'samtools', '{0}-samtools.vcf'.format(bam_root))

    return '{samtools} mpileup -uD -d 10000000 -f {reference} {bam} | {samtools} view -ceg {samtools_args} - > {vcf}'.format(**{
        'samtools': path,
        'reference': shlex.quote(reference),
        'bam': shlex.quote(bam),
        'samtools_args': ' '.join(map(shlex.quote, shlex.split(args))),
        'vcf': shlex.quote(vcf)
    })

    
def _run_gatk(nickname, bam_file, snpcaller, job_submitter, aligner_job_id, reference, output_folder):
    import os

    (path, args, job_parms) = snpcaller[1:4]
    snpcaller_name = "gatk"
    ncpus = job_parms['num_cpus']
    memory = job_parms['mem_requested']
    vcf_nickname = "%s-%s" % (nickname, snpcaller_name)
    work_dir = os.path.join(output_folder, snpcaller_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    command = "java -Xmx%sG -jar %s -T UnifiedGenotyper -dt NONE -glm BOTH -I %s -R %s -nt %s -o %s.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE %s" % (
        memory, path, bam_file, reference, ncpus, vcf_nickname, args)
    final_file = os.path.join(work_dir, "%s.vcf" % vcf_nickname)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return vcf_nickname, job_id, final_file


def _run_solsnp(nickname, bam_file, snpcaller, job_submitter, aligner_job_id, reference, output_folder):
    import os

    (path, args, job_parms) = snpcaller[1:4]
    snpcaller_name = "solsnp"
    memory = job_parms['mem_requested']
    vcf_nickname = "%s-%s" % (nickname, snpcaller_name)
    work_dir = os.path.join(output_folder, snpcaller_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.vcf" % vcf_nickname)
    bam_link = os.path.join(work_dir, os.path.splitext(os.path.basename(bam_file))[0])
    os.symlink(bam_file, bam_link)
    command = "java -Xmx%sG -jar %s INPUT=%s REFERENCE_SEQUENCE=%s OUTPUT=%s SUMMARY=true CALCULATE_ALLELIC_BALANCE=true MINIMUM_COVERAGE=1 PLOIDY=Haploid STRAND_MODE=None OUTPUT_FORMAT=VCF OUTPUT_MODE=AllCallable %s" % (
        memory, path, bam_link, reference, final_file, args)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return vcf_nickname, job_id, final_file


def _run_varscan(nickname, bam_file, snpcaller, samtools, job_submitter, aligner_job_id, reference, output_folder):
    import os
    import re

    (path, args, job_parms) = snpcaller[1:4]
    sampath = samtools[1]
    snpcaller_name = "varscan"
    memory = job_parms['mem_requested']
    vcf_nickname = "%s-%s" % (nickname, snpcaller_name)
    test = re.match('^(.*)-[a-z]*$', nickname, re.IGNORECASE)
    read_nickname = test.group(1) if test else nickname
    work_dir = os.path.join(output_folder, snpcaller_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    sample_list = os.path.join(work_dir, "%s.txt" % nickname)
    pileup_file = os.path.join(os.path.dirname(bam_file), "%s.mpileup" % nickname)
    final_file = os.path.join(work_dir, "%s.vcf" % vcf_nickname)
    command_parts = ["echo %s > %s" % (read_nickname, sample_list),
                     "%s mpileup -B -d 10000000 -f %s %s > %s" % (sampath, reference, bam_file, pileup_file),
                     "java -Xmx%sG -jar %s mpileup2cns %s --output-vcf 1 --vcf-sample-list %s > %s %s" % (
                         memory, path, pileup_file, sample_list, final_file, args)]
    command = "\n".join(command_parts)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return vcf_nickname, job_id, final_file


def _run_samtools(nickname, bam_file, snpcaller, samtools, job_submitter, aligner_job_id, reference, output_folder):
    import os

    (path, args, job_parms) = snpcaller[1:4]
    sampath = samtools[1]
    snpcaller_name = "samtools"
    vcf_nickname = "%s-%s" % (nickname, snpcaller_name)
    work_dir = os.path.join(output_folder, snpcaller_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.vcf" % vcf_nickname)
    command_parts = ["%s mpileup -uD -d 10000000 -f %s %s" % (sampath, reference, bam_file),
                     "%s view -ceg %s - > %s" % (path, args, final_file)]
    command = " | ".join(command_parts)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return vcf_nickname, job_id, final_file


def _find_dups(configuration, index_job_id, reference):
    import os

    (name, path, _, job_parms) = configuration["dup_finder"]
    command = "find_duplicates --nucmerpath %s --reference %s" % (path, reference)
    work_dir = os.path.dirname(reference)
    if not os.path.exists(work_dir):
        # NOTE: This directory is implicitly created by _index_reference.
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "duplicates.txt")
    job_parms['name'] = "nasp_%s" % name
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
    return job_id, final_file


def _convert_external_genome(assembly, configuration, index_job_id, reference):
    import os

    (tool, path, args, job_parms) = configuration["assembly_importer"]
    (nucmer_path, nucmer_args) = configuration["dup_finder"][1:3]
    (name, fasta) = assembly
    work_dir = os.path.join(configuration["output_folder"], "external")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    new_fasta = os.path.join(work_dir, os.path.basename(fasta))
    command_parts = ["format_fasta --inputfasta %s --outputfasta %s" % (fasta, new_fasta),
                     "convert_external_genome --nucmerpath %s --nucmerargs \'%s\' --deltafilterpath %s --deltafilterargs \'%s\' --reference %s --external %s --name %s" % (
                         nucmer_path, nucmer_args, path, args, reference, new_fasta, name)]
    command = "\n".join(command_parts)
    final_file = os.path.join(work_dir, "%s.frankenfasta" % name)
    job_parms['name'] = "nasp_%s_%s" % (tool, name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
    return job_id, final_file


# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
def _trimmomatic_command(path, args, ncpu, output_folder, sample_name, read1, read2=None):

    parameters = {
        'trimmomatic': path,
        'ncpu': shlex.quote(str(ncpu)),
        'input_read1': shlex.quote(read1),
        'input_read2': shlex.quote(read2),
        'output_pe_read1_paired': shlex.quote(os.path.join(output_folder, 'trimmomatic', sample_name + '_R1_paired.fastq.gz')),
        'output_pe_read1_unpaired': shlex.quote(os.path.join(output_folder, 'trimmomatic', sample_name+ '_R1_unpaired.fastq.gz')),
        'output_pe_read2_paired': shlex.quote(os.path.join(output_folder, 'trimmomatic', sample_name+ '_R2_paired.fastq.gz')),
        'output_pe_read2_unpaired': shlex.quote(os.path.join(output_folder, 'trimmomatic', sample_name+ '_R2_unpaired.fastq.gz')),
        'output_se_trimmed': shlex.quote(os.path.join(output_folder, 'trimmomatic', sample_name+ '_trimmed.fastq.gz')),
        'trimmomatic_args': ' '.join(map(shlex.quote, shlex.split(args)))
    }

    if read2:
        # Paired End
        return 'java -jar {trimmomatic} PE -threads {ncpu} {input_read1} {input_read2} {output_pe_read1_paired} {output_pe_read2_unpaired} {output_pe_read2_paired} {output_pe_read2_unpaired} {trimmomatic_args}'.format(**parameters)

    # Single End
    return 'java -jar {trimmomatic} SE -threads {ncpu} {input_read1} {output_se_trimmed} {trimmomatic_args}'.format(**parameters)


def _trim_adapters(read_tuple, configuration):
    import os

    (_, path, args, job_parms) = configuration["read_trimmer"]
    (name, read1) = read_tuple[0:2]
    read2 = read_tuple[2] if len(read_tuple) >= 3 else None
    trim_dir = os.path.join(configuration["output_folder"], 'trimmed')
    if not os.path.exists(trim_dir):
        os.makedirs(trim_dir)
    #job_params = {'queue':'', 'mem_requested':6, 'num_cpus':4, 'walltime':8, 'args':''}
    job_parms['name'] = "nasp_trim_%s" % name
    job_parms['work_dir'] = trim_dir
    if read2:
        out_reads1 = [name+"_R1_trimmed.fastq", name+"_R1_unpaired.fastq"]
        out_reads2 = [name+"_R2_trimmed.fastq", name+"_R2_unpaired.fastq"]
        out_reads = [os.path.join(trim_dir, out_reads1[0]), os.path.join(trim_dir, out_reads2[0])]
        command = "java -jar %s PE -threads %s %s %s %s %s %s %s %s" % (path, job_parms['num_cpus'], read1, read2, out_reads1[0], out_reads1[1], out_reads2[0], out_reads2[1], args)
    else:
        out_reads = [os.path.join(trim_dir, name+"_trimmed.fastq")]
        command = "java -jar %s SE -threads %s %s %s %s" % (path, job_parms['num_cpus'], read1, out_reads[0], args)
    jobid = _submit_job(configuration["job_submitter"], command, job_parms)
    return (tuple([name] + out_reads), jobid)


def _align_reads(read_tuple, configuration, index_job_id, reference):
    import re
    import os

    output_folder = configuration['output_folder']
    samtools_path = configuration['samtools'][1]
    job_submitter = configuration['job_submitter']
    aligner_output = []

    for aligner in map(App._make, configuration['aligners']):
        aligner_name = ''
        align_command = ''

        if re.search('bwa', aligner.name, re.IGNORECASE):
            if re.search('mem', aligner.name, re.IGNORECASE):
                aligner_name = 'bwamem'
                align_command = _bwamem_command(aligner.path, aligner.args, aligner.job_params['num_cpus'], reference, *read_tuple)
            else:
                aligner_name = 'bwa'
                align_command = _bwa_command(aligner.path, aligner.args, aligner.job_params['num_cpus'], reference, output_folder, *read_tuple)
        elif re.search('b(ow)?t(ie)?2', aligner.name, re.IGNORECASE):
            aligner_name = 'bowtie2'
            align_command = _bowtie2_command(aligner.path, aligner.args, aligner.job_params['num_cpus'], reference, *read_tuple)
        elif re.search('novo', aligner.name, re.IGNORECASE):
            aligner_name = 'novo'
            align_command = _novoalign_command(aligner.path, aligner.args, aligner.job_params['num_cpus'], reference, *read_tuple)
        elif re.search('snap', aligner.name, re.IGNORECASE):
            aligner_name = 'snap'
            align_command = _snap_command(aligner.path, aligner.args, aligner.job_params['num_cpus'], reference, output_folder, *read_tuple)
        else:
            print("Unknown aligner \'{0}\' found, don't know what to do. Skipping...".format(aligner.name))
            continue

        output_bam = os.path.join(output_folder, aligner_name, '{sample}-{aligner}.bam'.format(sample=read_tuple[0], aligner=aligner_name))
        
        command = align_command + ' | ' + _samtools_view_sort_index_pipe_command(samtools_path, output_bam)

        work_dir = os.path.join(output_folder, aligner_name)
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)

        outfile = os.path.join(work_dir, "{bam_prefix}.bam".format(bam_prefix=bam_prefix))
        aligner.job_params['name'] = "nasp_{aligner}_{sample}".format(aligner=aligner_name, sample=read_tuple[0])
        aligner.job_params['work_dir'] = work_dir

        job_id = _submit_job(job_submitter, align_command, aligner.job_params, (index_job_id,))
        if job_id:
            aligner_output.append((bam_prefix, job_id, outfile, aligner_name))

    return aligner_output


def _call_snps(aligner_output, configuration, reference):
    import re

    # TODO: skip task if the file exists and a --force flag was not used
    # TODO: implement a GNU make style --force flag
    snpcaller_output = []
    for (nickname, aligner_job_id, bam_file, aligner_name) in aligner_output:
        if aligner_job_id:
            for snpcaller in configuration["snpcallers"]:
                name = snpcaller[0]
                if re.search('gatk', name, re.IGNORECASE):
                    (vcf_nickname, job_id, final_file) = _run_gatk(nickname, bam_file, snpcaller,
                                                                   configuration["job_submitter"], aligner_job_id,
                                                                   reference, configuration["output_folder"])
                    if job_id:
                        snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
                elif re.search('solsnp', name, re.IGNORECASE):
                    (vcf_nickname, job_id, final_file) = _run_solsnp(nickname, bam_file, snpcaller,
                                                                     configuration["job_submitter"], aligner_job_id,
                                                                     reference, configuration["output_folder"])
                    if job_id:
                        snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
                elif re.search('varscan', name, re.IGNORECASE):
                    (vcf_nickname, job_id, final_file) = _run_varscan(nickname, bam_file, snpcaller,
                                                                      configuration["samtools"],
                                                                      configuration["job_submitter"], aligner_job_id,
                                                                      reference, configuration["output_folder"])
                    if job_id:
                        snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
                elif re.search('samtools', name, re.IGNORECASE):
                    (vcf_nickname, job_id, final_file) = _run_samtools(nickname, bam_file, snpcaller,
                                                                       configuration["samtools"],
                                                                       configuration["job_submitter"], aligner_job_id,
                                                                       reference, configuration["output_folder"])
                    if job_id:
                        snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
                else:
                    print("Unknown SNP caller \'%s\' found, don't know what to do. Skipping..." % name)
    return snpcaller_output


def _index_bams(configuration, index_job_id):
    import os

    alignments = configuration["alignments"]
    output_folder = configuration["output_folder"]
    job_parms = configuration["bam_index"][3]
    sampath = configuration["samtools"][1]
    picard_path = configuration["picard"][1] or ""
    bam_folder = os.path.join(output_folder, "bams")
    if not os.path.exists(bam_folder):
        os.makedirs(bam_folder)
    bam_files = []
    command_parts = []
    for (name, bam) in alignments:
        new_file = os.path.join(bam_folder, "%s.bam" % name)
        # GATK requires bams have a ReadGroup header. If one if not present,
        # assign a default value.
        command_parts.append((
            "if ! {samtools} view -H {in_bam} | grep \"^@RG\"; then\n"
            " java -jar {picard} AddOrReplaceReadGroups"
            " INPUT={in_bam}"
            " OUTPUT={out_bam}"
            " SORT_ORDER=coordinate"
            " ID=1"
            " SM={sample_name}"
            " LB={sample_name}"
            " PL=illumina"
            " PU=1\n"
            "else\n"
            " ln -s -f {in_bam} {out_bam}\n"
            "fi"
        ).format(samtools=sampath, picard=picard_path, out_bam=new_file, in_bam=bam, sample_name=os.path.splitext(bam)[0]))
        command_parts.append("%s index %s" % (sampath, new_file))
        command = "\n".join(command_parts)
        job_parms['work_dir'] = bam_folder
        job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
        bam_files.append((name, new_file, job_id))
    return bam_files


def _create_matrices(configuration, reference, dups_file, vcf_files, franken_fastas, job_ids):
    import nasp.matrix_DTO as matrix_DTO
    import os

    output_dir = configuration['output_folder']
    path = configuration["matrix_generator"][1]
    job_parms = configuration["matrix_generator"][3]
    matrix_parms = {'reference-fasta': reference, 'reference-dups': dups_file}
    if "coverage_filter" in configuration:
        matrix_parms['minimum-coverage'] = configuration['coverage_filter']
    if "proportion_filter" in configuration:
        matrix_parms['minimum-proportion'] = configuration['proportion_filter']
    matrix_parms['matrix-folder'] = os.path.join(output_dir, 'matrices')
    if not os.path.exists(matrix_parms['matrix-folder']):
        os.makedirs(matrix_parms['matrix-folder'])
    matrix_parms['stats-folder'] = os.path.join(output_dir, 'statistics')
    if not os.path.exists(matrix_parms['stats-folder']):
        os.makedirs(matrix_parms['stats-folder'])
    if 'filter_matrix_format' in configuration:
        matrix_parms['filter-matrix-format'] = configuration['filter_matrix_format']
    dto_file = os.path.join(output_dir, "matrix_dto.xml")
    matrix_DTO.write_dto(matrix_parms, franken_fastas, vcf_files, dto_file)
    jobs_to_wait_for = (":".join(job_ids), 'afterany') if job_ids else None
    command = "%s matrix --dto-file %s --num-threads %s" % (path, dto_file, job_parms['num_cpus'])
    job_parms['work_dir'] = output_dir
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, jobs_to_wait_for, notify=True)
    return job_id


def _export_matrices(configuration, matrix_job_id):
    import os

    gonasp_path = configuration["matrix_generator"][1]
    matrix_folder = os.path.join(configuration['output_folder'], 'matrices')
    job_parms = configuration["matrix_generator"][3]
    job_parms['name'] = "nasp_export"
    job_parms['work_dir'] = matrix_folder
    commands = []

    # The command will be of the following form:
    #   gonasp export --type fasta bestsnp.tsv > bestsnp.fasta &
    #   gonasp export --type vcf bestsnp.tsv > bestsnp.vcf &
    #   wait
    # The '&' will launch each export commands as background tasks then the 'wait' will wait for them to finish
    for file_type in ['vcf', 'fasta']:
        for exported_matrix in ['bestsnp', 'missingdata']:
            commands.append("{0} export --type {1} {2}.tsv > {2}.{1}".format(gonasp_path, file_type, exported_matrix))

    command = '; '.join(commands)

    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (matrix_job_id,), notify=False)


def begin(configuration):
    (index_job_id, reference) = _index_reference(configuration)
    if not index_job_id:
        print("Failed to submit the index job, there is no point in continuing. Please try again.")
        raise SystemExit()
    dups_file = None
    job_ids = []
    vcf_files = []
    franken_fastas = []
    if configuration["find_dups"] != "False":
        (job_id, dups_file) = _find_dups(configuration, index_job_id, reference)
        if job_id:
            job_ids.append(job_id)
    for assembly in configuration["assemblies"]:
        (job_id, final_file) = _convert_external_genome(assembly, configuration, index_job_id, reference)
        if job_id:
            job_ids.append(job_id)
            franken_fastas.append((assembly[0], "nucmer", final_file))
    if configuration["alignments"]:
        pre_aligned = []
        (bam_files) = _index_bams(configuration, index_job_id)
        for (name, bam, bamindex_job_id) in bam_files:
            pre_aligned.append((name, bamindex_job_id, bam, "pre-aligned"))
        snpcaller_output = _call_snps(pre_aligned, configuration, reference)
        for (vcf_nickname, job_id, final_file, aligner, snpcaller) in snpcaller_output:
            if job_id:
                job_ids.append(job_id)
                vcf_files.append((vcf_nickname, aligner, snpcaller, final_file))
    for read_tuple in configuration["reads"]:
        dependent_job_id = None
        if "trim_reads" in configuration and configuration["trim_reads"] == "True":
            (read_tuple, dependent_job_id) = _trim_adapters(read_tuple, configuration)
        dependencies = index_job_id
        if dependent_job_id:
            dependencies += ":"+dependent_job_id
        aligner_output = _align_reads(read_tuple, configuration, dependencies, reference)
        snpcaller_output = _call_snps(aligner_output, configuration, reference)
        for (vcf_nickname, job_id, final_file, aligner, snpcaller) in snpcaller_output:
            if job_id:
                job_ids.append(job_id)
                vcf_files.append((vcf_nickname, aligner, snpcaller, final_file))
    for (name, vcf) in configuration["vcfs"]:
        vcf_files.append((name, "pre-aligned", "pre-called", vcf))

    matrix_job_id = _create_matrices(configuration, reference, dups_file, vcf_files, franken_fastas, job_ids)
    _export_matrices(configuration, matrix_job_id)
    _release_hold(configuration["job_submitter"], index_job_id)


def main():
    import nasp.configuration_parser as configuration_parser

    commandline_args = _parse_args()
    configuration = configuration_parser.parse_config(commandline_args.config)
    begin(configuration)


if __name__ == "__main__":
    main()
