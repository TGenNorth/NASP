#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.5"
__email__ = "dlemmer@tgen.org"

'''
Created on Mar 4, 2014

@author: dlemmer
'''

import logging
import configuration_parser

def _parse_args():
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--config", required=True, help="Path to the configuration xml file." )
    return parser.parse_args()

def _submit_job( job_submitter, command, job_parms, waitfor_id=None ):
    import subprocess
    import re
    output = jobid = None
    if job_submitter == "qsub":
        waitfor = ""
        if waitfor_id:
            dependency_string = waitfor_id[1] if len(waitfor_id) > 1 else 'afterok'
            waitfor = "-W depend=%s:%s" % (dependency_string, waitfor_id[0])
        submit_command = "qsub -d \'%s\' -w \'%s\' -l ncpus=%s,mem=%sgb,walltime=%s:00:00 -m a -N \'%s\' %s" % (job_parms["work_dir"], job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'], job_parms['walltime'], job_parms['name'], waitfor)
        print("command = %s" % command)
        output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
        jobid = re.search('^(\d+)\..*$', output).group(1)
    else:
        pass
    print("jobid = %s" % jobid)
    return(jobid)

def _index_reference( configuration ):
    import os
    import re
    output_folder = configuration["output_folder"]
    ref_path = configuration["reference"][1]
    ref_folder = os.path.join(output_folder, "reference")
    if not os.path.exists(ref_folder):
        os.makedirs(ref_folder)
    #Copy the reference as $output_folder/reference/reference.fasta, verifying its format first. Replace it if it already exists.
    reference = os.path.join(ref_folder, "reference.fasta")
    if os.path.exists(reference):
        os.remove(reference)
    index_commands = ["format_fasta.py --inputfasta %s --outputfasta %s" % (ref_path, reference)] 
    
    #Gather all of the index commands that need to be run
    bwa_done = False
    for aligner in configuration["aligners"]:
        (name, path) = aligner[0:2]
        if re.search('bwa', name, re.IGNORECASE):
            if not bwa_done:
                index_commands.append("%s index %s" % (path, reference))
                bwa_done = True
        elif re.search('novo', name, re.IGNORECASE):
            novopath = os.path.split(path)[0]
            novoindex_path = os.path.join(novopath, "novoindex")
            index_commands.append("%s %s.idx %s" % (novoindex_path, reference, reference))
        elif re.search('snap', name, re.IGNORECASE):
            index_commands.append("%s index %s %s" % (path, reference, os.path.join(ref_folder, "snap")))
        else:
            print("Unknown aligner \'%s\' found, don't know how to index the reference for it. Skipping..." % name)
    
    #if we are using GATK, we also need to create a Sequence Dictionary and samtools index of the reference        
    if next((v for i,v in enumerate(configuration["snpcallers"]) if re.search('gatk', v[0], re.IGNORECASE)), None):
        picard_path = configuration["picard"][1] or ""
        picard_memory = 2
        if configuration["picard"][3]:
            picard_memory = configuration["picard"][3]['mem_requested'] or 2
        dict_generator = os.path.join(picard_path, "CreateSequenceDictionary.jar")
        samtools_path = configuration["samtools"][1] or "samtools"
        out_file = os.path.join(ref_folder, "reference.dict")
        index_commands.append("java -Xmx%sG -jar %s R=%s O=%s" % (picard_memory, dict_generator, reference, out_file))
        index_commands.append("%s faidx %s" % (samtools_path, reference))
    
    command = "\n".join(index_commands)
    job_parms = {'name':'nasp_index', 'num_cpus':'1', 'mem_requested':'2', 'walltime':'4', 'work_dir':ref_folder}
    job_id = _submit_job(configuration["job_submitter"], command, job_parms)
    return (job_id, reference)

def _run_bwa(read_tuple, aligner, samtools, job_submitter, index_job_id, reference, output_folder):
    import re
    import os
    (name, read1) = read_tuple[0:2]
    read2 = read_tuple[2] if len(read_tuple) >= 3 else ""
    bam_string = "\'@RG\\tID:%s\\tSM:%s\'" % (name, name)
    sampath = samtools[1]
    (aligner_name, path, args, job_parms) = aligner
    ncpus = job_parms['num_cpus']
    aligner_command = ""
    if re.search('mem', aligner_name, re.IGNORECASE):
        aligner_name = "bwamem"
        aligner_command = "%s mem -R %s %s -t %s %s %s %s" % (path, bam_string, args, ncpus, reference, read1, read2)
    else:
        aligner_name = "bwa"
        old_format_string = "-I" if re.search('(?:.*\/)?[^\/]+?_[12]_sequence\.txt(?:\.gz)?$', read1, re.IGNORECASE) else ""
        work_dir = os.path.join(output_folder, aligner_name)
        command_parts = []
        if read2:
            output_file_1 = os.path.join(work_dir, "%s-R1.sai" % name)
            output_file_2 = os.path.join(work_dir, "%s-R2.sai" % name)
            command_parts.append("%s aln %s %s %s -t %s -f %s %s" % (path, old_format_string, reference, read1, ncpus, output_file_1, args))
            command_parts.append("%s aln %s %s %s -t %s -f %s %s" % (path, old_format_string, reference, read2, ncpus, output_file_2, args))
            command_parts.append("%s sampe -r %s %s %s %s %s %s %s" % (path, bam_string, reference, output_file_1, output_file_2, read1, read2, args))
        else:
            output_file = os.path.join(work_dir, "%s.sai" % name)
            command_parts.append("%s aln %s %s %s -t %s -f %s %s" % (path, old_format_string, reference, read1, ncpus, output_file, args))
            command_parts.append("%s samse -r %s %s %s %s %s %s %s" % (path, bam_string, reference, output_file, read1, args))
        aligner_command = "\n".join(command_parts)
    bam_nickname = "%s-%s" % (name, aligner_name)    
    samview_command = "%s view -S -b -h -" % (sampath)
    samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s | %s | %s \n %s" % (aligner_command, samview_command, samsort_command, samindex_command)
    work_dir = os.path.join(output_folder, aligner_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)    
    job_parms['name'] = "nasp_%s_%s" % (aligner_name, name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (index_job_id,))
    return (bam_nickname, job_id, final_file)

def _run_novoalign(read_tuple, aligner, samtools, job_submitter, index_job_id, reference, output_folder):
    import os
    (name, read1) = read_tuple[0:2]
    read2 = read_tuple[2] if len(read_tuple) >= 3 else ""
    paired_string = "-i PE 500,100" if read2 else ""
    bam_string = "\'@RG\\tID:%s\\tSM:%s\'" % (name, name)
    sampath = samtools[1]
    (path, args, job_parms) = aligner[1:4]
    aligner_name = "novo"
    ncpus = job_parms['num_cpus']
    aligner_command = "%s -f %s %s %s -c %s -o SAM %s -d %s.idx %s" % (path, read1, read2, paired_string, ncpus, bam_string, reference, args)
    bam_nickname = "%s-%s" % (name, aligner_name)    
    samview_command = "%s view -S -b -h -" % (sampath)
    samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s | %s | %s \n %s" % (aligner_command, samview_command, samsort_command, samindex_command)
    work_dir = os.path.join(output_folder, aligner_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)    
    job_parms['name'] = "nasp_%s_%s" % (aligner_name, name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (index_job_id,))
    return (bam_nickname, job_id, final_file)

def _run_snap(read_tuple, aligner, samtools, job_submitter, index_job_id, reference, output_folder):
    import re
    import os
    (name, read1) = read_tuple[0:2]
    read2 = read_tuple[2] if len(read_tuple) >= 3 else ""
    paired_string = "paired" if read2 else "single"
    sampath = samtools[1]
    (path, args, job_parms) = aligner[1:4]
    aligner_name = "snap"
    unzip = []
    remove_temp = []
    reads = []
    for read in (read1, read2):
        match = re.match('/^(.+)\.gz$/', read, re.IGNORECASE)
        if match:
            decompressed = os.path.join(os.path.join(output_folder, aligner_name), match.group(0))
            unzip.append("zcat %s > %s" % (read, decompressed))
            reads.append(decompressed)
            remove_temp.append("rm %s" % decompressed)
        else:
            reads.append(read)
    read_string = " ".join(reads)
    unzip_command = "\n".join(unzip)
    remove_temp_command = "\n".join(remove_temp)
    ref_dir = os.path.join(os.path.join(output_folder, "reference"), aligner_name)
    aligner_command = "%s %s %s %s -o %s-%s.sam %s" % (path, paired_string, ref_dir, read_string, name, aligner_name, args)
    bam_nickname = "%s-%s" % (name, aligner_name)    
    samview_command = "%s view -S -b -h %s.sam" % (sampath, bam_nickname)
    samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s \n %s \n %s | %s \n %s \n %s" % (unzip_command, aligner_command, samview_command, samsort_command, samindex_command, remove_temp_command)
    work_dir = os.path.join(output_folder, aligner_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)    
    job_parms['name'] = "nasp_%s_%s" % (aligner_name, name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (index_job_id,))
    return (bam_nickname, job_id, final_file)

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
    command = "java -Xmx%sG -jar %s -T UnifiedGenotyper -dt NONE -glm BOTH -I %s -R %s -nt %s -ploidy 1 -o %s.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE %s" % (memory, path, bam_file, reference, ncpus, vcf_nickname, args)
    final_file = os.path.join(work_dir, "%s.vcf" % vcf_nickname)    
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return (vcf_nickname, job_id, final_file)    

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
    command = "java -Xmx%sG -jar %s INPUT=%s REFERENCE_SEQUENCE=%s OUTPUT=%s SUMMARY=true CALCULATE_ALLELIC_BALANCE=true MINIMUM_COVERAGE=1 PLOIDY=Haploid STRAND_MODE=None OUTPUT_FORMAT=VCF OUTPUT_MODE=AllCallable %s" % (memory, path, bam_file, reference, final_file, args)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return (vcf_nickname, job_id, final_file)    

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
    command_parts = ["echo %s > %s" % (read_nickname, sample_list)]
    command_parts.append("%s mpileup -B -d 10000000 -f %s %s > %s" % (sampath, reference, bam_file, pileup_file))
    command_parts.append("java -Xmx%sG -jar %s mpileup2cns %s --output-vcf 1 --vcf-sample-list %s > %s %s" % (memory, path, pileup_file, sample_list, final_file, args))
    command = "\n".join(command_parts)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return (vcf_nickname, job_id, final_file)    

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
    command_parts = ["%s mpileup -uD -d 10000000 -f %s %s" % (sampath, reference, bam_file)]
    command_parts.append("%s view -ceg %s - > %s" % (path, args, final_file))
    command = " | ".join(command_parts)
    job_parms['name'] = "nasp_%s_%s" % (snpcaller_name, nickname)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(job_submitter, command, job_parms, (aligner_job_id,))
    return (vcf_nickname, job_id, final_file)    

def _find_dups( configuration, index_job_id, reference ):
    import os
    (name, path, args, job_parms) = configuration["dup_finder"]
    command = "find_duplicates.py --nucmerpath %s --reference %s" % (path, reference)
    work_dir = os.path.dirname(reference)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "duplicates.txt")
    job_parms['name'] = "nasp_%s" % (name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
    return (job_id, final_file)

def _convert_external_genome( assembly, configuration, index_job_id, reference ):
    import os
    (tool, path, args, job_parms) = configuration["assembly_importer"]
    (nucmer_path, nucmer_args) = configuration["dup_finder"][1:3]
    (name, fasta) = assembly
    extraargs = "\'%s\'" % args
    work_dir = os.path.join(configuration["output_folder"], "external")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    new_fasta = os.path.join(work_dir, os.path.basename(fasta))
    command_parts = ["format_fasta.py --inputfasta %s --outputfasta %s" % (fasta, new_fasta)]
    command_parts.append("convert_external_genome.py --nucmerpath %s --nucmerargs %s --deltafilterpath %s --reference %s --external %s --name %s" % (nucmer_path, extraargs, path, reference, fasta, name))
    command = "\n".join(command_parts)
    final_file = os.path.join(work_dir, "%s.frankenfasta" % name)
    job_parms['name'] = "nasp_%s_%s" % (tool, name)
    job_parms['work_dir'] = work_dir
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
    return (job_id, final_file)

def _align_reads( read_tuple, configuration, index_job_id, reference ):
    import re
    aligner_output = []
    for aligner in configuration["aligners"]:
        name = aligner[0]
        if re.search('bwa', name, re.IGNORECASE):
            (bam_nickname, job_id, final_file) = _run_bwa(read_tuple, aligner, configuration["samtools"], configuration["job_submitter"], index_job_id, reference, configuration["output_folder"])
            aligner_output.append((bam_nickname, job_id, final_file, name))
        elif re.search('novo', name, re.IGNORECASE):
            (bam_nickname, job_id, final_file) = _run_novoalign(read_tuple, aligner, configuration["samtools"], configuration["job_submitter"], index_job_id, reference, configuration["output_folder"])
            aligner_output.append((bam_nickname, job_id, final_file, name))
        elif re.search('snap', name, re.IGNORECASE):
            (bam_nickname, job_id, final_file) = _run_snap(read_tuple, aligner, configuration["samtools"], configuration["job_submitter"], index_job_id, reference, configuration["output_folder"])
            aligner_output.append((bam_nickname, job_id, final_file, name))
        else:
            print("Unknown aligner \'%s\' found, don't know what to do. Skipping..." % name)
    return aligner_output

def _call_snps( aligner_output, configuration, reference ):
    import re
    snpcaller_output = []
    print(aligner_output)
    print(len(aligner_output))
    for (nickname, aligner_job_id, bam_file, aligner_name) in aligner_output:
        for snpcaller in configuration["snpcallers"]:
            name = snpcaller[0]
            if re.search('gatk', name, re.IGNORECASE):
                (vcf_nickname, job_id, final_file) = _run_gatk(nickname, bam_file, snpcaller, configuration["job_submitter"], aligner_job_id, reference, configuration["output_folder"])
                snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
            elif re.search('solsnp', name, re.IGNORECASE):
                (vcf_nickname, job_id, final_file) = _run_solsnp(nickname, bam_file, snpcaller, configuration["job_submitter"], aligner_job_id, reference, configuration["output_folder"])
                snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
            elif re.search('varscan', name, re.IGNORECASE):
                (vcf_nickname, job_id, final_file) = _run_varscan(nickname, bam_file, snpcaller, configuration["samtools"], configuration["job_submitter"], aligner_job_id, reference, configuration["output_folder"])
                snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
            elif re.search('samtools', name, re.IGNORECASE):
                (vcf_nickname, job_id, final_file) = _run_samtools(nickname, bam_file, snpcaller, configuration["samtools"], configuration["job_submitter"], aligner_job_id, reference, configuration["output_folder"])
                snpcaller_output.append((vcf_nickname, job_id, final_file, aligner_name, name))
            else:
                print("Unknown SNP caller \'%s\' found, don't know what to do. Skipping..." % name)
    return snpcaller_output

def _index_bams( configuration, index_job_id ):
    import os
    alignments = configuration["alignments"]
    output_folder = configuration["output_folder"]
    sampath = configuration["samtools"][1]
    bam_folder = os.path.join(output_folder, "bams")
    if not os.path.exists(bam_folder):
        os.makedirs(bam_folder)
    bam_files = []
    command_parts = []
    for (name, bam) in alignments:
        new_file = os.path.join(bam_folder, "%s.bam" % name)
        bam_files.append((name, new_file))
        command_parts.append("ln -s -f %s %s" % (bam, new_file))
        command_parts.append("%s index %s" % (sampath, new_file))
    command = "\n".join(command_parts)
    job_parms = {'name':'nasp_bamindex', 'num_cpus':'1', 'mem_requested':'2', 'walltime':'4', 'work_dir':bam_folder}
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (index_job_id,))
    return (bam_files, job_id)

def _create_matrices( configuration, reference, dups_file, vcf_files, franken_fastas, job_ids ):
    import matrix_DTO
    import os
    output_dir = configuration['output_folder']
    matrix_parms = {'reference-fasta':reference, 'reference-dups':dups_file}
    matrix_parms['minimum-coverage'] = configuration['coverage_filter']
    matrix_parms['minimum-proportion'] = configuration['proportion_filter']
    matrix_parms['master-matrix'] = os.path.join(output_dir, 'master_matrix.tsv')
    matrix_parms['filter-matrix'] = os.path.join(output_dir, 'filter_matrix.tsv')
    matrix_parms['general-stats'] = os.path.join(output_dir, 'general_stats.tsv')
    matrix_parms['contig-stats'] = os.path.join(output_dir, 'contig_stats.tsv')
    dto_file = os.path.join(output_dir, "matrix_dto.xml")
    matrix_DTO.write_dto(matrix_parms, franken_fastas, vcf_files, dto_file)
    jobs_to_wait_for = ":".join(job_ids)
    command = "vcf_to_matrix.py --mode xml --dto-file %s --num-threads %s" % (dto_file, '10')
    job_parms = {'name':'nasp_matrix', 'num_cpus':'10', 'mem_requested':'4', 'walltime':'96', 'work_dir':output_dir}
    job_id = _submit_job(configuration["job_submitter"], command, job_parms, (jobs_to_wait_for, 'afterany'))    
    return job_id

def main():
    commandline_args = _parse_args()
    configuration = configuration_parser.parse_config( commandline_args.config )
    (index_job_id, reference) = _index_reference( configuration )
    dups_file = None
    job_ids = []
    vcf_files = []
    franken_fastas = []
    if configuration["find_dups"]:
        (job_id, dups_file) = _find_dups( configuration, index_job_id, reference )
    for assembly in configuration["assemblies"]:
        (job_id, final_file) = _convert_external_genome( assembly, configuration, index_job_id, reference )
        job_ids.append(job_id)
        franken_fastas.append((assembly[0], "nucmer", final_file))
    if configuration["alignments"]:
        pre_aligned = []
        (bam_files, bamindex_job_id) = _index_bams(configuration, index_job_id)
        for (name, bam) in bam_files:
            pre_aligned.append((name, bamindex_job_id, bam, "pre-aligned"))
        snpcaller_output = _call_snps( pre_aligned, configuration, reference )
        for (vcf_nickname, job_id, final_file, aligner, snpcaller) in snpcaller_output:
            job_ids.append(job_id)
            vcf_files.append((vcf_nickname, aligner, snpcaller, final_file))
    for read_tuple in configuration["reads"]:
        aligner_output = _align_reads( read_tuple, configuration, index_job_id, reference )
        snpcaller_output = _call_snps( aligner_output, configuration, reference )
        for (vcf_nickname, job_id, final_file, aligner, snpcaller) in snpcaller_output:
            job_ids.append(job_id)
            vcf_files.append((vcf_nickname, aligner, snpcaller, final_file))
    for (name, vcf) in configuration["vcfs"]:
        vcf_files.append((name, "pre-aligned", "pre-called", vcf))
        
    _create_matrices(configuration, reference, dups_file, vcf_files, franken_fastas, job_ids)
    

if __name__ == "__main__": main()
