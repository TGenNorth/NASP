#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.7"
__email__ = "dlemmer@tgen.org"

'''
Created on Jun 11, 2014

@author: dlemmer
'''

nasp_version = __version__
import logging


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(prog="nasp",
                                     description="This is the experimental \"Northern Arizona SNP Pipeline\", version %s" % nasp_version)
    parser.add_argument("reference_fasta", nargs="?", default="", help="Path to the reference fasta.")
    parser.add_argument("output_folder", nargs="?", default="", help="Folder to store the output files.")
    parser.add_argument("--config", help="Path to the configuration xml file.")
    return parser.parse_args()


def _expand_path(path):
    import os
    import re

    user_match = re.match('^(~)(.*)$', path)
    if user_match:
        path = os.path.expanduser(path)
    return os.path.abspath(path)


def _create_file_tuple(path):
    import os

    name = os.path.splitext(os.path.basename(path))[0]
    return name, path


def _find_files(path, extension):
    import os
    import re

    file_list = []
    for file in os.listdir(path):
        is_type = re.search('(.*)(\.%s)$' % extension, file, re.IGNORECASE)
        if is_type:
            sample_name = is_type.group(1)
            file_list.append((sample_name, os.path.join(path, file)))
            logging.info((sample_name, os.path.join(path, file)))
    return file_list


def _find_executable(application):  # This method is not OS-independent. Should work on a better way
    import re
    import subprocess

    executable = ""
    try:
        executable = subprocess.check_output(["which", application])
    except subprocess.CalledProcessError:
        try:
            executable = subprocess.check_output("find ~ -name %s" % application, shell=True)
        except subprocess.CalledProcessError as cpe:
            executable = cpe.output
    match = re.search('^(.*)\n.*$', str(executable, "utf-8"))
    if match:
        executable = match.group(1)
    return executable


def _find_reads(path):
    import os
    import re

    read_list = []
    for file in os.listdir(path):
        is_read = re.search('(.*)(\.fastq(?:\.gz)?)$', file, re.IGNORECASE)
        if is_read:
            sample_name = is_read.group(1)
            is_paired = re.search('^(.*)(_[R]?)([12])(.*)$', sample_name, re.IGNORECASE)
            if is_paired:
                if is_paired.group(3) == '1':  # If paired, only process read 1, so we don't double count the pair
                    sample_name = is_paired.group(1)
                    read1 = file
                    read2 = "%s%s2%s%s" % (is_paired.group(1), is_paired.group(2), is_paired.group(4), is_read.group(2))
                    if os.path.exists(os.path.join(path, read2)):
                        read = (sample_name, os.path.join(path, read1), os.path.join(path, read2))
                        read_list.append(read)
                        logging.info(read)
                    else:
                        logging.warning("Cannot find %s, the matching read to %s. Skipping..." % (read2, read1))
            else:
                read = (sample_name, os.path.join(path, file))
                read_list.append(read)
                logging.info(read)
    return read_list


def _get_bams(cwd):
    import re

    bam_list = []
    response = input("\nDo you have pre-aligned SAM/BAM files you wish to include [N]? ")
    if re.match('^[Yy]', response):
        path = input("Where are these files located [%s]? " % cwd)
        path = _expand_path(path) if path else cwd
        logging.info("Looking for sams/bams in %s...", path)
        bam_list = _find_files(path, "bam")
        bam_list.extend(_find_files(path, "sam"))
    return bam_list


def _get_vcfs(cwd):
    import re

    vcf_list = []
    response = input("\nDo you have pre-called VCFfiles you wish to include [N]? ")
    if re.match('^[Yy]', response):
        path = input("Where are these files located [%s]? " % cwd)
        path = _expand_path(path) if path else cwd
        logging.info("Looking for vcfs in %s...", path)
        vcf_list = _find_files(path, "vcf")
    return vcf_list


def _get_external_fastas(cwd, exclude):
    import re

    fasta_list = []
    response = input("\nDo you have fasta files for external genomes you wish to include [Y]? ")
    if not re.match('^[Nn]', response):
        path = input("Where are these files located [%s]? " % cwd)
        path = _expand_path(path) if path else cwd
        logging.info("Looking for external fastas in %s...", path)
        fasta_list = _find_files(path, "fasta")
        # remove the reference from the fasta_list if it is in there
        for fasta_tuple in fasta_list:
            fasta = fasta_tuple[1]
            if fasta == exclude:
                fasta_list.remove(fasta_tuple)
    return fasta_list


def _get_reads(cwd):
    import re

    read_list = []
    response = input("\nDo you have read files you wish to include [Y]? ")
    if not re.match('^[Nn]', response):
        path = input("Where are these files located [%s]? " % cwd)
        path = _expand_path(path) if path else cwd
        logging.info("Looking for read files in %s...", path)
        read_list = _find_reads(path)
    return read_list


def _get_application_path(application):
    import os

    app_path = _find_executable(application)
    if not app_path:
        app_path = input("\nUnable to find '%s', please enter the full path to '%s': " % (application, application))
    while not os.access(app_path, os.X_OK):
        app_path = input(
            "\n'%s' either does not exist or you don't have permission to run it, please enter the full path to '%s': " % (
            app_path, application))
    return app_path


def _get_java_path(jarfile):
    import os
    import fnmatch

    paths = ['/usr/share/java/']
    paths.extend(os.environ['PATH'].split(os.pathsep))
    for path in paths:
        if os.path.exists(os.path.join(path, jarfile)):
            return os.path.join(path, jarfile)
    # Didn't find it in path, check user's home directory
    for path, dirs, files in os.walk(os.path.expanduser("~")):
        for filename in fnmatch.filter(files, jarfile):
            return os.path.join(path, filename)
    #Didn't find it there, recursively check current directory
    for path, dirs, files in os.walk(os.getcwd()):
        for filename in fnmatch.filter(files, jarfile):
            return os.path.join(path, filename)
    #Let's ask the user
    jar_path = input("\nUnable to find '%s', please enter the full path to '%s': " % (jarfile, jarfile))
    while not os.access(jar_path, os.R_OK):
        jar_path = input(
            "\n'%s' either does not exist or you don't have permission to run it, please enter the full path to '%s': " % (
            jar_path, jarfile))
    return jar_path


def _get_advanced_settings(app_name, app_path, app_args, job_parms):
    import re
    import os

    response = input("Would you like to set advanced %s settings [N]? " % app_name)
    if re.match('^[Yy]', response):
        alt_version = input("  Would you like to use an alternate %s version [N]? " % app_name)
        if re.match('^[Yy]', alt_version):
            path = input("    What is the full path to the %s runtime you wish to use [%s]? " % (app_name, app_path))
            path = path if path else app_path
            while not os.access(path, os.X_OK):
                path = input(
                    "    '%s' either does not exist or you don't have permission to run it, please enter the full path to '%s': " % (
                    path, app_name))
            app_path = path
        args = input("  What additional arguments would you like to pass to %s [%s]? " % (app_name, app_args))
        if args:
            app_args = args
        queue = input("  What queue/partition should %s run on [%s]? " % (app_name, job_parms['queue']))
        if queue:
            job_parms['queue'] = queue
        mem = input("  How much memory (GB) will %s require [%s]? " % (app_name, job_parms['mem_requested']))
        if re.match("^[1-9][0-9]*$", mem):
            job_parms['mem_requested'] = mem
        cpus = input("  How many CPUs do you want %s to use [%s]? " % (app_name, job_parms['num_cpus']))
        if re.match("^[1-9][0-9]*$", cpus):
            job_parms['num_cpus'] = cpus
        hours = input("  How many hours will %s take to run [%s]? " % (app_name, job_parms['walltime']))
        if re.match("^[1-9][0-9]*$", hours):
            job_parms['walltime'] = hours
    return app_name, app_path, app_args, job_parms


def _get_aligners(queue, args):
    import re

    aligner_list = []
    bwa_path = ""
    print(
        "\nThis pipeline currently supports three aligners: BWA, Novoalign, and SNAP.\nYou can also provide pre-aligned BAM files, and you can choose as many options as you want.")
    response = input("\nWould you like to run BWA samp/se [N]?* ")
    if re.match('^[Yy]', response):
        bwa_path = _get_application_path("bwa")
        bwa_sampe_settings = _get_advanced_settings("BWA-sampe", bwa_path, "",
                                                    {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                     'queue': queue, 'args': args})
        aligner_list.append(bwa_sampe_settings)
        logging.info(bwa_sampe_settings)
    response = input("\nWould you like to run BWA mem [Y]? ")
    if not re.match('^[Nn]', response):
        if not bwa_path:
            bwa_path = _get_application_path("bwa")
        bwa_mem_settings = _get_advanced_settings("BWA-mem", bwa_path, "",
                                                  {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                   'queue': queue, 'args': args})
        aligner_list.append(bwa_mem_settings)
        logging.info(bwa_mem_settings)
    response = input("\nWould you like to run Bowtie2 [Y]? ")
    if not re.match('^[Nn]', response):
        bt2_path = _get_application_path("bowtie2")
        bt2_settings = _get_advanced_settings("Bowtie2", bt2_path, "",
                                              {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36', 'queue': queue,
                                               'args': args})
        aligner_list.append(bt2_settings)
        logging.info(bt2_settings)
    response = input("\nWould you like to run Novoalign [Y]? ")
    if not re.match('^[Nn]', response):
        novo_path = _get_application_path("novoalign")
        novo_settings = _get_advanced_settings("Novoalign", novo_path, "-r all",
                                               {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                'queue': queue, 'args': args})
        aligner_list.append(novo_settings)
        logging.info(novo_settings)
    response = input("\nWould you like to run SNAP [N]?* ")
    if re.match('^[Yy]', response):
        snap_path = _get_application_path("snap")
        snap_settings = _get_advanced_settings("SNAP", snap_path, "",
                                               {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                'queue': queue, 'args': args})
        aligner_list.append(snap_settings)
        logging.info(snap_settings)
    return aligner_list


def _get_snpcallers(queue, args):
    import re

    snpcaller_list = []
    using_gatk = False
    print(
        "\nThis pipeline currently supports four SNP callers: GATK, SolSNP, VarScan, and SAMtools, and you can provide VCF files.\nYou can choose as many options as you want.")
    response = input("\nWould you like to run GATK [Y]? ")
    if not re.match('^[Nn]', response):
        gatk_path = _get_java_path("GenomeAnalysisTK.jar")
        gatk_settings = _get_advanced_settings("GATK", gatk_path, "-stand_call_conf 100 -stand_emit_conf 100 -ploidy 1",
                                               {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                'queue': queue, 'args': args})
        snpcaller_list.append(gatk_settings)
        logging.info(gatk_settings)
        using_gatk = True
    response = input("\nWould you like to run SolSNP [N]? ")
    if re.match('^[Yy]', response):
        solsnp_path = _get_java_path("SolSNP.jar")
        solsnp_settings = _get_advanced_settings("SolSNP", solsnp_path, "",
                                                 {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                  'queue': queue, 'args': args})
        snpcaller_list.append(solsnp_settings)
        logging.info(solsnp_settings)
    response = input("\nWould you like to run VarScan [Y]? ")
    if not re.match('^[Nn]', response):
        varscan_path = _get_java_path("VarScan.jar")
        varscan_settings = _get_advanced_settings("VarScan", varscan_path, "",
                                                  {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                   'queue': queue, 'args': args})
        snpcaller_list.append(varscan_settings)
        logging.info(varscan_settings)
    response = input("\nWould you like to run SAMtools [Y]? ")
    if not re.match('^[Nn]', response):
        samtools_path = _get_application_path("bcftools")
        samtools_settings = _get_advanced_settings("SAMtools", samtools_path, "",
                                                   {'num_cpus': '4', 'mem_requested': '10', 'walltime': '36',
                                                    'queue': queue, 'args': args})
        snpcaller_list.append(samtools_settings)
        logging.info(samtools_settings)
    return snpcaller_list, using_gatk


def _get_job_submitter():
    import re

    job_submitter = "invalid"
    queue = ""
    args = ""
    response = input(
        "\nWhat system do you use for job management (PBS/TORQUE, SLURM, SGE*, and 'none' are currently supported) [PBS]? ")
    while job_submitter == "invalid":
        if re.match('^(PBS|Torque|qsub)$', response, re.IGNORECASE) or response == "":
            job_submitter = "PBS"
        elif re.match('^(SLURM|sbatch)$', response, re.IGNORECASE):
            job_submitter = "SLURM"
        elif re.match('^SGE', response, re.IGNORECASE):
            job_submitter = "SGE"
        elif re.match('^none$', response, re.IGNORECASE):
            job_submitter = "NONE"
        else:
            response = input("  %s is not a valid job management system, please enter another [PBS]? " % response)
    if job_submitter != "NONE":
        queue = input(
            "  Would you like to specify a queue/partition to use for all jobs (leave blank to use default queue) []? ")
        args = input("  What additional arguments do you need to pass to the job management system []? ")
    return job_submitter, queue, args


def _get_user_input(reference, output_folder):
    import os
    import re
    import sys

    configuration = {}
    cwd = os.getcwd()

    print("Welcome to the experimental python NASP version %s." % nasp_version)
    print("* Starred features are less well tested, and may not work.")

    if not output_folder:
        output_folder = input("\nWhere would you like output files to be written [nasp_results]? ")
        if not output_folder:
            output_folder = "nasp_results"
    output_folder = _expand_path(output_folder)
    if os.path.exists(output_folder):
        response = input(
            "\nOutput folder %s already exists!\nFiles in it may be overwritten!\nShould we continue anyway [N]? " % output_folder)
        if not re.match('^[Yy]', response):
            print("Operation cancelled!")
            quit()
    else:
        os.makedirs(output_folder)
    configuration["output_folder"] = output_folder

    logfile = os.path.join(output_folder, "runlog.txt")
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S',
                        filename=logfile,
                        filemode='w')

    if not reference:
        reference = input("\nWhere is the reference fasta file you would like to use? ")
    reference = _expand_path(reference)
    if not os.path.exists(reference):
        print("\nCannot continue because reference file %s does not seem to exist!" % reference)
        quit()
    configuration["reference"] = _create_file_tuple(reference)
    logging.info("Reference = %s", configuration["reference"])

    response = input(
        "\nDo you want to check the reference for duplicated regions\nand skip SNPs that fall in those regions [Y]? ")
    configuration["find_dups"] = "False" if re.match('^[Nn]', response) else "True"
    logging.info("FindDups = %s", configuration["find_dups"])

    (job_submitter, queue, args) = _get_job_submitter()
    configuration["job_submitter"] = job_submitter
    logging.info("JobSubmitter = %s", configuration["job_submitter"])

    name_match = re.search('^.*/(.*)$',
                           output_folder)  # Warning, not OS-independent! Should find a better way to do this.
    configuration["run_name"] = name_match.group(
        1)  # Temporary: setting the run name to be whatever the the output folder is named. Should ask user.
    logging.info("RunName = %s", configuration["run_name"])

    samtools_path = _find_executable("samtools")
    configuration["samtools"] = ("Samtools", samtools_path, "", {})
    logging.info("Samtools = %s", configuration["samtools"])

    run_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    configuration["index"] = ("Index", run_path, "",
                              {'name': 'nasp_index', 'num_cpus': '1', 'mem_requested': '2', 'walltime': '4',
                               'queue': queue, 'args': args})
    logging.info("Index = %s", configuration["index"])

    fasta_list = _get_external_fastas(cwd, reference)
    configuration["assemblies"] = fasta_list

    if configuration["find_dups"] or len(fasta_list) > 0:
        nucmer_path = _get_application_path("nucmer")
        nucmer_args = ""
        if len(fasta_list) > 0:
            deltafilter_path = _get_application_path("delta-filter")
            deltafilter_args = ""
            response = input("  Would you like to set advanced NUCmer settings [N]? ")
            if re.match('^[Yy]', response):
                nucmer_args = input(
                    "  What additional arguments would you like to pass to 'nucmer' while importing external genomes? ")
                deltafilter_args = input(
                    "  What additional arguments would you like to pass to 'delta-filter' while importing external genomes? ")
            configuration["assembly_importer"] = ("AssemblyImporter", deltafilter_path, deltafilter_args,
                                                  {'num_cpus': '1', 'mem_requested': '2', 'walltime': '1',
                                                   'queue': queue, 'args': args})
            logging.info("AssemblyImporter = %s", configuration["assembly_importer"])
        configuration["dup_finder"] = ("DupFinder", nucmer_path, nucmer_args,
                                       {'num_cpus': '1', 'mem_requested': '2', 'walltime': '1', 'queue': queue,
                                        'args': args})
        logging.info("DupFinder = %s", configuration["dup_finder"])

    read_list = _get_reads(cwd)
    configuration["reads"] = read_list

    if len(read_list) > 0:
        logging.info("Getting Aligners...")
        configuration["aligners"] = _get_aligners(queue, args)
    else:
        configuration["aligners"] = []

    bam_list = _get_bams(cwd)
    configuration["alignments"] = bam_list
    if len(bam_list) > 0:
        configuration["bam_index"] = ("BamIndex", run_path, "",
                                      {'name': 'nasp_bamindex', 'num_cpus': '1', 'mem_requested': '2', 'walltime': '4',
                                       'queue': queue, 'args': args})
        logging.info("BamIndex = %s", configuration["bam_index"])

    if len(read_list) > 0 or len(bam_list) > 0:
        logging.info("Getting SNP Callers...")
        (configuration["snpcallers"], using_gatk) = _get_snpcallers(queue, args)

        if using_gatk:
            picard_path = _get_java_path("CreateSequenceDictionary.jar")
            configuration["picard"] = ("Picard", os.path.dirname(picard_path), "", {})
            logging.info("Picard = %s", configuration["picard"])
    else:
        configuration["snpcallers"] = []

    vcf_list = _get_vcfs(cwd)
    configuration["vcfs"] = vcf_list

    if len(read_list) > 0 or len(bam_list) > 0 or len(vcf_list) > 0:
        coverage_filter = input(
            "\nThis pipeline can do filtering based on coverage.\nIf you do not want filtering based on coverage, enter 0.\nWhat is your minimum coverage threshold [10]? ")
        if not coverage_filter:
            coverage_filter = 10
        configuration["coverage_filter"] = str(coverage_filter)
        logging.info("CoverageFilter = %s", configuration["coverage_filter"])

        proportion_filter = input(
            "\nThis pipeline can do filtering based on the proportion of reads that match the call made by the SNP caller.\nIf you do not want filtering based on proportion, enter 0.\nWhat is the minimum acceptable proportion [0.9]? ")
        if not proportion_filter:
            proportion_filter = 0.9
        configuration["proportion_filter"] = str(proportion_filter)
        logging.info("ProportionFilter = %s", configuration["proportion_filter"])

    matrix_path = os.path.join(run_path, "vcf_to_matrix")
    if not os.path.exists(matrix_path):
        matrix_path = "vcf_to_matrix"
    print();
    matrix_settings = _get_advanced_settings("MatrixGenerator", matrix_path, "", {'name':'nasp_matrix', 'num_cpus':'12', 'mem_requested':'45', 'walltime':'48', 'queue':queue, 'args':args})
    configuration["matrix_generator"] = matrix_settings
    logging.info("MatrixGenerator = %s", configuration["matrix_generator"])

#    include_missing = input("\nDo you want to allow uncalled and filtered positions in the filtered matrix [N]? ")
#    if re.match('^[Yy]', include_missing):
#        configuration["filter_matrix_format"] = "missingdata"
#        logging.info("FilterMatrixFormat = %s", configuration["filter_matrix_format"])

    return configuration


def main():
    import nasp.dispatcher as dispatcher
    import nasp.configuration_parser as configuration_parser
    import os
    import re

    commandline_args = _parse_args()
    if commandline_args.config:
        configuration = configuration_parser.parse_config(commandline_args.config)
        output_folder = configuration["output_folder"]
        if os.path.exists(output_folder):
            response = input(
                "\nOutput folder %s already exists!\nFiles in it may be overwritten!\nShould we continue anyway [N]? " % output_folder)
            if not re.match('^[Yy]', response):
                print("Operation cancelled!")
                quit()
        else:
            os.makedirs(output_folder)
        logfile = os.path.join(output_folder, "runlog.txt")
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S',
                            filename=logfile,
                            filemode='w')
    else:
        configuration = _get_user_input(commandline_args.reference_fasta, commandline_args.output_folder)
    configuration_parser.write_config(configuration)
    dispatcher.begin(configuration)


if __name__ == "__main__":
    main()
