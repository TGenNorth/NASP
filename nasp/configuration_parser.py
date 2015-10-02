#!/usr/bin/env python3
# coding=utf-8

__author__ = "Darrin Lemmer"
__version__ = "1.0.0"
__email__ = "dlemmer@tgen.org"

'''
Created on Feb 18, 2014

@author: dlemmer
'''

import logging
from xml.etree import ElementTree

configuration = {}
read_list = []
fasta_list = []
bam_list = []
vcf_list = []
aligner_list = []
snpcaller_list = []


def _parse_args():
    """
    Returns:
    """
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--config", required=True, help="Path to the configuration xml file.")
    return parser.parse_args()


def _parse_options(options_node):
    """
    Args:
        options_node (xml.etree.ElementTree.Element):
    """
    configuration["run_name"] = options_node.findtext('RunName')
    configuration["output_folder"] = options_node.findtext('OutputFolder')
    reference_node = options_node.find('Reference')
    configuration["reference"] = (reference_node.get('name'), reference_node.get('path'))
    configuration["find_dups"] = reference_node.findtext('FindDups')
    filter_node = options_node.find('Filters')
    if filter_node.find('CoverageFilter') is not None:
        configuration["coverage_filter"] = filter_node.findtext('CoverageFilter')
    if filter_node.find('ProportionFilter') is not None:
        configuration["proportion_filter"] = filter_node.findtext('ProportionFilter')
    configuration["job_submitter"] = options_node.findtext('JobSubmitter')
    if options_node.find('FilterMatrixFormat') is not None:
        configuration["filter_matrix_format"] = options_node.findtext('FilterMatrixFormat')
    if options_node.find('TrimReads') is not None:
        configuration["trim_reads"] = options_node.findtext('TrimReads')


def _find_reads(folder, filepath):
    """
    Args:
        folder (xml.etree.ElementTree.Element): The Folder parent Element
        filepath (str): Path to the folder where the fasta files are located

    Returns:
        int: Number of unique read files found, not including paired files
    """
    import os
    import re

    num_reads = 0
    for file in os.listdir(filepath):
        is_read = re.search('(.*)(\.fastq(?:\.gz)?)$', file, re.IGNORECASE)
        if is_read:
            sample_name = is_read.group(1)
            is_paired = re.search('^(.*)(_[R]?)([12])(.*)$', sample_name, re.IGNORECASE)
            if is_paired:
                if is_paired.group(3) == '1':  # If paired, only process read 1, so we don't double count the pair
                    sample_name = is_paired.group(1)
                    read1 = file
                    read2 = "%s%s2%s%s" % (is_paired.group(1), is_paired.group(2), is_paired.group(4), is_read.group(2))
                    if os.path.exists(os.path.join(filepath, read2)):
                        read_node = ElementTree.SubElement(folder, "ReadPair", {'sample': sample_name})
                        read1_node = ElementTree.SubElement(read_node, "Read1Filename")
                        read1_node.text = read1
                        read2_node = ElementTree.SubElement(read_node, "Read2Filename")
                        read2_node.text = read2
                        read_list.append((sample_name, os.path.join(filepath, read1), os.path.join(filepath, read2)))
                        num_reads += 1
                    else:
                        print("Cannot find %s, the matching read to %s. Skipping..." % (read2, read1))
            else:
                read_node = ElementTree.SubElement(folder, "Read", {'sample': sample_name})
                read_node.text = file
                read_list.append((sample_name, os.path.join(filepath, file)))
                num_reads += 1
    return num_reads


def _find_files(folder, filepath, filetype, extension):
    """
    Args:
        folder:
        filepath:
        filetype:
        extension:

    Returns:
        int:
    """
    import os
    import re

    num_files = 0
    for file in os.listdir(filepath):
        is_type = re.search('(.*)(\.%s)$' % extension, file, re.IGNORECASE)
        if is_type:
            sample_name = is_type.group(1)
            file_node = ElementTree.SubElement(folder, filetype, {'sample': sample_name})
            file_node.text = file
            if filetype == 'Assembly':
                fasta_list.append((sample_name, os.path.join(filepath, file)))
            if filetype == 'Alignment':
                bam_list.append((sample_name, os.path.join(filepath, file)))
            if filetype == 'VCFFile':
                vcf_list.append((sample_name, os.path.join(filepath, file)))
            num_files += 1
    return num_files


def _get_reads(folder):
    """
    Args:
        folder (xml.etree.ElementTree.Element):

    Returns:
        int:
    """
    from os import path

    filepath = folder.get('path')
    num_reads = 0
    for read in folder.findall('Read'):
        SEread = (read.get('sample'), path.join(filepath, read.text))
        read_list.append(SEread)
        num_reads += 1
    for readpair in folder.findall('ReadPair'):
        read1 = path.join(filepath, readpair.find('Read1Filename').text)
        read2 = path.join(filepath, readpair.find('Read2Filename').text)
        PEread = (readpair.get('sample'), read1, read2)
        read_list.append(PEread)
        num_reads += 1
    if num_reads == 0:
        num_reads = _find_reads(folder, filepath)
    return num_reads


def _get_fastas(folder):
    """
    Args:
        folder (xml.etree.ElementTree.Element):

    Returns:
        int:
    """
    from os import path

    filepath = folder.get('path')
    num_files = 0
    for assembly in folder.findall('Assembly'):
        fasta = (assembly.get('sample'), path.join(filepath, assembly.text))
        fasta_list.append(fasta)
        num_files += 1
    if num_files == 0:
        num_files = _find_files(folder, filepath, 'Assembly', 'fasta')
    return num_files


def _get_bams(folder):
    """
    Args:
        folder (xml.etree.ElementTree.Element):

    Returns:
        int:
    """
    from os import path

    filepath = folder.get('path')
    num_files = 0
    for alignment in folder.findall('Alignment'):
        bam = (alignment.get('sample'), path.join(filepath, alignment.text))
        bam_list.append(bam)
        num_files += 1
    if num_files == 0:
        num_files = _find_files(folder, filepath, 'Alignment', 'bam')
    return num_files


def _get_vcfs(folder):
    """
    Args:
        folder (xml.etree.ElementTree.Element):

    Returns:
        int:
    """
    from os import path

    filepath = folder.get('path')
    num_files = 0
    for vcffile in folder.findall('VCFFile'):
        vcf = (vcffile.get('sample'), path.join(filepath, vcffile.text))
        vcf_list.append(vcf)
        num_files += 1
    if num_files == 0:
        num_files = _find_files(folder, filepath, 'VCFFile', 'vcf')
    return num_files


def _parse_files(files_node):
    """
    Args:
        files_node (xml.etree.ElementTree.Element):
    """
    for folder in files_node.findall('ReadFolder'):
        _get_reads(folder)
    configuration["reads"] = read_list
    for folder in files_node.findall('AssemblyFolder'):
        _get_fastas(folder)
    configuration["assemblies"] = fasta_list
    for folder in files_node.findall('AlignmentFolder'):
        _get_bams(folder)
    configuration["alignments"] = bam_list
    for folder in files_node.findall('VCFFolder'):
        _get_vcfs(folder)
    configuration["vcfs"] = vcf_list


def _get_application(app_node, name=None):
    """
    Args:
        app_node (xml.etree.ElementTree.Element):
        name (str):

    Returns:
        (str, str, str, dict): job name, command path, command arguments, dispatcher parameters
    """
    name = name or app_node.get('name')
    path = app_node.get('path')
    args = app_node.findtext('AdditionalArguments', default="")
    job_parms = {}
    job_node = app_node.find('JobParameters')
    if job_node is not None:
        job_parms['name'] = job_node.get("name") if job_node.get("name") else ""
        job_parms['mem_requested'] = job_node.findtext('MemRequested', default="")
        job_parms['num_cpus'] = job_node.findtext('NumCPUs', default="")
        job_parms['walltime'] = job_node.findtext('Walltime', default="")
        job_parms['queue'] = job_node.findtext('Queue', default="")
        job_parms['args'] = job_node.findtext('JobSubmitterArgs', default="")
    return name, path, args, job_parms


def _parse_applications(applications_node):
    """
    Args:
        applications_node (xml.etree.ElementTree.Element):
    """
    if applications_node.find('Picard'):
        configuration["picard"] = _get_application(applications_node.find('Picard'), "Picard")
    configuration["samtools"] = _get_application(applications_node.find('Samtools'), "Samtools")
    if applications_node.find('ReadTrimmer'):
        configuration["read_trimmer"] = _get_application(applications_node.find('ReadTrimmer'), "ReadTrimmer")
    if applications_node.find('DupFinder'):
        configuration["dup_finder"] = _get_application(applications_node.find('DupFinder'), "DupFinder")
    if applications_node.find('Index'):
        configuration["index"] = _get_application(applications_node.find('Index'), "Index")
    if applications_node.find('BamIndex'):
        configuration["bam_index"] = _get_application(applications_node.find('BamIndex'), "BamIndex")
    configuration["matrix_generator"] = _get_application(applications_node.find('MatrixGenerator'), "MatrixGenerator")
    if applications_node.find('AssemblyImporter'):
        configuration["assembly_importer"] = _get_application(applications_node.find('AssemblyImporter'),
                                                              "AssemblyImporter")
    for aligner in applications_node.findall('Aligner'):
        aligner_list.append(_get_application(aligner))
    configuration["aligners"] = aligner_list
    for snpcaller in applications_node.findall('SNPCaller'):
        snpcaller_list.append(_get_application(snpcaller))
    configuration["snpcallers"] = snpcaller_list


def _write_reads(node, read_list):
    """
    Args:
        node:
        read_list:

    Returns:
    """
    import os
    from collections import defaultdict

    read_dict = defaultdict(list)
    for read_tuple in read_list:
        (name, read1) = read_tuple[0:2]
        read2 = read_tuple[2] if len(read_tuple) >= 3 else ""
        folder = os.path.dirname(read1)
        read1 = os.path.basename(read1)
        if read2:
            read2 = os.path.basename(read2)
        read_dict[folder].append((name, read1, read2))
    for folder in read_dict:
        folder_node = ElementTree.SubElement(node, "ReadFolder", {'path': folder})
        for (name, read1, read2) in read_dict[folder]:
            if read2:  # We have paired reads
                read_node = ElementTree.SubElement(folder_node, "ReadPair", {'sample': name})
                read1_node = ElementTree.SubElement(read_node, "Read1Filename")
                read1_node.text = read1
                read2_node = ElementTree.SubElement(read_node, "Read2Filename")
                read2_node.text = read2
            else:
                read_node = ElementTree.SubElement(folder_node, "Read", {'sample': name})
                read_node.text = read1
    return node


def _write_files(node, file_list, foldernode, filenode):
    """
    Args:
        node:
        file_list:
        foldernode:
        filenode:

    Returns:
    """
    import os
    from collections import defaultdict

    file_dict = defaultdict(list)
    for (name, file) in file_list:
        folder = os.path.dirname(file)
        file = os.path.basename(file)
        file_dict[folder].append((name, file))
    for folder in file_dict:
        folder_node = ElementTree.SubElement(node, foldernode, {'path': folder})
        for (name, file) in file_dict[folder]:
            file_node = ElementTree.SubElement(folder_node, filenode, {'sample': name})
            file_node.text = file
    return node


def _write_application(node, details, app_type):
    """
    Args:
        node:
        details:
        app_type:

    Returns:
    """
    (name, path, args, job_parms) = details
    app_node = ElementTree.SubElement(node, app_type, {'name': name, 'path': path})
    arg_node = ElementTree.SubElement(app_node, "AdditionalArguments")
    arg_node.text = args
    if len(job_parms) > 0:
        job_node = ElementTree.SubElement(app_node, "JobParameters")
        if "name" in job_parms:
            job_node.set('name', job_parms["name"])
        ElementTree.SubElement(job_node, "MemRequested").text = job_parms[
            "mem_requested"] if "mem_requested" in job_parms else ""
        ElementTree.SubElement(job_node, "NumCPUs").text = job_parms["num_cpus"] if "num_cpus" in job_parms else ""
        ElementTree.SubElement(job_node, "Walltime").text = job_parms["walltime"] if "walltime" in job_parms else ""
        ElementTree.SubElement(job_node, "Queue").text = job_parms["queue"] if "queue" in job_parms else ""
        ElementTree.SubElement(job_node, "JobSubmitterArgs").text = job_parms["args"] if "args" in job_parms else ""
    return node


def _write_config_node(root, xml_file):
    """
    Args:
        root:
        xml_file:

    Returns:
    """
    from xml.dom import minidom

    dom = minidom.parseString(ElementTree.tostring(root, 'utf-8'))
    output = open(xml_file, 'w')
    output.write(dom.toprettyxml(indent="    ", newl="\n"))
    output.close()
    return xml_file


def write_config(configuration):
    """
    Args:
        configuration:
    """
    import os

    root = ElementTree.Element("NaspInputData")

    # Create the Options section
    options_node = ElementTree.SubElement(root, "Options")

    run_name = configuration["run_name"]
    node = ElementTree.SubElement(options_node, "RunName")
    node.text = run_name

    output_folder = configuration["output_folder"]
    node = ElementTree.SubElement(options_node, "OutputFolder")
    node.text = output_folder

    xml_file = os.path.join(output_folder, "%s-config.xml" % run_name)

    (name, path) = configuration["reference"]
    ref_node = ElementTree.SubElement(options_node, "Reference", {'name': name, 'path': path})
    node = ElementTree.SubElement(ref_node, "FindDups")
    node.text = configuration["find_dups"]

    filter_node = ElementTree.SubElement(options_node, "Filters")
    if "proportion_filter" in configuration:
        proportion = ElementTree.SubElement(filter_node, "ProportionFilter")
        proportion.text = configuration["proportion_filter"]
    if "coverage_filter" in configuration:
        coverage = ElementTree.SubElement(filter_node, "CoverageFilter")
        coverage.text = configuration["coverage_filter"]

    node = ElementTree.SubElement(options_node, "JobSubmitter")
    node.text = configuration["job_submitter"]

    if "filter_matrix_format" in configuration:
        node = ElementTree.SubElement(options_node, "FilterMatrixFormat")
        node.text = configuration["filter_matrix_format"]

    if "trim_reads" in configuration:
        node = ElementTree.SubElement(options_node, "TrimReads")
        node.text = configuration["trim_reads"]

    #Create the Files section
    files_node = ElementTree.SubElement(root, "Files")

    read_list = configuration["reads"]
    if len(read_list) > 0:
        _write_reads(files_node, read_list)

    fasta_list = configuration["assemblies"]
    if len(fasta_list) > 0:
        _write_files(files_node, fasta_list, "AssemblyFolder", "Assembly")

    bam_list = configuration["alignments"]
    if len(bam_list) > 0:
        _write_files(files_node, bam_list, "AlignmentFolder", "Alignment")

    vcf_list = configuration["vcfs"]
    if len(fasta_list) > 0:
        _write_files(files_node, vcf_list, "VCFFolder", "VCFFile")

    #Create the ExternalApplications section
    applications_node = ElementTree.SubElement(root, "ExternalApplications")

    _write_application(applications_node, configuration["index"], "Index")
    if "read_trimmer" in configuration:
        _write_application(applications_node, configuration["read_trimmer"], "ReadTrimmer")
    if "bam_index" in configuration:
        _write_application(applications_node, configuration["bam_index"], "BamIndex")
    _write_application(applications_node, configuration["matrix_generator"], "MatrixGenerator")
    if "picard" in configuration:
        _write_application(applications_node, configuration["picard"], "Picard")
    _write_application(applications_node, configuration["samtools"], "Samtools")
    if "dup_finder" in configuration:
        _write_application(applications_node, configuration["dup_finder"], "DupFinder")
    if "assembly_importer" in configuration:
        _write_application(applications_node, configuration["assembly_importer"], "AssemblyImporter")

    for aligner in configuration["aligners"]:
        _write_application(applications_node, aligner, "Aligner")

    for snpcaller in configuration["snpcallers"]:
        _write_application(applications_node, snpcaller, "SNPCaller")

    _write_config_node(root, xml_file)


def parse_config(config_file):
    """

    Args:
        config_file (str): path to an XML formatted configuration file

    Returns:
        dict: configuration dictionary
    """
    xmltree = ElementTree.parse(config_file)
    root = xmltree.getroot()
    _parse_options(root.find('Options'))
    _parse_files(root.find('Files'))
    _parse_applications(root.find('ExternalApplications'))
    # Commented part might be nicer way to do this here, but ElementTree keeps whitespace with each node, thus repeated reading
    # and pretty printing the same file keeps adding additional carriage returns. So we need to create new tree each time instead.
    #
    #    output_folder = configuration['output_folder']
    #    xml_file = os.path.join(output_folder, config_file)
    #    run_name = configuration['run_name']
    #    if run_name:
    #        xml_file = os.path.join(output_folder, "%s-config.xml" % run_name)
    #    _write_config_node( root, xml_file )
    #    write_config(configuration)
    return configuration


def main():
    commandline_args = _parse_args()
    parse_config(commandline_args.config)


if __name__ == "__main__":
    main()
