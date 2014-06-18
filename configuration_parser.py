#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.5"
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

# from UserDict import UserDict
# 
# class JobParameters(UserDict):
#     "store job scheduling parameters"
#     def __init__(self, scheduler="qsub"):
#         UserDict.__init__(self)
#         self["scheduler"] = scheduler
#     def __setitem__(self, key, item):
#         self.data[key] = item
#     def __getitem__(self, key):
#         return self.data[key]
    
def _parse_args():
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--config", required=True, help="Path to the configuration xml file." )
    return parser.parse_args()

def _parse_options( options_node ):
    configuration["run_name"] = options_node.find('RunName').text
    configuration["output_folder"] = options_node.find('OutputFolder').text
    reference_node = options_node.find('Reference')
    configuration["reference"] = (reference_node.get('name'), reference_node.get('path'))
    configuration["find_dups"] = reference_node.find('FindDups').text
    filter_node = options_node.find('Filters')
    configuration["coverage_filter"] = filter_node.find('CoverageFilter').text
    configuration["proportion_filter"] = filter_node.find('ProportionFilter').text
    configuration["job_submitter"] = options_node.find('JobSubmitter').text
    if options_node.find('FilterMatrixFormat'):
        configuration["filter_matrix_format"] = options_node.find('FilterMatrixFormat').text

def _find_reads( folder, filepath ):
    import os
    import re
    num_reads = 0
    folder.text = "\n\t\t\t"
    for file in os.listdir(filepath):
        is_read = re.search('(.*)(\.fastq(?:\.gz)?)$', file, re.IGNORECASE)
        if is_read:
            sample_name = is_read.group(1)
            is_paired = re.search('^(.*)(_[R]?)([12])(.*)$', sample_name, re.IGNORECASE)
            if is_paired:
                if is_paired.group(3) == '1': #If paired, only process read 1, so we don't double count the pair
                    sample_name = is_paired.group(1)
                    read1 = file
                    read2 = "%s%s2%s%s" % (is_paired.group(1), is_paired.group(2), is_paired.group(4), is_read.group(2))
                    if os.path.exists(os.path.join(filepath, read2)):
                        read_node = ElementTree.SubElement(folder, "ReadPair", {'sample':sample_name})
                        read_node.text = "\n\t\t\t\t"
                        read_node.tail = "\n\t\t\t"
                        read1_node = ElementTree.SubElement(read_node, "Read1Filename")
                        read1_node.text = read1
                        read1_node.tail = "\n\t\t\t\t"
                        read2_node = ElementTree.SubElement(read_node, "Read2Filename")
                        read2_node.text = read2                    
                        read2_node.tail = "\n\t\t\t"
                        read_list.append((sample_name, os.path.join(filepath, read1), os.path.join(filepath, read2)))
                        num_reads += 1
                    else:
                        print("Cannot find %s, the matching read to %s. Skipping..." % (read2, read1))
            else:
                read_node = ElementTree.SubElement(folder, "Read", {'sample':sample_name})
                read_node.text = file
                read_node.tail = "\n\t\t\t"
                read_list.append((sample_name, os.path.join(filepath, file)))
                num_reads += 1
    return num_reads

def _find_files( folder, filepath, filetype, extension ):
    import os
    import re
    num_files = 0
    folder.text = "\n\t\t\t"
    for file in os.listdir(filepath):
        is_type = re.search('(.*)(\.%s)$' % extension, file, re.IGNORECASE)
        if is_type:
            sample_name = is_type.group(1)
            file_node = ElementTree.SubElement(folder, filetype, {'sample':sample_name})
            file_node.text = file
            file_node.tail = "\n\t\t\t"
            if filetype == 'Assembly':
                fasta_list.append((sample_name, os.path.join(filepath, file)))
            if filetype == 'Alignment':
                bam_list.append((sample_name, os.path.join(filepath, file)))
            if filetype == 'VCFFile':
                vcf_list.append((sample_name, os.path.join(filepath, file)))
            num_files += 1
    return num_files

def _get_reads( folder ):
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

def _get_fastas( folder ):
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

def _get_bams( folder ):
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

def _get_vcfs( folder ):
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

def _parse_files( files_node ):
    for folder in files_node.findall('ReadFolder'):
        _get_reads( folder )
    configuration["reads"] = read_list
    for folder in files_node.findall('AssemblyFolder'):
        _get_fastas( folder )
    configuration["assemblies"] = fasta_list
    for folder in files_node.findall('AlignmentFolder'):
        _get_bams( folder )
    configuration["alignments"] = bam_list
    for folder in files_node.findall('VCFFolder'):
        _get_vcfs( folder )
    configuration["vcfs"] = vcf_list

def _get_application( app_node, name=None ):
    name = name or app_node.get('name')
    path = app_node.get('path')
    args = app_node.find('AdditionalArgs').text or ""
    job_parms = {}
    job_node = app_node.find('JobParameters')
    if ( job_node is not None ):
        job_parms['name'] = job_node.get("name")
        job_parms['mem_requested'] = job_node.find('MemRequested').text
        job_parms['num_cpus'] = job_node.find('NumCPUs').text
        job_parms['walltime'] = job_node.find('Walltime').text
    return (name, path, args, job_parms)

def _parse_applications( applications_node ):
    configuration["picard"] = _get_application(applications_node.find('Picard'), "picard")
    configuration["samtools"] = _get_application(applications_node.find('Samtools'), "samtools")
    configuration["dup_finder"] = _get_application(applications_node.find('DupFinder'), "dupFinder")
    configuration["assembly_importer"] = _get_application(applications_node.find('AssemblyImporter'), "assemblyImporter")
    for aligner in applications_node.findall('Aligner'):
        aligner_list.append(_get_application(aligner))
    configuration["aligners"] = aligner_list
    for snpcaller in applications_node.findall('SNPCaller'):
        snpcaller_list.append(_get_application(snpcaller))
    configuration["snpcallers"] = snpcaller_list
        
def parse_config( config_file ):
    xmltree = ElementTree.parse( config_file )
    root = xmltree.getroot()
    _parse_options(root.find('Options'))
    _parse_files(root.find('Files'))
    _parse_applications(root.find('ExternalApplications'))
    ElementTree.ElementTree(root).write("%s_temp.xml" % config_file)
    return configuration

def main():
    commandline_args = _parse_args()
    parse_config( commandline_args.config )

if __name__ == "__main__": main()
