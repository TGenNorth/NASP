#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.2"
__email__ = "dlemmer@tgen.org"

'''
Created on Feb 18, 2014

@author: dlemmer
'''

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

def _get_reads( folder ):
    from os import path
    filepath = folder.get('path')
    for read in folder.findall('Read'):
        SEread = (read.get('sample'), path.join(filepath, read.text))
        read_list.append(SEread)
    for readpair in folder.findall('ReadPair'):
        read1 = path.join(filepath, readpair.find('Read1Filename').text)
        read2 = path.join(filepath, readpair.find('Read2Filename').text)
        PEread = (readpair.get('sample'), read1, read2)
        read_list.append(PEread)

def _get_fastas( folder ):
    from os import path
    filepath = folder.get('path')
    for assembly in folder.findall('Assembly'):
        fasta = (assembly.get('sample'), path.join(filepath, assembly.text))
        fasta_list.append(fasta)

def _get_bams( folder ):
    from os import path
    filepath = folder.get('path')
    for alignment in folder.findall('Alignment'):
        bam = (alignment.get('sample'), path.join(filepath, alignment.text))
        bam_list.append(bam)

def _get_vcfs( folder ):
    from os import path
    filepath = folder.get('path')
    for vcffile in folder.findall('VCFFile'):
        vcf = (vcffile.get('sample'), path.join(filepath, vcffile.text))
        vcf_list.append(vcf)

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
    args = app_node.find('AdditionalArgs').text
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
    configuration['aligners'] = aligner_list
    for snpcaller in applications_node.findall('SNPCaller'):
        snpcaller_list.append(_get_application(snpcaller))
    configuration['snpcallers'] = snpcaller_list
    
    
def parse_config( config_file ):
    from xml.etree import ElementTree
    xmltree = ElementTree.parse( config_file )
    root = xmltree.getroot()
    _parse_options(root.find('Options'))
    _parse_files(root.find('Files'))
    _parse_applications(root.find('ExternalApplications'))
    return configuration


def main():
    commandline_args = _parse_args()
    parse_config( commandline_args.config )


if __name__ == "__main__": main()
