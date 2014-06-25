#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.5"
__email__ = "dsmith@tgen.org"

import logging


def _parse_args():
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--mode", required=True, choices=[ 'commandline', 'xml' ], help="Data passing mode, must be set to 'commandline' or 'xml'." )
    parser.add_argument( "--reference-fasta", help="Path to input reference fasta file." )
    parser.add_argument( "--reference-dups", help="Path to input reference dups file." )
    parser.add_argument( "--input-files", nargs="+", help="Path to input VCF/fasta files for matrix conversion." )
    parser.add_argument( "--master-matrix", default="master_matrix.tsv", help="Name of master matrix to create." )
    parser.add_argument( "--filter-matrix", default="filter_matrix.tsv", help="Name of custom matrix to create." )
    parser.add_argument( "--filter-matrix-format", help="String describing the custom format of the filter matrix." )
    parser.add_argument( "--general-stats", default="general_stats.tsv", help="Name of general statistics file to create." )
    parser.add_argument( "--sample-stats", default="sample_stats.tsv", help="Name of sample statistics file to create." )
    parser.add_argument( "--minimum-coverage", type=int, default=10, help="Minimum coverage depth at a position." )
    parser.add_argument( "--minimum-proportion", type=float, default=0.9, help="Minimum proportion of reads that must match the call at a position." )
    parser.add_argument( "--num-threads", type=int, default=1, help="Number of threads to use when processing input." )
    parser.add_argument( "--dto-file", help="Path to a matrix_dto XML file that defines all the parameters." )
    return parser.parse_args()

def _parse_input_config(commandline_args):
    import matrix_DTO
    (matrix_parms, input_files) = matrix_DTO.parse_dto(commandline_args.dto_file)
    commandline_args.reference_fasta = matrix_parms['reference-fasta']
    commandline_args.reference_dups = matrix_parms['reference-dups']
    commandline_args.master_matrix = matrix_parms['master-matrix']
    commandline_args.filter_matrix = matrix_parms['filter-matrix']
    commandline_args.general_stats = matrix_parms['general-stats']
    commandline_args.contig_stats = matrix_parms['contig-stats']
    commandline_args.minimum_coverage = int(matrix_parms['minimum-coverage'])
    commandline_args.minimum_proportion = float(matrix_parms['minimum-proportion'])
    if "filter-matrix-format" in matrix_parms:
        commandline_args.filter_matrix_format = matrix_parms['filter-matrix-format']
    commandline_args.input_files = input_files
    return commandline_args

def import_reference( reference, reference_path, dups_path ):
    reference.import_fasta_file( reference_path )
    if dups_path is not None:
        reference.import_dups_file( dups_path )
    #from sys import stdout
    #reference._genome._send_to_fasta_handle( stdout )
    #reference._dups._send_to_fasta_handle( stdout )

def import_external_fasta( input_file ):
    from nasp_objects import FastaGenome
    genome = FastaGenome()
    set_genome_metadata( genome, input_file )
    genome.import_fasta_file( genome.file_path(), "franken::" )
    #from sys import stdout
    #genome._genome._send_to_fasta_handle( stdout )
    return [ genome ]

# FIXME split into a larger number of smaller more testable functions
# FIXME This belongs in VCFGenome object perhaps?
def read_vcf_file( reference, min_coverage, min_proportion, input_file ):
    genomes = {}
    file_path = get_file_path( input_file )
    with open( file_path, 'r' ) as vcf_filehandle:
        from nasp_objects import VCFGenome, Genome, ReferenceCallMismatch, VCFRecord
        #import vcf
        vcf_record = VCFRecord( file_path )
        #vcf_data_handle = vcf.Reader( vcf_filehandle )
        vcf_samples = vcf_record.get_samples()
        #print( vcf_samples )
        for vcf_sample in vcf_samples:
            genomes[vcf_sample] = VCFGenome()
            set_genome_metadata( genomes[vcf_sample], input_file )
            genomes[vcf_sample].set_nickname( vcf_sample )
        while vcf_record.fetch_next_record():
            current_contig = vcf_record.get_contig()
            current_pos = vcf_record.get_position()
            if current_pos <= reference.get_contig_length( current_contig ):
                reference_call = reference.get_call( current_pos, None, current_contig )
                simplified_refcall = Genome.simple_call( reference_call )
                if ( simplified_refcall != 'N' ) and ( simplified_refcall != Genome.simple_call( vcf_record.get_reference_call()[0] ) ):
                    raise ReferenceCallMismatch( reference_call, vcf_record.get_reference_call(), file_path, current_contig, current_pos )
                for vcf_sample in vcf_samples:
                    sample_info = vcf_record.get_sample_info( vcf_sample )
                    # FIXME indels
                    if sample_info['call'] is not None:
                        genomes[vcf_sample].set_call( sample_info['call'], current_pos, 'X', current_contig )
                    if sample_info['was_called']:
                        genomes[vcf_sample].set_was_called( 'Y', current_pos, current_contig )
                    if sample_info['coverage'] is not None:
                        if sample_info['coverage'] >= min_coverage:
                            genomes[vcf_sample].set_coverage_pass( 'Y', current_pos, current_contig )
                        else:
                            genomes[vcf_sample].set_coverage_pass( 'N', current_pos, current_contig )
                    if sample_info['proportion'] is not None:
                        if sample_info['proportion'] >= min_proportion:
                            genomes[vcf_sample].set_proportion_pass( 'Y', current_pos, current_contig )
                        else:
                            genomes[vcf_sample].set_proportion_pass( 'N', current_pos, current_contig )
                    elif not sample_info['is_a_snp']:
                        genomes[vcf_sample].set_proportion_pass( '-', current_pos, current_contig )
    #from sys import stdout
    #for genome in genomes:
    #    genomes[genome]._genome._send_to_fasta_handle( stdout )
    return genomes.values()

# FIXME These three functions should be combined?
def determine_file_type( input_file ):
    import re
    file_type = None
    filename_match = re.match( r'^((?:[^,:]+,)+)::.*$', input_file )
    if filename_match:
        generator_array = filename_match.group(1).split( ',' )
        file_type = generator_array[0]
    #print( file_type )
    return file_type

def get_file_path( input_file ):
    import re
    file_path = None
    filename_match = re.match( r'^(?:[^,:]+,)+::(.*)$', input_file )
    if filename_match:
        file_path = filename_match.group(1)
    else:
        file_path = input_file
    #print( file_path )
    return file_path

def set_genome_metadata( genome, input_file ):
    import re
    filename_match = re.match( r'^((?:[^,:]+,)+)::(.*)$', input_file )
    if filename_match:
        generator_array = filename_match.group(1).split( ',' )
        genome.set_file_type( generator_array[0] )
        genome.add_generators( generator_array[1:-1] )
        genome.set_file_path( filename_match.group(2) )
    else:
        genome.set_file_path( input_file )
    #print( genome.identifier() )

def manage_input_thread( reference, min_coverage, min_proportion, input_q, output_q ):
    input_file = input_q.get()
    while input_file is not None:
        try:
            new_genomes = []
            file_type = determine_file_type( input_file )
            if file_type == "frankenfasta":
                new_genomes = import_external_fasta( input_file )
            elif file_type == "vcf":
                new_genomes = read_vcf_file( reference, min_coverage, min_proportion, input_file )
            for new_genome in new_genomes:
                output_q.put( new_genome )
        except:
            logging.exception( "Unable to read in data from '{0}'!".format( get_file_path( input_file ) ) )
        input_file = input_q.get()
    output_q.put( None )

def parse_input_files( input_files, num_threads, genomes, min_coverage, min_proportion ):
    from multiprocessing import Process, Queue
    #from queue import Queue
    from time import sleep
    input_q = Queue()
    output_q = Queue()
    for input_file in input_files:
        input_q.put( input_file )
    if num_threads > input_q.qsize():
        num_threads = input_q.qsize()
    sleep( 1 )
    thread_list = []
    for current_thread in range( num_threads ):
        input_q.put( None )
        current_thread = Process( target=manage_input_thread, args=[ genomes.reference(), min_coverage, min_proportion, input_q, output_q ] )
        current_thread.start()
        #manage_input_thread( genomes.reference(), min_coverage, min_proportion, input_q, output_q )
        thread_list.append( current_thread )
    sleep( 1 )
    while num_threads > 0:
        new_genome = output_q.get()
        if new_genome is None:
            num_threads -= 1
        else:
            genomes.add_genome( new_genome )
    sleep( 1 )
    for current_thread in thread_list:
        current_thread.join()

def write_output_matrices( genomes, master_matrix, filter_matrix, matrix_format ):
    genomes.write_to_matrices( master_matrix, filter_matrix, matrix_format )

def write_stats_data( genomes, general_stats, sample_stats ):
    genomes.write_to_stats_files( general_stats, sample_stats )


def main():
    commandline_args = _parse_args()
    if(commandline_args.dto_file):
        commandline_args = _parse_input_config(commandline_args)
    logging.basicConfig( level=logging.WARNING )
    from nasp_objects import ReferenceGenome, GenomeCollection
    reference = ReferenceGenome()
    import_reference( reference, commandline_args.reference_fasta, commandline_args.reference_dups )
    genomes = GenomeCollection()
    genomes.set_reference( reference )
    parse_input_files( commandline_args.input_files, commandline_args.num_threads, genomes, commandline_args.minimum_coverage, commandline_args.minimum_proportion )
    write_output_matrices( genomes, commandline_args.master_matrix, commandline_args.filter_matrix, commandline_args.filter_matrix_format )
    write_stats_data( genomes, commandline_args.general_stats, commandline_args.sample_stats )

if __name__ == "__main__": main()


