#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.2"
__email__ = "dsmith@tgen.org"


def _parse_args():
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--reference-fasta", required=True, help="Path to input reference fasta file." )
    parser.add_argument( "--reference-dups", help="Path to input reference dups file." )
    parser.add_argument( "--input-files", nargs="+", required=True, help="Path to input VCF/fasta files for matrix conversion." )
    parser.add_argument( "--master-matrix", default="master_matrix.tsv", help="Name of master matrix to create." )
    parser.add_argument( "--general-stats", default="general_stats.tsv", help="Name of general statistics file to create." )
    parser.add_argument( "--contig-stats", default="contig_stats.tsv", help="Name of contig statistics file to create." )
    parser.add_argument( "--minimum-coverage", type=int, default=10, help="Minimum coverage depth at a position." )
    parser.add_argument( "--minimum-proportion", type=float, default=0.9, help="Minimum proportion of reads that must match the call at a position." )
    parser.add_argument( "--num-threads", type=int, default=1, help="Number of threads to use when processing input." )
    return( parser.parse_args() )

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

def read_vcf_file( input_file ):
    genomes = {}
    file_path = get_file_path( input_file )
    with open( file_path, 'r' ) as vcf_filehandle:
        from nasp_objects import VCFGenome
        import vcf
        vcf_data_handle = vcf.Reader( vcf_filehandle )
        vcf_samples = vcf_data_handle.samples
        #print( vcf_samples )
        for vcf_sample in vcf_samples:
            genomes[vcf_sample] = VCFGenome()
            set_genome_metadata( genomes[vcf_sample], input_file )
            genomes[vcf_sample].set_nickname( vcf_sample )
        for vcf_record in vcf_data_handle:
            pass
            #print( vcf_record, vcf_record.INFO )
    return genomes.values()

def determine_file_type( input_file ):
    import re
    file_type = None
    filename_match = re.match( r'^((?:[A-Za-z0-9._-]+,)+)::.*$', input_file )
    if filename_match:
        generator_array = filename_match.group(1).split( ',' )
        file_type = generator_array[0]
    #print( file_type )
    return file_type

def get_file_path( input_file ):
    import re
    file_path = None
    filename_match = re.match( r'^(?:[A-Za-z0-9._-]+,)+::(.*)$', input_file )
    if filename_match:
        file_path = filename_match.group(1)
    else:
        file_path = input_file
    #print( file_path )
    return file_path

def set_genome_metadata( genome, input_file ):
    import re
    filename_match = re.match( r'^((?:[A-Za-z0-9._-]+,)+)::(.*)$', input_file )
    if filename_match:
        generator_array = filename_match.group(1).split( ',' )
        genome.set_file_type( generator_array[0] )
        genome.add_generators( generator_array[1:-1] )
        genome.set_file_path( filename_match.group(2) )
    else:
        genome.set_file_path( input_file )
    #print( genome.identifier() )

def manage_input_thread( input_q, output_q ):
    while not input_q.empty():
        input_file = input_q.get()[0]
        new_genomes = []
        file_type = determine_file_type( input_file )
        if file_type == "frankenfasta":
            new_genomes = import_external_fasta( input_file )
        elif file_type == "vcf":
            new_genomes = read_vcf_file( input_file )
        for new_genome in new_genomes:
            output_q.put( [ new_genome ] )

def parse_input_files( input_files, num_threads, genomes ):
    from multiprocessing import Process, Queue
    input_q = Queue()
    output_q = Queue()
    for input_file in input_files:
        input_q.put( [ input_file ] )
    if num_threads > input_q.qsize():
        num_threads = input_q.qsize()
    thread_list = []
    for current_thread in range( num_threads ):
        current_thread = Process( target=manage_input_thread, args=[ input_q, output_q ] )
        thread_list.append( current_thread )
        current_thread.start()
    for current_thread in thread_list:
        current_thread.join()
    while not output_q.empty():
        genomes.add_genome( output_q.get()[0] )


def main():
    commandline_args = _parse_args()
    from nasp_objects import ReferenceGenome, GenomeCollection
    reference = ReferenceGenome()
    import_reference( reference, commandline_args.reference_fasta, commandline_args.reference_dups )
    genomes = GenomeCollection()
    genomes.set_reference( reference )
    parse_input_files( commandline_args.input_files, commandline_args.num_threads, genomes )
    #write_allcallable_matrix( genomes, commandline_args.master_matrix )

if __name__ == "__main__": main()


