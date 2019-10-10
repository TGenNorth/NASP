#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.9"
__email__ = "dsmith@tgen.org"

import logging


def _parse_args():
    """
    Parse command line arguments and return an object whose attributes are the
    settings provided, as per standard argparse behavior.
    For backward compatibility when we were changing dispatchers, the options
    could be directly provided on the command line, or a single argument could
    specify the path to the XML file with the remaining arguments in it.  As
    the old dispatcher is deprecated, so perhaps is every option below except
    '--dto-file', which would then be made required.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--mode", required=True, choices=['commandline', 'xml'],
                        help="Data passing mode, must be set to 'commandline' or 'xml'.")
    parser.add_argument("--reference-fasta", help="Path to input reference fasta file.")
    parser.add_argument("--reference-dups", help="Path to input reference dups file.")
    parser.add_argument("--input-files", nargs="+", help="Path to input VCF/fasta files for matrix conversion.")
    parser.add_argument("--matrix-folder", default="matrices", help="Name of folder to write output matries to.")
    parser.add_argument("--stats-folder", default="statistics", help="Name of folder to write statistics files to.")
    parser.add_argument("--minimum-coverage", type=int, default=10, help="Minimum coverage depth at a position.")
    parser.add_argument("--minimum-proportion", type=float, default=0.9,
                        help="Minimum proportion of reads that must match the call at a position.")
    parser.add_argument("--num-threads", type=int, default=1, help="Number of threads to use when processing input.")
    parser.add_argument("--dto-file", help="Path to a matrix_dto XML file that defines all the parameters.")
    return parser.parse_args()


def _parse_input_config(commandline_args):
    """
    Set the configuration options from the XML file and add them as attributes
    to the object returned by _parse_args().
    Since getting these options directly from the command line is perhaps
    deprecated, we should probably store the configuration options from the
    XML file in a separate dictionary, instead of mixing the two data sources
    together as now.
    """
    import nasp.matrix_DTO as matrix_DTO

    (matrix_parms, input_files) = matrix_DTO.parse_dto(commandline_args.dto_file)
    commandline_args.reference_fasta = matrix_parms['reference-fasta']
    commandline_args.reference_dups = matrix_parms.get('reference-dups')
    commandline_args.matrix_folder = matrix_parms['matrix-folder']
    commandline_args.stats_folder = matrix_parms['stats-folder']
    commandline_args.minimum_coverage = int(
        matrix_parms['minimum-coverage']) if "minimum-coverage" in matrix_parms else 0
    commandline_args.minimum_proportion = float(
        matrix_parms['minimum-proportion']) if "minimum-proportion" in matrix_parms else 0
    if "filter-matrix-format" in matrix_parms:
        commandline_args.filter_matrix_format = matrix_parms['filter-matrix-format']
    else:
        commandline_args.filter_matrix_format = None
    commandline_args.input_files = input_files
    return commandline_args


def import_reference(reference, reference_path, dups_path):
    """
    Take an empty reference object and populate it with the data from a
    reference file, and a dups file if any.
    Does not return anything, as the passed-in object is modified.
    """
    reference.import_fasta_file(reference_path)
    if dups_path is not None:
        reference.import_dups_file(dups_path)
        # from sys import stdout
        #reference._genome._send_to_fasta_handle( stdout )
        #reference._dups._send_to_fasta_handle( stdout )


def import_external_fasta(input_file):
    """
    Create a FastaGenome object, set its metadata, and populate it with the
    data from a fasta file.
    Must return the data as the single item in an array, because other file
    formats potentially contain several genomes.
    """
    from nasp.nasp_objects import FastaGenome

    genome = FastaGenome()
    set_genome_metadata(genome, input_file)
    genome.import_fasta_file(genome.file_path(), "franken::")
    # from sys import stdout
    #genome._genome._send_to_fasta_handle( stdout )
    return [genome]


# FIXME split into a larger number of smaller more testable functions
# FIXME This belongs in VCFGenome object perhaps?
def read_vcf_file(reference, min_coverage, min_proportion, input_file):
    """
    Submit VCF to be read in to VCF parser, populate genome data and filter
    data from the parsed VCF data, return a list of the read-in genomes.
    """
    genomes = {}
    file_path = get_file_path(input_file)
    with open(file_path, 'r') as vcf_filehandle:
        from nasp.nasp_objects import VCFGenome, Genome, ReferenceCallMismatch, VCFRecord

        vcf_record = VCFRecord(file_path)
        vcf_samples = vcf_record.get_samples()
        for vcf_sample in vcf_samples:
            genomes[vcf_sample] = VCFGenome()
            set_genome_metadata(genomes[vcf_sample], input_file)
            genomes[vcf_sample].set_nickname(vcf_sample)
        while vcf_record.fetch_next_record():
            current_contig = vcf_record.get_contig()
            current_pos = vcf_record.get_position()
            # Skip if position isn't in reference; maybe user truncated reference to exclude an uninteresting region.
            if current_pos <= reference.get_contig_length(current_contig):
                reference_call = reference.get_call(current_pos, None, current_contig)
                simplified_refcall = Genome.simple_call(reference_call)
                if ( simplified_refcall != 'N' ) and (
                    simplified_refcall != Genome.simple_call(vcf_record.get_reference_call()[0]) ):
                    # Reference call from reference fasta differs from reference call in VCF file at the same position.
                    raise ReferenceCallMismatch(reference_call, vcf_record.get_reference_call(), file_path,
                                                current_contig, current_pos)
                for vcf_sample in vcf_samples:
                    sample_info = vcf_record.get_sample_info(vcf_sample)
                    # FIXME indels
                    if sample_info['call'] is not None:
                        genomes[vcf_sample].set_call(sample_info['call'], current_pos, 'X', current_contig)
                    if sample_info['was_called']:
                        genomes[vcf_sample].set_was_called('Y', current_pos, current_contig)
                    if sample_info['coverage'] is not None:
                        if sample_info['coverage'] == 'PASS' or sample_info['coverage'] >= min_coverage:
                            genomes[vcf_sample].set_coverage_pass('Y', current_pos, current_contig)
                        else:
                            genomes[vcf_sample].set_coverage_pass('N', current_pos, current_contig)
                    if sample_info['proportion'] is not None:
                        if sample_info['proportion'] == 'PASS' or sample_info['proportion'] >= min_proportion:
                            genomes[vcf_sample].set_proportion_pass('Y', current_pos, current_contig)
                        else:
                            genomes[vcf_sample].set_proportion_pass('N', current_pos, current_contig)
                    elif not sample_info['is_a_snp']:
                        # Some big SNP callers, like GATK, do not provide proportion information when
                        # the position is called reference.  We cannot filter these positions.
                        genomes[vcf_sample].set_proportion_pass('-', current_pos, current_contig)
    # from sys import stdout
    #for genome in genomes:
    #    genomes[genome]._genome._send_to_fasta_handle( stdout )
    return genomes.values()


# FIXME The information these three functions capture is in the XML file.
# FIXME These should all be deprecated and the XML data used instead.
def determine_file_type(input_file):
    """ Get the file type from the packed input filename string """
    import re

    file_type = None
    # Format is "type,aligner,caller,::path"
    filename_match = re.match(r'^((?:[^,:]+,)+)::.*$', input_file)
    if filename_match:
        generator_array = filename_match.group(1).split(',')
        file_type = generator_array[0]
    # print( file_type )
    return file_type


def get_file_path(input_file):
    """ Get the file path from the packed input filename string """
    import re

    file_path = None
    # Format is "type,aligner,caller,::path"
    filename_match = re.match(r'^(?:[^,:]+,)+::(.*)$', input_file)
    if filename_match:
        file_path = filename_match.group(1)
    else:
        file_path = input_file
    # print( file_path )
    return file_path


def set_genome_metadata(genome, input_file):
    """ Get all the metadata from the packed input filename string """
    import re
    # Format is "type,aligner,caller,::path"
    filename_match = re.match(r'^((?:[^,:]+,)+)::(.*)$', input_file)
    if filename_match:
        generator_array = filename_match.group(1).split(',')
        genome.set_file_type(generator_array[0])
        genome.add_generators(generator_array[1:-1])
        genome.set_file_path(filename_match.group(2))
    else:
        genome.set_file_path(input_file)
        # print( genome.identifier() )


def manage_input_thread(reference, min_coverage, min_proportion, input_q, output_q):
    """
    Manage one input file worker thread, for reading the data from the file.
    Input filenames are pulled one at a time from the input queue, and the
    genome data from the read-in files is placed on the output queue.
    When an input filename of "None" appears, we know we're done and put
    "None" on the output queue so the controlling thread knows we won't be
    adding more data.
    """
    input_file = input_q.get()
    while input_file is not None:
        try:
            new_genomes = []
            file_type = determine_file_type(input_file)
            if file_type == "frankenfasta":
                new_genomes = import_external_fasta(input_file)
            elif file_type == "vcf":
                new_genomes = read_vcf_file(reference, min_coverage, min_proportion, input_file)
            for new_genome in new_genomes:
                output_q.put(new_genome)
        except:
            failed_file_path = get_file_path(input_file)
            logging.exception("Unable to read in data from '{0}'!".format(failed_file_path))
            # If the genome on the output queue is just a string, its contents
            # are the name of the file we failed to read in.
            output_q.put(failed_file_path)
        input_file = input_q.get()
    output_q.put(None)


def parse_input_files(input_files, num_threads, genomes, min_coverage, min_proportion):
    """
    Use a pool of worker threads to, in parallel, read in the input files.
    Populate the genome collection with the read-in data.
    This is the "poison pill" thread management algorithm, where threads are
    each given a "you can stop now" task once the actual queue of tasks is
    complete.
    """
    # Lines below marked "Single-thread" can be uncommented and replace
    # the lines marked "Multi-thread" to get single-thread behavior.
    from multiprocessing import Process, Queue  # Multi-thread
    # from queue import Queue  # Single-thread
    from time import sleep

    input_q = Queue()
    output_q = Queue()
    for input_file in input_files:
        input_q.put(input_file)
    # If the number of jobs is already smaller than the thread pool...
    if num_threads > input_q.qsize():
        num_threads = input_q.qsize()
    sleep(1)
    thread_list = []
    for current_thread in range(num_threads):
        input_q.put(None)
        current_thread = Process(target=manage_input_thread,
                                 args=[genomes.reference(), min_coverage, min_proportion, input_q,
                                       output_q])  # Multi-thread
        current_thread.start()  # Multi-thread
        #manage_input_thread( genomes.reference(), min_coverage, min_proportion, input_q, output_q )  # Single-thread
        thread_list.append(current_thread)  # Multi-thread
    sleep(1)
    while num_threads > 0:
        new_genome = output_q.get()
        if new_genome is None:
            num_threads -= 1
        elif isinstance(new_genome, str):
            # Reading this file in failed.  We only know the filename.
            genomes.add_failed_genome(new_genome)
        else:
            genomes.add_genome(new_genome)
    sleep(1)
    for current_thread in thread_list:  # Multi-thread
        current_thread.join()  # Multi-thread


def write_output_matrices(genomes, matrix_folder, matrix_format_choices):
    """
    Write matrices from genome collection data.
    Defines the matrix types for future expansion of custom matrix options.
    This information eventually should come from the user interface and be
    included in the XML configuration file, rather than hardcoded here.
    The matrix_format_choices option comes from the XML to here.
    """
    matrix_formats = [
        {
            'filename': '' + matrix_folder + '/master_matrix.tsv',
            'dataformat': 'matrix',
            'filter': 'allcallable'
        },
        {
            'filename': '' + matrix_folder + '/bestsnp_matrix.tsv',
            'dataformat': 'matrix',
            'filter': 'bestsnp'
        },
        {
            'filename': '' + matrix_folder + '/bestsnp_matrix.snpfasta',
            'dataformat': 'fasta',
            'filter': 'bestsnp'
        },
        {
            'filename': '' + matrix_folder + '/bestsnp_matrix.vcf',
            'dataformat': 'vcf',
            'filter': 'bestsnp'
        },
        {
            'filename': '' + matrix_folder + '/missingdata_matrix.tsv',
            'dataformat': 'matrix',
            'filter': 'missingdata'
        },
        {
            'filename': '' + matrix_folder + '/missingdata_matrix.snpfasta',
            'dataformat': 'fasta',
            'filter': 'missingdata'
        },
        {
            'filename': '' + matrix_folder + '/missingdata_matrix.vcf',
            'dataformat': 'vcf',
            'filter': 'missingdata'
        }]
    # The XML format transferring this information should be more elaborate someday.
    if matrix_format_choices is not None and matrix_format_choices == "include_allref_pos":
        matrix_formats.extend( [
            {
                'filename': '' + matrix_folder + '/withallrefpos_matrix.tsv',
                'dataformat': 'matrix',
                'filter': 'includeref'
            },
            {
                'filename': '' + matrix_folder + '/withallrefpos_matrix.snpfasta',
                'dataformat': 'fasta',
                'filter': 'includeref'
            },
            {
                'filename': '' + matrix_folder + '/withallrefpos_matrix.vcf',
                'dataformat': 'vcf',
                'filter': 'includeref'
            }] )
    genomes.write_to_matrices(matrix_formats)


def write_stats_data(genomes, stats_folder):
    """ Write stats data from genome collection to preset filenames """
    general_stats_file = '' + stats_folder + '/general_stats.tsv'
    sample_stats_file = '' + stats_folder + '/sample_stats.tsv'
    genomes.write_to_stats_files(general_stats_file, sample_stats_file)


def main():
    """
    Main flow control function for the script.  Is sequential over the
    following steps:
    1.  Parse arguments
    2.  Read in reference
    3.  Read in all query genomes in parallel
    4.  Write output matrices
    5.  Write stats files
    """
    commandline_args = _parse_args()
    if commandline_args.dto_file:
        commandline_args = _parse_input_config(commandline_args)
    logging.basicConfig(level=logging.WARNING)
    from nasp.nasp_objects import ReferenceGenome, GenomeCollection

    reference = ReferenceGenome()
    import_reference(reference, commandline_args.reference_fasta, commandline_args.reference_dups)
    genomes = GenomeCollection()
    genomes.set_reference(reference)
    parse_input_files(commandline_args.input_files, commandline_args.num_threads, genomes,
                      commandline_args.minimum_coverage, commandline_args.minimum_proportion)
    write_output_matrices(genomes, commandline_args.matrix_folder, commandline_args.filter_matrix_format)
    write_stats_data(genomes, commandline_args.stats_folder)


if __name__ == "__main__":
    main()


