#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.5"
__email__ = "dsmith@tgen.org"

import logging


class GenomeStatus:

    def __init__( self ):
        self._status_data = {}
        self._current_contig = None

    def set_current_contig( self, contig_name ):
        if contig_name not in self._status_data:
            raise InvalidContigName()
        else:
            self._current_contig = contig_name

    def add_contig( self, contig_name, change_current_contig = True ):
        if contig_name is not None:
            if change_current_contig:
                self._current_contig = contig_name
            if contig_name not in self._status_data:
                self._status_data[contig_name] = ""
        else:
            raise InvalidContigName()

    def append_contig( self, genome_data, contig_name = None, change_current_contig = True ):
        if contig_name is None:
            contig_name = self._current_contig
        self.add_contig( contig_name, change_current_contig )
        self._status_data[contig_name] = '' + self._status_data[contig_name] + genome_data

    def get_contigs( self ):
        return sorted( self._status_data.keys() )

    def extend_contig( self, new_length, missing_range_filler, contig_name = None, change_current_contig = True ):
        if contig_name is None:
            contig_name = self._current_contig
        self.add_contig( contig_name, change_current_contig )
        if len( self._status_data[contig_name] ) < new_length:
            self._status_data[contig_name] = '' + self._status_data[contig_name] + ( missing_range_filler * ( new_length - len( self._status_data[contig_name] ) ) )

    def set_value( self, new_data, first_position, missing_range_filler = "!", contig_name = None, change_current_contig = False ):
        if contig_name is None:
            contig_name = self._current_contig
        self.add_contig( contig_name, change_current_contig )
        first_position = first_position - 1
        self.extend_contig( first_position, missing_range_filler, contig_name )
        self._status_data[contig_name] = '' + self._status_data[contig_name][:first_position] + new_data + self._status_data[contig_name][( first_position + len( new_data ) ):]

    def get_value( self, first_position, last_position = None, contig_name = None, filler_value = None ):
        if contig_name is None:
            contig_name = self._current_contig
        queried_value = filler_value
        if contig_name in self.get_contigs():
            if last_position == -1:
                last_position = len( self._status_data[contig_name] )
            elif ( last_position is None ) or ( last_position < first_position ):
                last_position = first_position
            first_position = first_position - 1
            if first_position < len( self._status_data[contig_name] ):
                queried_value = self._status_data[contig_name][first_position:last_position]
        elif filler_value == None:
            raise InvalidContigName()
        return queried_value

    def get_contig_length( self, contig_name = None ):
        if contig_name is None:
            contig_name = self._current_contig
        return len( self._status_data[contig_name] )

    def _send_to_fasta_handle( self, output_handle, contig_prefix = "", max_chars_per_line = 80 ):
        for current_contig in self.get_contigs():
            output_handle.write( ">" + contig_prefix + current_contig + "\n" )
            if max_chars_per_line > 0:
                i = 0
                while ( max_chars_per_line * i ) < len( self._status_data[current_contig] ):
                    output_handle.write( '' + self._status_data[current_contig][( max_chars_per_line * i ):( max_chars_per_line * ( i + 1 ) )] + "\n" )
                    i = i + 1
            else:
                output_handle.write( '' + self._status_data[current_contig] + "\n" )

    def write_to_fasta_file( self, output_filename, contig_prefix = "", max_chars_per_line = 80 ):
        output_handle = open( output_filename, 'w' )
        self._send_to_fasta_handle( output_handle, contig_prefix, max_chars_per_line )
        output_handle.close()


class Genome:

    def __init__( self ):
        self._genome = GenomeStatus()

    def set_current_contig( self, contig_name ):
        self._genome.set_current_contig( contig_name )

    def append_contig( self, genome_data, contig_name = None ):
        self._genome.append_contig( genome_data, contig_name )

    def add_contig( self, contig_name, change_current_contig = True ):
        self._genome.add_contig( contig_name, change_current_contig )

    def get_contigs( self ):
        return self._genome.get_contigs()

    def extend_contig( self, new_length, missing_range_filler, contig_name = None, change_current_contig = False ):
        self._genome.extend_contig( new_length, missing_range_filler, contig_name, change_current_contig )

    def set_call( self, new_data, first_position, missing_range_filler = "X", contig_name = None, change_current_contig = False ):
        self._genome.set_value( new_data, first_position, missing_range_filler, contig_name, change_current_contig )

    def get_call( self, first_position, last_position = None, contig_name = None, filler_value = "X" ):
        return self._genome.get_value( first_position, last_position, contig_name, filler_value )

    def _import_fasta_line( self, line_from_fasta, contig_prefix = "" ):
        import re
        contig_match = re.match( r'^>' + re.escape( contig_prefix ) + r'([^\s]+)(?:\s|$)', line_from_fasta )
        if contig_match:
            self._genome.add_contig( contig_match.group(1) )
        else:
            data_match = re.match( r'^([A-Za-z.-]+)\s*$', line_from_fasta )
            if data_match:
                self._genome.append_contig( data_match.group(1) )

    def import_fasta_file( self, fasta_filename, contig_prefix = "" ):
        fasta_handle = open( fasta_filename, 'r' )
        for line_from_fasta in fasta_handle:
            self._import_fasta_line( line_from_fasta, contig_prefix )
        fasta_handle.close()

    def get_contig_length( self, contig_name = None ):
        return self._genome.get_contig_length( contig_name )

    def write_to_fasta_file( self, output_filename, contig_prefix = "", max_chars_per_line = 80 ):
        self._genome.write_to_fasta_file( output_filename, contig_prefix, max_chars_per_line )

    @staticmethod
    def reverse_complement( dna_string ):
        return dna_string.translate( ''.maketrans( 'ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr' ) )[::-1]

    @staticmethod
    def simple_call( dna_string, allow_x = False, allow_del = False ):
        simple_base = 'N'
        if len( dna_string ) > 0:
            simple_base = dna_string[0:1]
        simple_base = simple_base.upper()
        if simple_base == 'U':
            simple_base = 'T'
        if simple_base not in [ 'A', 'C', 'G', 'T', 'X', '.' ]:
            simple_base = 'N'
        if not allow_x and ( simple_base == 'X' ):
            simple_base = 'N'
        if not allow_del and ( simple_base == '.' ):
            simple_base = 'N'
        return simple_base


class GenomeMeta:

    def __init__( self ):
        self._nickname = None
        self._file_path = None
        self._file_type = None
        self._generators = [] # This should probably be a dictionary someday

    def set_file_path( self, file_path ):
        self._file_path = file_path
        if self._nickname == None:
            self._nickname = GenomeMeta.generate_nickname_from_filename( file_path )

    def set_file_type( self, type_string ):
        self._file_type = type_string

    def set_nickname( self, nickname ):
        self._nickname = nickname

    def add_generators( self, generator_array ):
        self._generators.extend( generator_array )

    def file_path( self ):
        return self._file_path

    def file_type( self ):
        return self._file_type

    def nickname( self ):
        return self._nickname

    def identifier( self ):
        identifier = self._nickname
        if len( self._generators ) > 0:
            identifier = '' + identifier + "::" + ( ','.join( self._generators ) )
        return identifier

    @staticmethod
    def generate_nickname_from_filename( filename ):
        import re
        import random
        filename_match = re.match( r'^(?:.*\/)?([^\/]+?)(?:\.(?:[Ff][Rr][Aa][Nn][Kk][Ee][Nn])?[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?|\.[Vv][Cc][Ff])?$', filename )
        if filename_match:
            nickname = filename_match.group(1)
        else:
            nickname = "file_" + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) )
        return nickname

    @staticmethod
    def reverse_complement( dna_string ):
        return dna_string.translate( ''.maketrans( 'ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr' ) )[::-1]


class IndelList:

    def __init__( self ):
        self._indels = {}


class ReferenceGenome( Genome ):

    def __init__( self ):
        Genome.__init__( self )
        self._dups = GenomeStatus()

    def get_dups_call( self, first_position, last_position = None, contig_name = None ):
        return self._dups.get_value( first_position, last_position, contig_name, "?" )

    def _import_dups_line( self, line_from_dups_file, contig_prefix = "" ):
        import re
        contig_match = re.match( r'^>' + re.escape( contig_prefix ) + r'([^\s]+)(?:\s|$)', line_from_dups_file )
        if contig_match:
            self.add_contig( contig_match.group(1), False )
            self._dups.add_contig( contig_match.group(1) )
        else:
            data_match = re.match( r'^([01-]+)\s*$', line_from_dups_file )
            if data_match:
                self._dups.append_contig( data_match.group(1) )

    def import_dups_file( self, dups_filename, contig_prefix = "" ):
        dups_handle = open( dups_filename, 'r' )
        for line_from_dups_file in dups_handle:
            self._import_dups_line( line_from_dups_file, contig_prefix )
        dups_handle.close()


class FastaGenome( Genome ):

    def __init__( self ):
        Genome.__init__( self )
        self._meta = GenomeMeta()
        self._indels = IndelList()

    def set_file_path( self, file_path ):
        self._meta.set_file_path( file_path )

    def set_file_type( self, type_string ):
        self._meta.set_file_type( type_string )

    def set_nickname( self, nickname ):
        self._meta.set_nickname( nickname )

    def add_generators( self, generator_name ):
        self._meta.add_generators( generator_name )

    def file_path( self ):
        return self._meta.file_path()

    def file_type( self ):
        return self._meta.file_type()

    def nickname( self ):
        return self._meta.nickname()

    def identifier( self ):
        return self._meta.identifier()

    # FIXME This data should be put into a real structure when the fasta code is pulled in
    def get_was_called( self, current_pos, contig_name = None ):
        return_value = "N"
        call_to_check = self.get_call( current_pos, None, contig_name, "X" )
        if call_to_check != "X" and call_to_check != "N":
            return_value = "Y"
        return return_value

    def get_coverage_pass( self, current_pos, contig_name = None ):
        return "-"

    def get_proportion_pass( self, current_pos, contig_name = None ):
        return "-"


class VCFGenome( Genome ):

    def __init__( self ):
        Genome.__init__( self )
        self._meta = GenomeMeta()
        self._indels = IndelList()
        self._was_called = GenomeStatus()
        self._passed_coverage = GenomeStatus()
        self._passed_proportion = GenomeStatus()

    def set_file_path( self, file_path ):
        self._meta.set_file_path( file_path )

    def set_file_type( self, type_string ):
        self._meta.set_file_type( type_string )

    def set_nickname( self, nickname ):
        self._meta.set_nickname( nickname )

    def add_generators( self, generator_name ):
        self._meta.add_generators( generator_name )

    def file_path( self ):
        return self._meta.file_path()

    def file_type( self ):
        return self._meta.file_type()

    def nickname( self ):
        return self._meta.nickname()

    def identifier( self ):
        return self._meta.identifier()

    def set_was_called( self, pass_value, current_pos, contig_name = None, change_current_contig = False ):
        self._was_called.set_value( pass_value, current_pos, "N", contig_name, change_current_contig )

    def set_coverage_pass( self, pass_value, current_pos, contig_name = None, change_current_contig = False ):
        self._passed_coverage.set_value( pass_value, current_pos, "?", contig_name, change_current_contig )

    def set_proportion_pass( self, pass_value, current_pos, contig_name = None, change_current_contig = False ):
        self._passed_proportion.set_value( pass_value, current_pos, "?", contig_name, change_current_contig )

    def get_was_called( self, current_pos, contig_name = None ):
        return self._was_called.get_value( current_pos, None, contig_name, "N" )

    def get_coverage_pass( self, current_pos, contig_name = None ):
        return self._passed_coverage.get_value( current_pos, None, contig_name, "?" )

    def get_proportion_pass( self, current_pos, contig_name = None ):
        return self._passed_proportion.get_value( current_pos, None, contig_name, "?" )


class CollectionStatistics:

    def __init__( self ):
        self._contig_stats = {}
        self._sample_stats = {}
        self._cumulative_cache = {}

    def _increment_by_contig( self, stat_id, contig_name ):
        if ( stat_id, contig_name ) not in self._contig_stats:
            self._contig_stats[( stat_id, contig_name )] = 0
        self._contig_stats[( stat_id, contig_name )] += 1

    def increment_contig_stat( self, stat_id, contig_name = None ):
        self._increment_by_contig( stat_id, contig_name )
        if contig_name is not None:
            self._increment_by_contig( stat_id, None )

    def get_contig_stat( self, stat_id, contig_name = None ):
        return_value = 0 
        if ( stat_id, contig_name ) in self._contig_stats:
            return_value = self._contig_stats[( stat_id, contig_name )]
        return return_value

    def _cache_cumulative_stats( self, stat_id, sample_nickname, did_pass ):
        if ( stat_id, sample_nickname, 't' ) not in self._cumulative_cache:
            self._cumulative_cache[( stat_id, sample_nickname, 't' )] = 0
            self._cumulative_cache[( stat_id, sample_nickname, 'p' )] = 0
        self._cumulative_cache[( stat_id, sample_nickname, 't' )] += 1
        if did_pass:
            self._cumulative_cache[( stat_id, sample_nickname, 'p' )] += 1

    def _increment_by_sample( self, stat_id, sample_nickname, sample_info, cum_type ):
        if ( stat_id, sample_nickname, sample_info, cum_type ) not in self._sample_stats:
            self._sample_stats[( stat_id, sample_nickname, sample_info, cum_type )] = 0
        self._sample_stats[( stat_id, sample_nickname, sample_info, cum_type )] += 1

    def record_sample_stat( self, stat_id, sample_nickname, sample_identifier, sample_path, did_pass ):
        if did_pass:
            self._increment_by_sample( stat_id, sample_nickname, ( sample_identifier, sample_path ), None )
        self._cache_cumulative_stats( stat_id, sample_nickname, did_pass )
        self._cache_cumulative_stats( stat_id, None, did_pass )

    def get_sample_stat( self, stat_id, sample_nickname, sample_identifier, sample_path ):
        return_value = 0 
        if ( stat_id, sample_nickname, ( sample_identifier, sample_path ), None ) in self._contig_stats:
            return_value = self._contig_stats[( stat_id, sample_nickname, ( sample_identifier, sample_path ), None )]
        return return_value

    def get_cumulative_stat( self, stat_id, cum_type, sample_nickname = None ):
        return_value = 0 
        if ( stat_id, sample_nickname, None, cum_type ) in self._contig_stats:
            return_value = self._contig_stats[( stat_id, sample_nickname, None, cum_type )]
        return return_value

    def flush_cumulative_stat_cache( self ):
        for ( stat_id, sample_nickname, stat_type ) in self._cumulative_cache:
            if stat_type == 'p':
                if self._cumulative_cache[( stat_id, sample_nickname, 'p' )] == self._cumulative_cache[( stat_id, sample_nickname, 't' )]:
                    self._increment_by_sample( stat_id, sample_nickname, None, 'all' )
                if self._cumulative_cache[( stat_id, sample_nickname, 'p' )] > 0:
                    self._increment_by_sample( stat_id, sample_nickname, None, 'any' )


class GenomeCollection:

    def __init__( self ):
        self._reference = None
        self._genomes = []
        self._stats = CollectionStatistics()

    # To preserve sample-analysis order across different runs
    @staticmethod
    def _get_key( genome ):
        return genome.identifier()

    def set_reference( self, reference ):
        self._reference = reference

    def reference( self ):
        return self._reference

    def get_dups_call( self, first_position, last_position = None, contig_name = None ):
        return self._reference.get_dups_call( first_position, last_position, contig_name )

    def add_genome( self, genome ):
        self._genomes.append( genome )
        self._genomes.sort( key=GenomeCollection._get_key )

    def set_current_contig( self, contig_name ):
        self._reference.set_current_contig( contig_name )
        for genome in self._genomes:
            genome.set_current_contig( contig_name )

    def get_contigs( self ):
        return self._reference.get_contigs()

    def get_genome_count( self ):
        return len( self._genomes )

    def increment_contig_stat( self, stat_id, contig_name = None ):
        self._stats.increment_contig_stat( stat_id, contig_name )

    def get_contig_stat( self, stat_id, contig_name = None ):
        return self._stats.get_contig_stat( stat_id, contig_name )

    def record_sample_stat( self, stat_id, sample_nickname, sample_identifier, sample_path, did_pass ):
        self._stats.record_sample_stat( stat_id, sample_nickname, sample_identifier, sample_path, did_pass )

    def get_sample_stat( self, stat_id, sample_nickname, sample_identifier, sample_path ):
        return self._stats.get_sample_stat( stat_id, sample_nickname, sample_identifier, sample_path )

    def get_cumulative_stat( self, stat_id, cum_type, sample_nickname = None ):
        return self._stats.get_cumulative_stat( stat_id, cum_type, sample_nickname )

    def flush_cumulative_stat_cache( self ):
        self._stats.flush_cumulative_stat_cache()

    # FIXME split into a larger number of smaller more testable functions
    # FIXME Some of this doesn't belong here
    def _format_matrix_line( self, current_contig, current_pos, matrix_format ):
        genome_count = self.get_genome_count()
        matrix_line = '' + current_contig + "::" + str( current_pos ) + "\t"
        reference_call = self._reference.get_call( current_pos, None, current_contig )
        simplified_refcall = Genome.simple_call( reference_call )
        self.increment_contig_stat( 'reference_length', current_contig )
        if simplified_refcall != 'N':
            self.increment_contig_stat( 'reference_clean', current_contig )
        matrix_line += '' + reference_call + "\t"
        custom_line = '' + matrix_line
        call_data = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'indel': 0, 'snpcall': 0, 'indelcall': 0, 'refcall': 0, 'callstring': '', 'covstring': '', 'propstring': '', 'called': 0, 'passcov': 0, 'passprop': 0 }
        consensus_check = {}
        # The expensive loop, single threaded and runs for every sample-analysis-contig-position
        for genome in self._genomes:
            sample_call = genome.get_call( current_pos, None, current_contig, 'X' )
            simplified_sample_call = Genome.simple_call( sample_call )
            matrix_line += '' + sample_call + "\t"
            call_data[simplified_sample_call] += 1
            genome_nickname = genome.nickname()
            genome_identifier = genome.identifier()
            genome_path = genome.file_path()
            was_called = genome.get_was_called( current_pos, current_contig )
            call_data['callstring'] += '' + was_called
            if was_called == 'Y':
                was_called = True
                call_data['called'] += 1
            else:
                was_called = False
            self.record_sample_stat( 'was_called', genome_nickname, genome_identifier, genome_path, was_called )
            passed_coverage = genome.get_coverage_pass( current_pos, current_contig )
            call_data['covstring'] += '' + passed_coverage
            if passed_coverage == 'Y' or passed_coverage == '-':
                passed_coverage = True
                call_data['passcov'] += 1
            else:
                passed_coverage = False
            self.record_sample_stat( 'passed_coverage_filter', genome_nickname, genome_identifier, genome_path, passed_coverage )
            passed_proportion = genome.get_proportion_pass( current_pos, current_contig )
            call_data['propstring'] += '' + passed_proportion
            if passed_proportion == 'Y' or passed_proportion == '-':
                passed_proportion = True
                call_data['passprop'] += 1
            else:
                passed_proportion = False
            self.record_sample_stat( 'passed_proportion_filter', genome_nickname, genome_identifier, genome_path, passed_proportion )
            if was_called and passed_coverage and passed_proportion:
                if genome_nickname in consensus_check:
                    if consensus_check[genome_nickname] != simplified_sample_call:
                        consensus_check[genome_nickname] = 'N'
                else:
                    consensus_check[genome_nickname] = simplified_sample_call
            else:
                consensus_check[genome_nickname] = 'N'
            if was_called and passed_coverage and passed_proportion and simplified_refcall != 'N':
                if simplified_refcall == simplified_sample_call:
                    call_data['refcall'] += 1
                elif simplified_sample_call != 'N':
                    call_data['snpcall'] += 1
            if matrix_format is None:
                custom_line += '' + sample_call + "\t"
            elif matrix_format == "missingdata":
                if was_called and passed_coverage and passed_proportion and simplified_sample_call != 'N':
                    custom_line += '' + sample_call + "\t"
                elif not was_called:
                    custom_line += "X\t"
                else:
                    custom_line += "N\t"
        for genome_nickname in consensus_check:
            if consensus_check[genome_nickname] != 'N':
                self.record_sample_stat( 'consensus', genome_nickname, None, None, True )
            else:
                self.record_sample_stat( 'consensus', genome_nickname, None, None, False )
        # FIXME The hard way? Why?
        matrix_line += '' + str( call_data['snpcall'] ) + "\t" + str( call_data['indelcall'] ) + "\t" + str( call_data['refcall'] ) + "\t"
        custom_line += '' + str( call_data['snpcall'] ) + "\t" + str( call_data['indelcall'] ) + "\t" + str( call_data['refcall'] ) + "\t"
        matrix_line += '' + str( call_data['called'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passcov'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passprop'] ) + "/" + str( genome_count ) + "\t"
        custom_line += '' + str( call_data['called'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passcov'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passprop'] ) + "/" + str( genome_count ) + "\t"
        matrix_line += '' + str( call_data['A'] ) + "\t" + str( call_data['C'] ) + "\t" + str( call_data['G'] ) + "\t" + str( call_data['T'] ) + "\t" + str( call_data['indel'] ) + "\t" + str( call_data['N'] ) + "\t"
        custom_line += '' + str( call_data['A'] ) + "\t" + str( call_data['C'] ) + "\t" + str( call_data['G'] ) + "\t" + str( call_data['T'] ) + "\t" + str( call_data['indel'] ) + "\t" + str( call_data['N'] ) + "\t"
        matrix_line += '' + current_contig + "\t" + str( current_pos ) + "\t"
        custom_line += '' + current_contig + "\t" + str( current_pos ) + "\t"
        dups_call = self._reference.get_dups_call( current_pos, None, current_contig )
        if dups_call == "1":
            dups_call = True
            self.increment_contig_stat( 'reference_duplicated', current_contig )
        else:
            dups_call = False
        if 'N' not in consensus_check.values():
            consensus_check = True
            self.increment_contig_stat( 'all_passed_consensus', current_contig )
        else:
            consensus_check = False
        if call_data['called'] == genome_count:
            self.increment_contig_stat( 'all_called', current_contig )
        if call_data['passcov'] == genome_count:
            self.increment_contig_stat( 'all_passed_coverage', current_contig )
        if call_data['passprop'] == genome_count:
            self.increment_contig_stat( 'all_passed_proportion', current_contig )
        if consensus_check and not dups_call and call_data['called'] == genome_count and call_data['passcov'] == genome_count and call_data['passprop'] == genome_count and call_data['N'] == 0:
            self.increment_contig_stat( 'quality_breadth', current_contig )
            if call_data['snpcall'] > 0:
                self.increment_contig_stat( 'best_snps', current_contig )
        if not dups_call and call_data['snpcall'] > 0:
            self.increment_contig_stat( 'any_snps', current_contig )
        matrix_line += '' + str( dups_call ) + "\t" + str( consensus_check ) + "\t"
        custom_line += '' + str( dups_call ) + "\t" + str( consensus_check ) + "\t"
        matrix_line += '' + str( call_data['callstring'] ) + "\t" + str( call_data['covstring'] ) + "\t" + str( call_data['propstring'] ) + "\n"
        if matrix_format is None:
            if call_data['snpcall'] == 0 or call_data['indelcall'] > 0 or call_data['snpcall'] + call_data['refcall'] < genome_count or dups_call or not consensus_check:
                custom_line = None
            else:
                custom_line += "\n"
        elif matrix_format == "missingdata":
            if call_data['snpcall'] == 0 or call_data['indelcall'] > 0 or dups_call:
                custom_line = None
            else:
                custom_line += '' + str( call_data['callstring'] ) + "\t" + str( call_data['covstring'] ) + "\t" + str( call_data['propstring'] ) + "\n"
        self.flush_cumulative_stat_cache()
        return ( matrix_line, custom_line )

    def _send_to_matrix_handles( self, master_handle, custom_handle, matrix_format ):
        master_handle.write( "LocusID\tReference\t" )
        custom_handle.write( "LocusID\tReference\t" )
        for genome in self._genomes:
            master_handle.write( '' + genome.identifier() + "\t" )
            custom_handle.write( '' + genome.identifier() + "\t" )
        master_handle.write( "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\n" )
        if matrix_format is None:
            custom_handle.write( "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\n" )
        elif matrix_format == "missingdata":
            custom_handle.write( "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\n" )
        for current_contig in self.get_contigs():
            for current_pos in range( 1, self._reference.get_contig_length( current_contig ) + 1 ):
                matrix_lines = self._format_matrix_line( current_contig, current_pos, matrix_format )
                master_handle.write( matrix_lines[0] )
                if matrix_lines[1] is not None:
                    custom_handle.write( matrix_lines[1] )

    def write_to_matrices( self, master_filename, custom_filename, matrix_format ):
        master_handle = open( master_filename, 'w' )
        custom_handle = open( custom_filename, 'w' )
        self._send_to_matrix_handles( master_handle, custom_handle, matrix_format )
        master_handle.close()
        custom_handle.close()

    def _write_general_stats( self, general_handle ):
        general_stat_array = [ 'reference_length', 'reference_clean', 'reference_duplicated', 'all_called', 'all_passed_coverage', 'all_passed_proportion', 'all_passed_consensus', 'quality_breadth', 'any_snps', 'best_snps' ]
        denominator_stat = 'reference_length'
        general_handle.write( "Contig\t" )
        for current_stat in general_stat_array:
            general_handle.write( '' + current_stat + "\t" )
            if current_stat != denominator_stat:
                general_handle.write( '' + current_stat + " (%)\t" )
        general_handle.write( "\n" )
        general_handle.write( "\tstat descriptions go here\n" )
        for current_contig in ( [ None ] + self.get_contigs() ):
            denominator_value = self.get_contig_stat( denominator_stat, current_contig )
            if current_contig is None:
                general_handle.write( "Whole Genome\t" )
            else:
                general_handle.write( '' + current_contig + "\t" )
            for current_stat in general_stat_array:
                general_handle.write( '' + str( self.get_contig_stat( current_stat, current_contig ) ) + "\t" )
                if current_stat != denominator_stat:
                    general_handle.write( "%.2f%%\t" % ( self.get_contig_stat( current_stat, current_contig ) / denominator_value * 100 ) )
            general_handle.write( "\n" )

    def _write_sample_stats( self, sample_handle ):
        pass

    def write_to_stats_files( self, general_filename, sample_filename ):
        general_handle = open( general_filename, 'w' )
        sample_handle = open( sample_filename, 'w' )
        self._write_general_stats( general_handle )
        self._write_sample_stats( sample_handle )
        general_handle.close()
        sample_handle.close()
        #print( self._stats._contig_stats )
        #print( self._stats._sample_stats )


class VCFRecord:

    def __init__( self, file_path ):
        self._file_path = file_path
        self._file_handle = open( self._file_path, 'r' )
        self._header_list = []
        self._sample_list = []
        self._current_record = {}
        self._get_header_map()

    def _get_header_map( self ):
        current_line = ''
        while current_line[0:6] != "#CHROM":
            current_line = self._file_handle.readline()
        header_list = current_line.rstrip()[1:].split( "\t" )
        sample_headers_started = False
        for current_header in header_list:
            self._header_list.append( current_header )
            if sample_headers_started:
                self._sample_list.append( current_header )
            elif current_header == "FORMAT":
                sample_headers_started = True
        if sample_headers_started == False:
            self._sample_list.append( 'vcf_sample' )
        #print( self._header_list )
        #print( self._sample_list )

    def _get_record_alts( self ):
        self._current_record['alts'] = [ self._current_record['global']['REF'] ]
        if self._current_record['global']['ALT'] != '.':
            self._current_record['alts'] += self._current_record['global']['ALT'].split( ',' )

    def _get_record_info( self ):
        self._current_record['info'] = {}
        info_strings = self._current_record['global']['INFO'].split( ';' )
        for info_string in info_strings:
            info_keyval = info_string.split( '=', 1 )
            if len( info_keyval ) > 1:
                self._current_record['info'][info_keyval[0]] = info_keyval[1]
            else:
                self._current_record['info'][info_keyval[0]] = None

    def _get_record_samples( self ):
        if 'FORMAT' in self._current_record['global']:
            self._current_record['samples'] = {}
            sample_header = self._current_record['global']['FORMAT'].split( ':' )
            for current_sample in self._sample_list:
                self._current_record['samples'][current_sample] = dict( zip( sample_header, self._current_record['global'][current_sample].split( ':' ) ) )

    def fetch_next_record( self ):
        current_line = "#"
        while current_line[0:1] == "#":
            current_line = self._file_handle.readline()
        return_value = False
        if current_line != '':
            record_list = current_line.rstrip().split( "\t" )
            self._current_record['global'] = dict( zip( self._header_list, record_list ) )
            self._get_record_alts()
            self._get_record_info()
            self._get_record_samples()
            return_value = True
        #print( self._current_record )
        return return_value

    def get_samples( self ):
        return self._sample_list

    def get_contig( self ):
        return self._current_record['global']['CHROM']

    def get_position( self ):
        return int( self._current_record['global']['POS'] )

    def get_reference_call( self ):
        return self._current_record['global']['REF']

    def get_sample_call( self, current_sample ):
        return_value = None
        if len( self._current_record['alts'] ) == 1:
            return_value = self._current_record['alts'][0]
        elif 'FORMAT' in self._current_record['global']:
            if 'GT' in self._current_record['samples'][current_sample]:
                alt_number = self._current_record['samples'][current_sample]['GT'].split( '/', 1 )[0].split( '|', 1 )[0]
                if alt_number.isdigit():
                    return_value = self._current_record['alts'][int( alt_number )]
                    # OMG varscan
                    if len( self._current_record['global']['REF'] ) > 1 and ( len( self._current_record['global']['REF'] ) - 1 ) == len( return_value ) and self._current_record['global']['REF'][:len( return_value )] != return_value and self._current_record['global']['REF'][-len( return_value ):] == return_value:
                        return_value = self._current_record['alts'][0]
        elif len( self._current_record['alts'] ) > 1:
            return_value = self._current_record['alts'][1]
        return return_value

    def get_coverage( self, current_sample ):
        sample_coverage = None
        if 'FORMAT' in self._current_record['global'] and 'DP' in self._current_record['samples'][current_sample] and self._current_record['samples'][current_sample]['DP'] is not None and self._current_record['samples'][current_sample]['DP'].isdigit():
            sample_coverage = int( self._current_record['samples'][current_sample]['DP'] )
        elif 'DP' in self._current_record['info'] and self._current_record['info']['DP'] is not None and self._current_record['info']['DP'].isdigit():
            sample_coverage = int( self._current_record['info']['DP'] ) / len( self._sample_list )
        elif 'ADP' in self._current_record['info'] and self._current_record['info']['ADP'] is not None and self._current_record['info']['ADP'].isdigit():
            sample_coverage = int( self._current_record['info']['ADP'] ) / len( self._sample_list )
        return sample_coverage
    
    def get_proportion( self, current_sample, sample_coverage, is_a_snp ):
        sample_proportion = None
        if 'FORMAT' in self._current_record['global'] and 'AD' in self._current_record['samples'][current_sample]:
            call_depths = self._current_record['samples'][current_sample]['AD'].split( ',' )
        # gatk, reliable and documented
            if len( call_depths ) > 1:
                alt_number = self._current_record['samples'][current_sample]['GT'].split( '/', 1 )[0].split( '|', 1 )[0]
                if alt_number.isdigit():
                   sample_proportion = int( call_depths[int( alt_number )] ) / sample_coverage
        # varscan, reliable and documented
            elif is_a_snp:
                sample_proportion = int( call_depths[0] ) / sample_coverage
            elif not is_a_snp and 'RD' in self._current_record['samples'][current_sample]:
                sample_proportion = int( self._current_record['samples'][current_sample]['RD'] ) / sample_coverage
        # solsnp, undocumented, no multi-sample support
        elif 'AR' in self._current_record['info']:
            sample_proportion = float( self._current_record['info']['AR'] )
            if not is_a_snp:
                sample_proportion = 1 - sample_proportion
        # samtools, estimate, dubious accuracy
        elif 'DP4' in self._current_record['info']:
            call_depths = self._current_record['info']['DP4'].split( ',' )
            if is_a_snp:
                sample_proportion = ( int( call_depths[2] ) + int( call_depths[3] ) ) / ( sample_coverage * len( self._sample_list ) )
            else:
                sample_proportion = ( int( call_depths[0] ) + int( call_depths[1] ) ) / ( sample_coverage * len( self._sample_list ) )
        return sample_proportion
    
    def get_sample_info( self, current_sample ):
        sample_info = {}
        sample_info['call'] = self.get_sample_call( current_sample )
        sample_info['was_called'] = False
        sample_info['is_a_snp'] = False
        if sample_info['call'] is not None and sample_info['call'][0] != 'N':
            sample_info['was_called'] = True
            if Genome.simple_call( sample_info['call'][0] ) != Genome.simple_call( self._current_record['global']['REF'] ):
                sample_info['is_a_snp'] = True
        sample_info['is_a_snp']
        # FIXME indels
        sample_info['is_an_insert'] = None
        sample_info['is_a_delete'] = None
        sample_info['coverage'] = self.get_coverage( current_sample )
        sample_info['proportion'] = self.get_proportion( current_sample, sample_info['coverage'], sample_info['is_a_snp'] )
        return sample_info


# FIXME user feedback
class InvalidContigName( Exception ):
    pass


# FIXME user feedback
class ReferenceCallMismatch( Exception ):
    pass


def main():
    print( "This is just a class definition, to be included from other python programs." )
    print( "  ...You sillyface." )

if __name__ == "__main__": main()


