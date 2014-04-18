#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.2"
__email__ = "dsmith@tgen.org"


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

    def get_coverage_pass( self, current_pos, contig_name = None ):
        return "-"

    def get_proportion_pass( self, current_pos, contig_name = None ):
        return "-"


class VCFGenome( Genome ):

    def __init__( self ):
        Genome.__init__( self )
        self._meta = GenomeMeta()
        self._indels = IndelList()
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

    def set_coverage_pass( self, pass_value, current_pos, contig_name = None, change_current_contig = False ):
        self._passed_coverage.set_value( pass_value, current_pos, "?", contig_name, change_current_contig )

    def set_proportion_pass( self, pass_value, current_pos, contig_name = None, change_current_contig = False ):
        self._passed_proportion.set_value( pass_value, current_pos, "?", contig_name, change_current_contig )

    def get_coverage_pass( self, current_pos, contig_name = None ):
        return self._passed_coverage.get_value( current_pos, None, contig_name, "?" )

    def get_proportion_pass( self, current_pos, contig_name = None ):
        return self._passed_proportion.get_value( current_pos, None, contig_name, "?" )


class GenomeCollection:

    def __init__( self ):
        self._reference = None
        self._genomes = []

    def set_reference( self, reference ):
        self._reference = reference

    def reference( self ):
        return self._reference

    def get_dups_call( self, first_position, last_position = None, contig_name = None ):
        return self._reference.get_dups_call( first_position, last_position, contig_name )

    def add_genome( self, genome ):
        self._genomes.append( genome )

    def set_current_contig( self, contig_name ):
        self._reference.set_current_contig( contig_name )
        for genome in self._genomes:
            genome.set_current_contig( contig_name )

    def get_contigs( self ):
        return self._reference.get_contigs()

    def get_genome_count( self ):
        return len( self._genomes )

    # FIXME split into a larger number of smaller more testable functions
    def _format_matrix_line( self, current_contig, current_pos, matrix_format ):
        genome_count = self.get_genome_count()
        matrix_line = '' + current_contig + "::" + str( current_pos ) + "\t"
        reference_call = self._reference.get_call( current_pos, None, current_contig )
        simplified_refcall = Genome.simple_call( reference_call )
        matrix_line += '' + reference_call + "\t"
        call_data = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'indel': 0, 'snpcall': 0, 'indelcall': 0, 'refcall': 0, 'callstring': '', 'covstring': '', 'propstring': '', 'called': 0, 'passcov': 0, 'passprop': 0 }
        consensus_check = {}
        # The expensive loop, single threaded and runs for every sample-chromosome-position
        for genome in self._genomes:
            sample_call = genome.get_call( current_pos, None, current_contig, 'X' )
            simplified_sample_call = Genome.simple_call( sample_call )
            matrix_line += '' + sample_call + "\t"
            call_data[simplified_sample_call] += 1
            genome_nickname = genome.nickname()
            if genome_nickname in consensus_check:
                if consensus_check[genome_nickname] != simplified_sample_call:
                    consensus_check[genome_nickname] = 'N'
            else:
                consensus_check[genome_nickname] = simplified_sample_call
            was_called = False
            if simplified_sample_call != 'N':
                call_data['callstring'] += 'Y'
                call_data['called'] += 1
                was_called = True
            else:
                call_data['callstring'] += 'N'
            passed_coverage = genome.get_coverage_pass( current_pos, current_contig )
            call_data['covstring'] += '' + passed_coverage
            if passed_coverage == 'Y' or passed_coverage == '-':
                passed_coverage = True
                call_data['passcov'] += 1
            else:
                passed_coverage = False
            passed_proportion = genome.get_proportion_pass( current_pos, current_contig )
            call_data['propstring'] += '' + passed_proportion
            if passed_proportion == 'Y' or passed_proportion == '-':
                passed_proportion = True
                call_data['passprop'] += 1
            else:
                passed_proportion = False
            if was_called and passed_coverage and passed_proportion and simplified_refcall != 'N':
                if simplified_refcall == simplified_sample_call:
                    call_data['refcall'] += 1
                elif simplified_sample_call != 'N':
                    call_data['snpcall'] += 1
        matrix_line += '' + str( call_data['snpcall'] ) + "\t" + str( call_data['indelcall'] ) + "\t" + str( call_data['refcall'] ) + "\t"
        matrix_line += '' + str( call_data['called'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passcov'] ) + "/" + str( genome_count ) + "\t" + str( call_data['passprop'] ) + "/" + str( genome_count ) + "\t"
        matrix_line += '' + str( call_data['A'] ) + "\t" + str( call_data['C'] ) + "\t" + str( call_data['G'] ) + "\t" + str( call_data['T'] ) + "\t" + str( call_data['indel'] ) + "\t" + str( call_data['N'] ) + "\t"
        matrix_line += '' + current_contig + "\t" + str( current_pos ) + "\t"
        dups_call = self._reference.get_dups_call( current_pos, None, current_contig )
        if dups_call == "1":
            dups_call = True
        else:
            dups_call = False
        if 'N' in consensus_check.values():
            consensus_check = False
        else:
            consensus_check = True
        matrix_line += '' + str( dups_call ) + "\t" + str( consensus_check ) + "\t"
        custom_line = None
        if call_data['snpcall'] > 0 and call_data['indelcall'] == 0 and call_data['snpcall'] + call_data['refcall'] == genome_count and dups_call != True and consensus_check:
            custom_line = '' + matrix_line + "\n"
        matrix_line += '' + str( call_data['callstring'] ) + "\t" + str( call_data['covstring'] ) + "\t" + str( call_data['propstring'] ) + "\n"
        return ( matrix_line, custom_line )

    def _send_to_matrix_handles( self, master_handle, custom_handle, matrix_format ):
        master_handle.write( "LocusID\tReference\t" )
        custom_handle.write( "LocusID\tReference\t" )
        for genome in self._genomes:
            master_handle.write( '' + genome.identifier() + "\t" )
            custom_handle.write( '' + genome.identifier() + "\t" )
        master_handle.write( "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tChromosome\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\n" )
        custom_handle.write( "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tChromosome\tPosition\tInDupRegion\tSampleConsensus\n" )
        for current_contig in self._reference.get_contigs():
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


