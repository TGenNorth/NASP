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
        self._status_data[contig_name] = self._status_data[contig_name] + genome_data

    def get_contigs( self ):
        return sorted( self._status_data.keys() )

    def extend_contig( self, new_length, missing_range_filler, contig_name = None, change_current_contig = True ):
        if contig_name is None:
            contig_name = self._current_contig
        self.add_contig( contig_name, change_current_contig )
        if len( self._status_data[contig_name] ) < new_length:
            self._status_data[contig_name] = self._status_data[contig_name] + ( missing_range_filler * ( new_length - len( self._status_data[contig_name] ) ) )

    def set_value( self, new_data, first_position, missing_range_filler = "!", contig_name = None, change_current_contig = False ):
        if contig_name is None:
            contig_name = self._current_contig
        self.add_contig( contig_name, change_current_contig )
        first_position = first_position - 1
        self.extend_contig( first_position, missing_range_filler, contig_name )
        self._status_data[contig_name] = self._status_data[contig_name][:first_position] + new_data + self._status_data[contig_name][( first_position + len( new_data ) ):]

    def get_value( self, first_position, last_position = None, contig_name = None ):
        if contig_name is None:
            contig_name = self._current_contig
        if contig_name not in self.get_contigs():
            raise InvalidContigName()
        if last_position == -1:
            last_position = len( self._status_data[contig_name] )
        elif ( last_position is None ) or ( last_position < first_position ):
            last_position = first_position
        first_position = first_position - 1
        return self._status_data[contig_name][first_position:last_position]

    def get_contig_length( self, contig_name = None ):
        if contig_name is None:
            contig_name = self._current_contig
        return len( self._status_data[contig_name] )

    def _send_to_filehandle( self, output_handle, contig_prefix = "", max_chars_per_line = 80 ):
        for current_contig in self.get_contigs():
            output_handle.write( ">" + contig_prefix + current_contig + "\n" )
            if max_chars_per_line > 0:
                i = 0
                while ( max_chars_per_line * i ) < len( self._status_data[current_contig] ):
                    output_handle.write( self._status_data[current_contig][( max_chars_per_line * i ):( max_chars_per_line * ( i + 1 ) )] + "\n" )
                    i = i + 1
            else:
                output_handle.write( self._status_data[current_contig] + "\n" )

    def write_to_file( self, output_filename, contig_prefix = "", max_chars_per_line = 80 ):
        output_handle = open( output_filename, 'w' )
        self._send_to_filehandle( output_handle, contig_prefix, max_chars_per_line )
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

    def get_call( self, first_position, last_position = None, contig_name = None ):
        return self._genome.get_value( first_position, last_position, contig_name )

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

    def write_to_file( self, output_filename, contig_prefix = "", max_chars_per_line = 80 ):
        self._genome.write_to_file( output_filename, contig_prefix, max_chars_per_line )

    @staticmethod
    def generate_nickname_from_fasta_filename( fasta_filename ):
        import re
        import random
        filename_match = re.match( r'^(?:.*\/)?([^\/]+?)(?:\.[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?)?$', fasta_filename )
        if filename_match:
            fasta_nickname = filename_match.group(1)
        else:
            fasta_nickname = "fasta_" + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) ) + str( random.randrange( 10 ) )
        return fasta_nickname

    @staticmethod
    def reverse_complement( dna_string ):
        return dna_string.translate( ''.maketrans( 'ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr' ) )[::-1]


class ReferenceGenome:

    def __init__( self ):
        self._genome = Genome()
        self._dups = GenomeStatus()

    def get_contigs( self ):
        return sorted( self._genome.keys )

    def get_call( self, contig_name, first_position, last_position = None ):
        return self._genome.get_call( contig_name, first_position, last_position )

    def get_dups_call( self, first_position, last_position = None, contig_name = None ):
        return self._dups.get_value( first_position, last_position, contig_name )

    def _import_dups_line( self, line_from_dups_file, contig_prefix = "" ):
        import re
        contig_match = re.match( r'^>' + re.escape( contig_prefix ) + r'([^\s]+)(?:\s|$)', line_from_dups_file )
        if contig_match:
            self._genome.add_contig( contig_match.group(1), False )
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


class ExternalGenome:

    def __init__( self ):
        self._genome = Genome()

    def get_contigs( self ):
        return self._genome.get_contigs()

    def set_call( self, new_data, first_position, missing_range_filler = "N", contig_name = None, change_current_contig = False ):
        self._genome.set_value( new_data, first_position, missing_range_filler, contig_name, change_current_contig )

    def get_call( self, first_position, last_position = None, contig_name = None ):
        return self._genome.get_value( first_position, last_position, contig_name )


class InvalidContigName( Exception ):
    pass


def main():
    print( "This is just a class definition, to be included from other python programs." )
    print( "  ...You sillyface." )

if __name__ == "__main__": main()


