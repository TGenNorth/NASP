#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.9.1"
__email__ = "dsmith@tgen.org"


class GenomeStatus:

    def __init__( self ):
        self._status = {}
        self._current_contig = None

    def append_contig( self, genome_data, contig_name = None ):
        if contig_name is not None:
            self._current_contig = contig_name
        if self._current_contig is not None:
            if self._current_contig in self._status:
                self._status[self._current_contig] = self._status[self._current_contig] + genome_data
            else:
                self._status[self._current_contig] = genome_data
        else:
            raise InvalidContigName()

    def add_contig( self, contig_name, change_current_contig = True ):
        if change_current_contig:
            self._current_contig = contig_name
        if contig_name not in self._status:
            self._status[contig_name] = ""

    def get_contigs( self ):
        return sorted( self._status.keys() )

    def set_value( self, new_data, contig_name, start_position, missing_range_filler = "?" ):
        start_position = start_position - 1
        if len( self._status[contig_name] ) < start_position:
            self._status[contig_name] = self._status[contig_name] + ( missing_range_filler * ( start_position - len( self._status[contig_name] ) ) )
        self._status[contig_name] = self._status[contig_name][:start_position] + new_data + self._status[contig_name][( start_position + len( new_data ) ):]

    def get_value( self, contig_name, start_position, end_position = 0 ):
        if ( end_position == 0 ) or ( end_position < start_position ):
            end_position = start_position
        start_position = start_position - 1
        return self._status[contig_name][start_position:end_position]

    def get_contig_length( self, contig_name = None ):
        if contig_name is None:
            contig_name = self._current_contig
        return len( self._status[contig_name] )

    def write_to_file( self, output_filename, max_chars_per_line = 80 ):
        output_handle = open( output_filename, 'w' )
        for current_contig in self.get_contigs():
            output_handle.write( ">" + current_contig + "\n" )
            if max_chars_per_line:
                i = 0
                while ( max_chars_per_line * i ) < len( self._status[current_contig] ):
                    output_handle.write( self._status[current_contig][( max_chars_per_line * i ):( max_chars_per_line * ( i + 1 ) )] + "\n" )
                    i = i + 1
            else:
                output_handle.write( self._status[current_contig] + "\n" )
        output_handle.close()


class Genome:

    def __init__( self ):
        self._genome = GenomeStatus()

    def append_contig( self, genome_data, contig_name = None ):
        self._genome.append_contig( genome_data, contig_name )

    def add_contig( self, contig_name, change_current_contig = True ):
        self._genome.add_contig( contig_name, change_current_contig )

    def get_contigs( self ):
        return self._genome.get_contigs()

    def set_call( self, new_data, contig_name, start_position, missing_range_filler = "X" ):
        self._genome.set_value( new_data, contig_name, start_position, missing_range_filler )

    def get_call( self, contig_name, start_position, end_position = 0 ):
        return self._genome.get_value( contig_name, start_position, end_position )

    def _import_fasta_line( self, line_from_fasta ):
        import re
        contig_match = re.match( '^>([^\s]+)(?:\s|$)', line_from_fasta )
        if contig_match:
            self._genome.add_contig( contig_match.group(1) )
        else:
            data_match = re.match( '^([A-Za-z.-]+)\s*$', line_from_fasta )
            if data_match:
                self._genome.append_contig( data_match.group(1) )

    def import_fasta_file( self, fasta_filename ):
        fasta_handle = open( fasta_filename, 'r' )
        for line_from_fasta in fasta_handle:
            self._import_fasta_line( line_from_fasta )
        fasta_handle.close()

    def get_contig_length( self, contig_name = None ):
        return self._genome.get_contig_length( contig_name )


class ReferenceGenome:

    def __init__( self ):
        self._genome = Genome()
        self._dups = GenomeStatus()

    def get_contigs( self ):
        return sorted( self._genome.keys )

    def get_call( self, contig_name, start_position, end_position = 0 ):
        return self._genome.get_call( contig_name, start_position, end_position )

    def get_dups_call( self, contig_name, start_position, end_position = 0 ):
        return self._dups.get_value( contig_name, start_position, end_position )

    def _import_dups_line( self, line_from_dups_file ):
        import re
        contig_match = re.match( '^>([^\s]+)(?:\s|$)', line_from_dups_file )
        if contig_match:
            self._genome.add_contig( contig_match.group(1), False )
            self._dups.add_contig( contig_match.group(1) )
        else:
            data_match = re.match( '^([01-]+)\s*$', line_from_dups_file )
            if data_match:
                self._dups.append_contig( data_match.group(1) )

    def import_dups_file( self, dups_filename ):
        dups_handle = open( dups_filename, 'r' )
        for line_from_dups_file in dups_handle:
            self._import_dups_line( line_from_dups_file )
        dups_handle.close()


class ExternalGenome:

    def __init__( self ):
        self._genome = Genome()

    def get_contigs( self ):
        return self._genome.get_contigs()

    def set_call( self, new_data, contig_name, start_position, missing_range_filler = "N" ):
        self._genome.set_value( new_data, contig_name, start_position, missing_range_filler )

    def get_call( self, contig_name, start_position, end_position = 0 ):
        return self._genome.get_value( contig_name, start_position, end_position )


class InvalidContigName( Exception ):
    pass


def main():
    print( "This is just a class definition, to be included from other python programs." )
    print( "  ...You sillyface." )

if __name__ == "__main__": main()


