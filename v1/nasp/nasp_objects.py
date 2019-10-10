#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1.0.0"
__email__ = "dsmith@tgen.org"

import logging


class GenomeStatus(object):
    """
    Contains and manipulates any generic data that is per-contig-position.
    This could be any single type of data, like actual bases, filter data,
    depth information, insertion data, pileups, etc.
    In the perl version, this object was originally a hash of lists, and
    converted to a hash of strings for performance.  In this version, it was
    originally a dictionary of strings, and converted to a dictionary of
    lists for performance and flexibility.
    Whether single characters, numbers, boolean values, or a mixture of data
    types are used does not seem to affect memory and performance.
    Storing lists per-contig-position with this class is a complicated
    affair, as several of the manipulation functions assume you mean to
    manipulate a continuous range of positions instead of a single position
    when you do that.
    """

    """
    The conventions used for what data is stored are as follows:
    Genomes:
        A, C, G, T, U, a, c, g, t, u:  The respective call.
        N, n:  Called "N" according to upstream analysis tools.
        X:  Not called by upstream analysis tools.
        . or empty string:  A deletion relative to reference.
        String of length >1:  An insertion relative to reference.
        Any other single letter:  A degeneracy.
        !:  An algorithmic failure.
    Duplicate region data:
        0:  Position not in a region that is duplicated within the reference.
        1:  Position is in a region that is duplicated.
        -:  Duplicate checking at this position was skipped by the user.
        !:  An algorithmic failure.
    Filters:
        Y:  This position passed its filter.
        N:  This position failed its filter.
        ?:  The filter could not be checked, and so the position is assumed
            to have failed.
        -:  The filter was not applicable, or skipped, or could not be checked
            for a known reason, and so is assumed to have passed.
        !:  An algorithmic failure.
    """

    # Arrays are zero-indexed, genome positions are one-indexed. Off-by-one errors? Never heard of 'em.
    def __init__(self):
        """
        Attributes:
            _status_data: The dictionary of lists that stores the actual genome.
            data. The keys of the dictionary are the contig names.  The lists
            correspond to the position data on that contig.  Genome position is
            list position + 1.
            _current_contig: Tracks the most recently-referenced contig, for
            convenience EG reading in fastas line-by-line.
        """
        self._status_data = {}
        self._current_contig = None

    # NOTE(jtravis): Does not raise exception if contig name is None or the empty string as documented
    def add_contig(self, contig_name):
        """
        Defines a new empty contig in the genome.
        By default, if an unrecognized contig is encountered, a new empty
        contig will be created and then acted upon.
        Otherwise, add_contig must be called on a new contig first, or an
        InvalidContigName will be thrown.

        Args:
            contig_name (str): Unique contig description.

        Raises:
            InvalidContigName: If contig_name is undefined.
        """
        if contig_name not in self._status_data:
            self._status_data[contig_name] = []
        self._current_contig = contig_name

    # NOTE(jtravis): unused parameter create_contig
    def set_current_contig(self, contig_name, create_contig=True):
        """
        Sets the most-recently-referenced contig without actually performing
        any action on the data.
        Can be called to return the current contig without changing it if
        given a contig_name of None.
        Will create the contig if it has not been encountered yet by
        default, or throw an InvalidContigName otherwise.

        Args:
            contig_name (str): Unique contig description or None to query the current contig name.
            create_contig (bool): If True and the contig does not exist, an empty contig will be created.

        Returns:
            str: Name of the last accessed contig or None.

        Raises:
            InvalidContigName: If create_contig is False and the contig does not exist.
        """
        if contig_name is None:
            contig_name = self._current_contig
        elif contig_name in self._status_data:
            self._current_contig = contig_name
        elif create_contig:
            self.add_contig(contig_name)
        else:
            raise InvalidContigName(contig_name, self.get_contigs())
        return contig_name

    def get_contigs(self):
        """
        Returns:
            list: Sorted list of contig names.
        """
        return sorted(self._status_data.keys())

    def append_contig(self, genome_data, contig_name=None):
        """
        Places the passed-in data at the position following the last
        defined position on the contig.  If passed a list, will give each
        item in the list its own position.

        Args:
            genome_data (list): List of nucleotide symbols.
            contig_name (str): Unique contig description.
        """
        contig_name = self.set_current_contig(contig_name)
        self._status_data[contig_name].extend(genome_data)

    def extend_contig(self, new_length, missing_range_filler, contig_name=None):
        """
        Ensures the contig is at least new_length positions long

        Args:
            new_length (int): Minimum contig length.
            missing_range_filler (str): Placeholder character for undefined areas at the end of the contig.
            contig_name (str): Unique contig description.
        """
        contig_name = self.set_current_contig(contig_name)
        if len(self._status_data[contig_name]) < new_length:
            self._status_data[contig_name].extend(
                [missing_range_filler] * ( new_length - len(self._status_data[contig_name]) ))

    def set_value(self, new_data, position_number, missing_range_filler="!", contig_name=None):
        """
        Sets the value at position_number on the contig.
        If passed a list, will change the continuous range of positions
        starting at position_number, one position per list item.
        Will extend the contig with missing_range_filler filling undefined
        values if the position to set is beyond the end of the contig.

        Args:
            new_data (str or list): Single or list of nucleotide symbols.
            position_number (int): 1-indexed contig position number.
            missing_range_filler (str): Filler for undefined regions before the set value. Modifies the data.
            contig_name (str): Unique contig description
        """
        contig_name = self.set_current_contig(contig_name)
        self.extend_contig(position_number, missing_range_filler, contig_name)
        if len(new_data) > 1:
            self._status_data[contig_name][position_number - 1:position_number - 1 + len(new_data)] = new_data
        else:
            self._status_data[contig_name][position_number - 1] = new_data

    def get_value(self, first_position, last_position=None, contig_name=None, filler_value=None):
        """
        Args:
            contig_name (str): Unique contig description.
            first_position (int): 1-indexed first position number.
            last_position (int): Optional last position to select a range or -1 to specify the end of the contig.
            filler_value (str): Optional filler for undefined regions beyond the genome data. Does not modify the data.

        Returns:
            Returns the nucleotide at first_position, list of values from
            first_position to last_position inclusive, or None.
        """
        contig_name = self.set_current_contig(contig_name)
        queried_value = filler_value
        if last_position is None:
            if first_position <= len(self._status_data[contig_name]):
                queried_value = self._status_data[contig_name][first_position - 1]
        else:
            queried_value = []
            if last_position == -1:
                last_position = len(self._status_data[contig_name])
            if last_position >= first_position and first_position <= len(self._status_data[contig_name]):
                queried_value = self._status_data[contig_name][first_position - 1:last_position]
                if filler_value is not None and len(queried_value) < last_position - first_position + 1:
                    queried_value.extend([filler_value] * ( last_position - first_position + 1 - len(queried_value) ))
        return queried_value

    def get_contig_length(self, contig_name=None):
        """
        Args:
            contig_name (str): Unique contig description.

        Returns:
            int: Number of positions defined in the contig
        """
        contig_name = self.set_current_contig(contig_name)
        return len(self._status_data[contig_name])

    def send_to_fasta_handle(self, output_handle, contig_prefix="", max_chars_per_line=80):
        """
        Assumes the genome data is in string format or stringifiable and one
        character per position, and then writes it in to the handle open for
        writing.  The file format is like a typical fasta were the genome
        data to be base calls (but no checks are performed).

        The contigs are sorted by name, not the order they were created.

        Args:
            output_handle (file object): File to append FASTA string.
            contig_prefix (str): Prefix for all contig names.
            max_chars_per_line (int): A positive value will limit the max chars per line.
        """
        for current_contig in self.get_contigs():
            output_handle.write(">" + contig_prefix + current_contig + "\n")
            if max_chars_per_line > 0:
                i = 0
                while ( max_chars_per_line * i ) < len(self._status_data[current_contig]):
                    output_handle.write(''.join(self._status_data[current_contig][
                                                ( max_chars_per_line * i ):( max_chars_per_line * ( i + 1 ) )]) + "\n")
                    i += 1
            else:
                output_handle.write(''.join(self._status_data[current_contig]) + "\n")

    def write_to_fasta_file(self, output_filename, contig_prefix="", max_chars_per_line=80):
        """
        Opens the passed filename and passes to send_to_fasta_handle.
        This is a separate function so that unit testing is easier, and
        file names or open file handles can be used as destinations.

        Args:
            output_filename (str): Output filename.
            contig_prefix (str): Prefix for all contig names.
            max_chars_per_line (int): A positive value will limit the max contig chars per line.
        """
        with open(output_filename, 'w') as output_handle:
            self.send_to_fasta_handle(output_handle, contig_prefix, max_chars_per_line)


class Genome(GenomeStatus):
    """
    A special type of GenomeStatus where the genome information being stored
    is always actual base calls, as strings.
    """

    def __init__(self):
        """
        Attributes:
            _genome is an alias of _status_data, provided for code clarity
            when working with an actual genome
        """
        GenomeStatus.__init__(self)
        self._genome = self._status_data

    def set_call(self, new_data, first_position, missing_range_filler="X", contig_name=None):
        """ Alias of set_value, for code clarity """
        self.set_value(new_data, first_position, missing_range_filler, contig_name)

    def get_call(self, first_position, last_position=None, contig_name=None, filler_value="X"):
        """ Alias of get_value, for code clarity """
        return self.get_value( first_position, last_position, contig_name, filler_value)

    # NOTE: contig_prefix is unused
    def _import_fasta_line(self, line_from_fasta, contig_prefix=""):
        """
        Assumes the string passed in is a line from a fasta file, and
        populates the genome with the information contained.  Not meant
        to be called on any data except a full fasta file in order by
        line top to bottom.

        Args:
            line_from_fasta (str): the current line to parse
            contig_prefix (str): the prefix will be removed from the parsed contig name
        """
        import re

        # Parse the contig name discarding the prefix and surrounding whitespace characters
        contig_match = re.match(r'^>' + re.escape(contig_prefix) + r'([^\s]+)(?:\s|$)', line_from_fasta)
        if contig_match:
            self.add_contig(contig_match.group(1))
        else:
            # Parse the contig sequence discarding trailing whitespace characters
            data_match = re.match(r'^([A-Za-z.-]+)\s*$', line_from_fasta)
            if data_match:
                self.append_contig(list(data_match.group(1)))

    # contig_prefix is used by vcf_to_matrix to discard the frankenfasta contig name prefix.
    def import_fasta_file(self, fasta_filename, contig_prefix=""):
        """ Read in a fasta file.

        Args:
            fasta_filename (str): fasta file to import
            contig_prefix (str): the prefix will be removed from the parsed contig names
        """
        with open(fasta_filename, 'r') as fasta_handle:
            for line_from_fasta in fasta_handle:
                self._import_fasta_line(line_from_fasta, contig_prefix)

    @staticmethod
    def reverse_complement(dna_string):
        """
        Args:
            dna_string (str): nucleotide sequence to reverse complement

        Returns:
            string: nucleotide sequence reverse complement
        """
        return dna_string.translate(
            ''.maketrans('ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr'))[::-1]

    @staticmethod
    def simple_call(dna_string, allow_x=False, allow_del=False):
        """
        Standardizes the DNA call assumed to be the base at position one.
        Discards insertion data, changes 'U' to 'T', and changes degeneracies to 'N'.
        'X' and deletes are changed to 'N' by default.

        Args:
            dna_string (str): only the first position is considered
            allow_x (bool):
            allow_del (bool):

        Returns:
            string: 'A', 'C', 'G', 'T', or 'N' with optional 'X' and '.'
        """
        simple_base = 'N'
        if len(dna_string) > 0:
            simple_base = dna_string[0].upper()
        elif allow_del:
            simple_base = '.'
        if simple_base == 'U':
            simple_base = 'T'
        if simple_base not in ['A', 'C', 'G', 'T', 'X', '.']:
            simple_base = 'N'
        if not allow_x and ( simple_base == 'X' ):
            simple_base = 'N'
        if not allow_del and ( simple_base == '.' ):
            simple_base = 'N'
        return simple_base


class GenomeMeta(object):
    """ Stores the metadata associated with a genome.  """

    def __init__(self):
        """
        Attributes:
            _nickname: is our best guess at the name information as supplied by the
            input files.  Its determination depends on filetype.  For example, the
            nickname of a genome from a fasta file is the filename minus the
            extension.  There is no expectation that this value is unique.  A
            ( _file_path, _nickname ) tuple or the like would stand a much higher
            chance of being unique.
            _filepath:
            _file_type:
            _generators: A list of analysis tools that have been run on input files to produce
            this data, from earliest to latest.
        """
        self._nickname = None
        self._file_path = None
        self._file_type = None
        self._generators = []  # TODO: This should probably be a dictionary someday

    # NOTE(jtravis): side effect alters nickname
    def set_file_path(self, file_path):
        """
        Args:
            file_path (str):
        """
        self._file_path = file_path
        # TODO: merge and replace if-statement with set_nickname()
        if self._nickname is None:
            self._nickname = GenomeMeta.generate_nickname_from_filename(file_path)

    def set_file_type(self, type_string):
        """
        Args:
            type_string (str):
        """
        self._file_type = type_string

    def set_nickname(self, nickname):
        """
        This value should be unique per-file, for multi-sample input files,
        but there is no expectation that this would be unique per run, and
        nothing to enforce even per-file uniqueness.

        Args:
            nickname (str):
        """
        self._nickname = nickname

    def add_generators(self, generator_array):
        """
        Args:
            generator_array: A list of analysis tools that have been run on input files to produce
            this data, from earliest to latest.
        """
        self._generators.extend(generator_array)

    def file_path(self):
        """
        Returns:
            string:
        """
        return self._file_path

    def file_type(self):
        """
        Returns:
            string:
        """
        return self._file_type

    def nickname(self):
        """
        Returns:
            string:
        """
        return self._nickname

    def identifier(self):
        """
        Returns:
            string: meant to be recognizable to the user to differentiate
            their sample-analyses.  There can be no assumption that this value is
            unique, because the nickname is often not unique.  As a result, this
            should not be used in-code (EG dictionary keys) to differentiate
            samples.  The best we can do is a ( _file_path, _nickname ) tuple, and
            even that may not be guaranteed to be unique.
        """
        identifier = self._nickname
        if len(self._generators) > 0:
            identifier = "{0}::{1}".format(identifier, ','.join(self._generators))
        return identifier

    @staticmethod
    def generate_nickname_from_filename(filename):
        """
        For single-sample input files that don't carry any sample name
        metadata within the file, generate a nickname for the sample by
        removing the extension.  If this fails, generate a random name in the
        format "file_XXXXXXXX" where X is an 8-digit random integer.

        Args:
            filename (str):

        Returns:
            string: filename sans extension or "file_XXXXXXXX" where X is an 8-digit random integer.
        """
        import re
        import random
        # Parse basename from fasta or vcf file
        filename_match = re.match(r'^(?:.*/)?([^/]+?)\.(?:(?:franken)?fas?(?:ta)?|vcf)?$', filename, re.IGNORECASE)
        if filename_match:
            nickname = filename_match.group(1)
        else:
            # postfix file with 8 digit random number where 0's are allowed
            nickname = "file_" + "%08d" % random.randrange(10 ** 8)
        return nickname

    @staticmethod
    def reverse_complement(dna_string):
        return dna_string.translate(
            ''.maketrans('ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr'))[::-1]


class IndelList(object):
    """
    For storing indel data separately from the reference-indexed position
    data.  Mostly a relic from when calls were stored as a long string with
    one character per position.  Might still have some use for indel
    implementation, but might be no longer useful.
    """
    def __init__(self):
        """
        Attributes:
            _indels (dict):
        """
        self._indels = {}


class ReferenceGenome(Genome):
    """
    A special type of genome that is to be used as our reference.
    Unlike other genomes, we know there will only be one, it carries duplicate
    region data, and we don't need to store any metadata about it.
    """

    def __init__(self):
        """
        Attributes:
            _dups (GenomeStatus): carries data about whether a particular
            region of the reference was found to be very similar to another region
            in the same reference.
        """
        Genome.__init__(self)
        self._dups = GenomeStatus()

    def get_dups_call(self, first_position, last_position=None, contig_name=None):
        """
        Args:
            first_position (int):
            last_position (int):
            contig_name (str):

        Returns:
            list:
        """
        return self._dups.get_value(first_position, last_position, contig_name, "?")

    def _import_dups_line(self, line_from_dups_file, contig_prefix=""):
        """
        Just like importing any other fasta-like file line-by-line, but
        specific to duplicate region data.

        Args:
            line_from_dups_file (str):
            contig_prefix (str):
        """
        import re
        contig_match = re.match(r'^>' + re.escape(contig_prefix) + r'([^\s]+)(?:\s|$)', line_from_dups_file)
        if contig_match:
            self.add_contig(contig_match.group(1))
            self._dups.add_contig(contig_match.group(1))
        else:
            data_match = re.match(r'^([01-]+)\s*$', line_from_dups_file)
            if data_match:
                self._dups.append_contig(list(data_match.group(1)))

    def import_dups_file(self, dups_filename, contig_prefix=""):
        """ Wrapper for _import_dups_line for flexibility and testing.

        Args:
            dups_filename (str):
            contig_prefix (str):
        """
        with open(dups_filename, 'r') as dups_handle:
            for line_from_dups_file in dups_handle:
                self._import_dups_line(line_from_dups_file, contig_prefix)

class FastaGenome(Genome, GenomeMeta):
    """
    A special type of genome where we know the data came from a fasta file,
    and so we can omit the depth and proportion filters.  Meant to mimic
    a VCFGenome object, with filter checks hard-coded.
    """

    def __init__(self):
        Genome.__init__(self)
        GenomeMeta.__init__(self)
        self._indels = IndelList()

    # FIXME This data should be put into a real structure when the fasta code is pulled in
    def get_was_called(self, current_pos, contig_name=None):
        """
        A stand-in for the was-called filter that just makes sure the position
        isn't an "N".
        """
        return_value = "N"
        call_to_check = self.get_call(current_pos, None, contig_name, "X")
        if call_to_check != "X" and call_to_check != "N":
            return_value = "Y"
        return return_value

    def get_coverage_pass(self, current_pos, contig_name=None):
        """ This filter is not applicable. """
        return "-"

    def get_proportion_pass(self, current_pos, contig_name=None):
        """ This filter is not applicable. """
        return "-"


class VCFGenome(Genome, GenomeMeta):
    """
    A standard sample for analysis.  Has genome data, metadata, and data for
    the three filters.
    """

    def __init__(self):
        Genome.__init__(self)
        GenomeMeta.__init__(self)
        self._indels = IndelList()
        self._was_called = GenomeStatus()
        self._passed_coverage = GenomeStatus()
        self._passed_proportion = GenomeStatus()

    def set_was_called(self, pass_value, current_pos, contig_name=None):
        self._was_called.set_value(pass_value, current_pos, "N", contig_name)

    def set_coverage_pass(self, pass_value, current_pos, contig_name=None):
        self._passed_coverage.set_value(pass_value, current_pos, "?", contig_name)

    def set_proportion_pass(self, pass_value, current_pos, contig_name=None):
        self._passed_proportion.set_value(pass_value, current_pos, "?", contig_name)

    def get_was_called(self, current_pos, contig_name=None):
        return self._was_called.get_value(current_pos, None, contig_name, "N")

    def get_coverage_pass(self, current_pos, contig_name=None):
        return self._passed_coverage.get_value(current_pos, None, contig_name, "?")

    def get_proportion_pass(self, current_pos, contig_name=None):
        return self._passed_proportion.get_value(current_pos, None, contig_name, "?")


class CollectionStatistics(object):
    """
    Stores a running tally for the statistics for the run.
    Stats are in two categories: per-contig and per-sample.
    Contig stats are basic counts, and percentages based on reference length.
    Sample stats are tallied as x out of y for each position, and then
    counts for each sample, all samples, and any samples, can be computed
    automatically.
    For this reason, the object needs to know when the run moves on to the
    next position, and the flush_cumulative_stat_cache function does this.
    """

    def __init__(self):
        """
        _cumulative_cache keeps a running tally of sample stats at the
        current position, and then writes the results to _sample_stats when
        flush_cumulative_stat_cache() is called.
        """
        self._contig_stats = {}
        self._sample_stats = {}
        self._cumulative_cache = {}

    def _increment_by_contig(self, stat_id, contig_name):
        """ Makes sure the contig stat is defined, then increments it. """
        if ( stat_id, contig_name ) not in self._contig_stats:
            self._contig_stats[( stat_id, contig_name )] = 0
        self._contig_stats[( stat_id, contig_name )] += 1

    def increment_contig_stat(self, stat_id, contig_name=None):
        """
        Increments the stat for the specified contig, and automatically
        increments the count for the all-contigs tally on the same stat.
        """
        self._increment_by_contig(stat_id, contig_name)
        if contig_name is not None:
            self._increment_by_contig(stat_id, None)

    def get_contig_stat(self, stat_id, contig_name=None):
        return_value = 0
        if ( stat_id, contig_name ) in self._contig_stats:
            return_value = self._contig_stats[( stat_id, contig_name )]
        return return_value

    def _cache_cumulative_stats(self, stat_id, sample_nickname, did_pass):
        """
        Writes data to the stat cache for the current position.  stat_type
        is either "p" for pass/positive or "t" for total.  Flushing the cache
        involves checking the p/t ratio.  When updating the cache,
        failed/negative samples just increment t, while pass/positive samples
        increment p and t.
        Cumulative stats are across all samples at a position, and then once
        each for all sample-analyses on a sample.
        """
        if ( stat_id, sample_nickname, 't' ) not in self._cumulative_cache:
            self._cumulative_cache[( stat_id, sample_nickname, 't' )] = 0
            self._cumulative_cache[( stat_id, sample_nickname, 'p' )] = 0
        self._cumulative_cache[( stat_id, sample_nickname, 't' )] += 1
        if did_pass:
            self._cumulative_cache[( stat_id, sample_nickname, 'p' )] += 1

    def _increment_by_sample(self, stat_id, sample_nickname, sample_info, cum_type):
        """ Defines if necessary, then increments, the sample stat. """
        if ( stat_id, sample_nickname, sample_info, cum_type ) not in self._sample_stats:
            self._sample_stats[( stat_id, sample_nickname, sample_info, cum_type )] = 0
        self._sample_stats[( stat_id, sample_nickname, sample_info, cum_type )] += 1

    def record_sample_stat(self, stat_id, sample_nickname, sample_identifier, sample_path, did_pass):
        """
        Updates the sample stat, and then the cumulative stat cache.
        Increments the cache for all samples, and for all analyses
        on this sample.
        """
        if did_pass:
            self._increment_by_sample(stat_id, sample_nickname, ( sample_identifier, sample_path ), None)
        self._cache_cumulative_stats(stat_id, sample_nickname, did_pass)
        self._cache_cumulative_stats(stat_id, None, did_pass)

    def get_sample_stat(self, stat_id, sample_nickname, sample_identifier, sample_path):
        return_value = 0
        if ( stat_id, sample_nickname, ( sample_identifier, sample_path ), None ) in self._sample_stats:
            return_value = self._sample_stats[( stat_id, sample_nickname, ( sample_identifier, sample_path ), None )]
        return return_value

    def get_cumulative_stat(self, stat_id, cum_type, sample_nickname=None):
        return_value = 0
        if ( stat_id, sample_nickname, None, cum_type ) in self._sample_stats:
            return_value = self._sample_stats[( stat_id, sample_nickname, None, cum_type )]
        return return_value

    def flush_cumulative_stat_cache(self):
        """
        Does the math on the p/t ratio for the stat cache for the current
        position, and writes that to the sample stats for the any/all counts.
        p > 0:  any++
        p = t:  all++
        """
        for ( stat_id, sample_nickname, stat_type ) in self._cumulative_cache:
            if stat_type == 'p':
                if self._cumulative_cache[( stat_id, sample_nickname, 'p' )] == self._cumulative_cache[
                    ( stat_id, sample_nickname, 't' )]:
                    self._increment_by_sample(stat_id, sample_nickname, None, 'all')
                if self._cumulative_cache[( stat_id, sample_nickname, 'p' )] > 0:
                    self._increment_by_sample(stat_id, sample_nickname, None, 'any')
        self._cumulative_cache = {}


class GenomeCollection(CollectionStatistics):
    """
    A "master matrix" object, of sorts.
    Carries all the data necessary to make a matrix, although some of it
    is computed on-the-fly as the matrix is actually written, for
    performance reasons.  Stats aren't available until after the matrix
    is written, for this reason.
    """

    def __init__(self):
        """
        _failed_genomes is just a list of filenames that had errors during
        import.  Added as blank columns to inform the user.
        """
        CollectionStatistics.__init__(self)
        self._reference = None
        self._genomes = []
        self._genome_identifiers = {}
        self._failed_genomes = []

    @staticmethod
    def _get_key(genome):
        """
        A wrapper getter for a genome's identifier, for the purpose of passing
        informative data to the sort function.  This is used to make sample
        order predictable and preserved between different runs.
        """
        return genome.identifier()

    def set_reference(self, reference):
        self._reference = reference

    def reference(self):
        return self._reference

    def get_dups_call(self, first_position, last_position=None, contig_name=None):
        return self._reference.get_dups_call(first_position, last_position, contig_name)

    def add_genome(self, genome):
        """
        Adds the genome to the collection, then makes sure the genome list is
        properly set and in order.
        """
        self._genomes.append(genome)
        genome_nickname = genome.nickname()
        if genome_nickname not in self._genome_identifiers:
            self._genome_identifiers[genome_nickname] = {}
        self._genome_identifiers[genome_nickname][( genome.identifier(), genome.file_path() )] = True
        self._genomes.sort(key=GenomeCollection._get_key)

    def add_failed_genome(self, genome_path):
        if genome_path not in self._failed_genomes:
            self._failed_genomes.append(genome_path)

    def set_current_contig(self, contig_name):
        contig_name = self._reference.set_current_contig(contig_name)
        for genome in self._genomes:
            genome.set_current_contig(contig_name)
        return contig_name

    def get_contigs(self):
        return self._reference.get_contigs()

    # FIXME this function is starting to become a bit of a cluster
    def _format_matrix_line(self, current_contig, current_pos, matrix_formats, pattern_data):
        """
        A matrix line represents all the sample data at one reference
        contig-position.
        This function contains a large portion of the ultimate NASP logic that
        goes into generating a matrix, like filter application and consensus
        checking.
        A section below must be executed linearly and single-threaded for all
        contig-position-sample-analyses.  This is the "expensive loop".
        The monolithic nature of this function, its massive size, and the fact
        that it makes little distinction between calculating something and
        writing something to a file, are all problems that will probably need
        to eventually be corrected.  I fear, however, that there is a
        significant tradeoff with it works / it's fast to run / it was fast to
        develop / it's easily maintained / it's beautiful.
        """
        all_matrices = []
        all_vcfs = []
        allcallable_matrices = []
        snp_matrices_best = []
        snp_matrices_md = []
        fasta_pending_data = {}
        vcf_pending_data = []
        current_pattern = ''
        # Key None stores next unused, value 1 reserved for reference
        current_pattern_legend = {None: 2}
        encountered_calls = []
        for matrix_format in matrix_formats:
            if matrix_format['dataformat'] == 'matrix':
                all_matrices.append(matrix_format)
                if matrix_format['filter'] == 'allcallable':
                    allcallable_matrices.append(matrix_format)
                elif matrix_format['filter'] == 'bestsnp' or matrix_format['filter'] == 'includeref':
                    snp_matrices_best.append(matrix_format)
                elif matrix_format['filter'] == 'missingdata':
                    snp_matrices_md.append(matrix_format)
            elif matrix_format['dataformat'] == 'vcf':
                all_vcfs.append(matrix_format)
        genome_count = len(self._genomes)
        failed_genome_count = len(self._failed_genomes)
        for matrix_format in all_matrices:
            matrix_format['linetowrite'] = "{0}::{1}\t".format(current_contig, str(current_pos))
        for matrix_format in all_vcfs:
            matrix_format['linetowrite'] = "{0}\t{1}\t.\t".format(current_contig, str(current_pos))
        reference_call = self._reference.get_call(current_pos, None, current_contig)
        simplified_refcall = Genome.simple_call(reference_call)
        fasta_pending_data['Reference'] = simplified_refcall
        if simplified_refcall == 'N':
            current_pattern += 'N'
        else:
            current_pattern_legend[simplified_refcall] = 1
            current_pattern += '1'
        self.increment_contig_stat('reference_length', current_contig)
        if simplified_refcall != 'N':
            self.increment_contig_stat('reference_clean', current_contig)
        for matrix_format in ( all_matrices + all_vcfs ):
            matrix_format['linetowrite'] += "{0}\t".format(reference_call)
        dups_call = self._reference.get_dups_call(current_pos, None, current_contig)
        if dups_call == "1":
            dups_call = True
            self.increment_contig_stat('reference_duplicated', current_contig)
        else:
            dups_call = False
        call_data = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'indel': 0, 'snpcall': 0, 'indelcall': 0, 'refcall': 0,
                     'callstring': '', 'covstring': '', 'propstring': '', 'called': 0, 'passcov': 0, 'passprop': 0}
        consensus_check = {}
        # The expensive loop, single threaded and runs for every sample-analysis-contig-position
        for genome in self._genomes:
            sample_call = genome.get_call(current_pos, None, current_contig, 'X')
            simplified_sample_call = Genome.simple_call(sample_call)
            call_data[simplified_sample_call] += 1
            genome_nickname = genome.nickname()
            genome_identifier = genome.identifier()
            genome_path = genome.file_path()
            was_called = genome.get_was_called(current_pos, current_contig)
            call_data['callstring'] += was_called
            if was_called == 'Y':
                was_called = True
                call_data['called'] += 1
            else:
                was_called = False
            self.record_sample_stat('was_called', genome_nickname, genome_identifier, genome_path, was_called)
            passed_coverage = genome.get_coverage_pass(current_pos, current_contig)
            call_data['covstring'] += passed_coverage
            if passed_coverage == 'Y' or passed_coverage == '-':
                passed_coverage = True
                call_data['passcov'] += 1
            else:
                passed_coverage = False
            self.record_sample_stat('passed_coverage_filter', genome_nickname, genome_identifier, genome_path,
                                    passed_coverage)
            passed_proportion = genome.get_proportion_pass(current_pos, current_contig)
            call_data['propstring'] += passed_proportion
            if passed_proportion == 'Y' or passed_proportion == '-':
                passed_proportion = True
                call_data['passprop'] += 1
            else:
                passed_proportion = False
            self.record_sample_stat('passed_proportion_filter', genome_nickname, genome_identifier, genome_path,
                                    passed_proportion)
            if was_called and passed_coverage and passed_proportion:
                if genome_nickname in consensus_check:
                    if consensus_check[genome_nickname] != simplified_sample_call:
                        consensus_check[genome_nickname] = 'N'
                else:
                    consensus_check[genome_nickname] = simplified_sample_call
            else:
                consensus_check[genome_nickname] = 'N'
            # FIXME indels
            if was_called and passed_coverage and passed_proportion and simplified_refcall != 'N':
                if not dups_call:
                    self.record_sample_stat('quality_breadth', genome_nickname, genome_identifier, genome_path, True)
                if simplified_refcall == simplified_sample_call:
                    call_data['refcall'] += 1
                    if not dups_call:
                        self.record_sample_stat('called_reference', genome_nickname, genome_identifier, genome_path,
                                                True)
                        self.record_sample_stat('called_snp', genome_nickname, genome_identifier, genome_path, False)
                        self.record_sample_stat('called_indel', genome_nickname, genome_identifier, genome_path, False)
                        self.record_sample_stat('called_degen', genome_nickname, genome_identifier, genome_path, False)
                elif simplified_sample_call != 'N':
                    call_data['snpcall'] += 1
                    if not dups_call:
                        self.record_sample_stat('called_reference', genome_nickname, genome_identifier, genome_path,
                                                False)
                        self.record_sample_stat('called_snp', genome_nickname, genome_identifier, genome_path, True)
                        self.record_sample_stat('called_indel', genome_nickname, genome_identifier, genome_path, False)
                        self.record_sample_stat('called_degen', genome_nickname, genome_identifier, genome_path, False)
                else:
                    if not dups_call:
                        self.record_sample_stat('called_reference', genome_nickname, genome_identifier, genome_path,
                                                False)
                        self.record_sample_stat('called_snp', genome_nickname, genome_identifier, genome_path, False)
                        self.record_sample_stat('called_indel', genome_nickname, genome_identifier, genome_path, False)
                        self.record_sample_stat('called_degen', genome_nickname, genome_identifier, genome_path, True)
            elif simplified_refcall != 'N':
                if not dups_call:
                    self.record_sample_stat('quality_breadth', genome_nickname, genome_identifier, genome_path, False)
            for matrix_format in ( allcallable_matrices + snp_matrices_best ):
                matrix_format['linetowrite'] += "{0}\t".format(sample_call)
            fasta_pending_data[genome_identifier] = simplified_sample_call
            vcf_current_data = { 'GT': '.', 'was_called': was_called, 'passed_coverage': passed_coverage, 'passed_proportion': passed_proportion }
            if was_called and passed_coverage and passed_proportion and simplified_sample_call != 'N':
                if simplified_sample_call not in current_pattern_legend:
                    current_pattern_legend[simplified_sample_call] = current_pattern_legend[None]
                    current_pattern_legend[None] += 1
                    # FIXME indels
                    encountered_calls.append(simplified_sample_call)
                current_pattern += str(current_pattern_legend[simplified_sample_call])
                vcf_current_data['GT'] = current_pattern_legend[simplified_sample_call] - 1
                for matrix_format in snp_matrices_md:
                    matrix_format['linetowrite'] += "{0}\t".format(sample_call)
            elif not was_called:
                current_pattern += 'N'
                fasta_pending_data[genome_identifier] = 'N'
                for matrix_format in snp_matrices_md:
                    matrix_format['linetowrite'] += "X\t"
            else:
                current_pattern += 'N'
                fasta_pending_data[genome_identifier] = 'N'
                for matrix_format in snp_matrices_md:
                    matrix_format['linetowrite'] += "N\t"
            vcf_pending_data.append(vcf_current_data)
        for matrix_format in all_matrices:
            matrix_format['linetowrite'] += "\t" * failed_genome_count
        for genome_nickname in consensus_check:
            if consensus_check[genome_nickname] != 'N':
                self.record_sample_stat('consensus', genome_nickname, None, None, True)
            else:
                self.record_sample_stat('consensus', genome_nickname, None, None, False)
        for matrix_format in all_matrices:
            matrix_format['linetowrite'] += "{0}\t{1}\t{2}\t".format(str(call_data['snpcall']),
                                                                     str(call_data['indelcall']),
                                                                     str(call_data['refcall']))
            matrix_format['linetowrite'] += "{0}/{3}\t{1}/{3}\t{2}/{3}\t".format(str(call_data['called']),
                                                                                 str(call_data['passcov']),
                                                                                 str(call_data['passprop']),
                                                                                 str(genome_count))
            matrix_format['linetowrite'] += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(str(call_data['A']),
                                                                                    str(call_data['C']),
                                                                                    str(call_data['G']),
                                                                                    str(call_data['T']),
                                                                                    str(call_data['indel']),
                                                                                    str(call_data['N']))
            matrix_format['linetowrite'] += "{0}\t{1}\t".format(current_contig, str(current_pos))
        if 'N' not in consensus_check.values():
            consensus_check = True
            self.increment_contig_stat('all_passed_consensus', current_contig)
        else:
            consensus_check = False
        if call_data['called'] == genome_count:
            self.increment_contig_stat('all_called', current_contig)
        if call_data['passcov'] == genome_count:
            self.increment_contig_stat('all_passed_coverage', current_contig)
        if call_data['passprop'] == genome_count:
            self.increment_contig_stat('all_passed_proportion', current_contig)
        if consensus_check and not dups_call and call_data['called'] == genome_count and call_data[
            'passcov'] == genome_count and call_data['passprop'] == genome_count and call_data['N'] == 0:
            self.increment_contig_stat('quality_breadth', current_contig)
            if call_data['snpcall'] > 0:
                self.increment_contig_stat('best_snps', current_contig)
        if not dups_call and call_data['snpcall'] > 0:
            self.increment_contig_stat('any_snps', current_contig)
        for matrix_format in all_matrices:
            matrix_format['linetowrite'] += "{0}\t{1}\t".format(str(dups_call), str(consensus_check))
        for matrix_format in ( allcallable_matrices + snp_matrices_md ):
            matrix_format['linetowrite'] += "{0}\t{1}\t{2}\t".format(str(call_data['callstring']), str(call_data['covstring']), str(call_data['propstring']))
        for matrix_format in all_matrices:
            if current_pattern not in pattern_data:
                pattern_data[current_pattern] = pattern_data[None]
                pattern_data[None] += 1
            matrix_format['linetowrite'] += "'{0}'\t{1}\n".format(current_pattern, str(pattern_data[current_pattern]))
        for matrix_format in all_vcfs:
            if len(encountered_calls) > 0:
                matrix_format['linetowrite'] += ",".join(encountered_calls)
            else:
                matrix_format['linetowrite'] += "."
            matrix_format['linetowrite'] += "\t.\tPASS\tAN={0};NS={1}\tGT:FT".format(len(encountered_calls)+1, str(call_data['snpcall']+call_data['indelcall']+call_data['refcall']))
            for vcf_current_data in vcf_pending_data:
                matrix_format['linetowrite'] += "\t{0}:".format(str(vcf_current_data['GT']))
                if not vcf_current_data['was_called']:
                    matrix_format['linetowrite'] += "NoCall"
                elif not vcf_current_data['passed_coverage']:
                    matrix_format['linetowrite'] += "CovFail"
                elif not vcf_current_data['passed_proportion']:
                    matrix_format['linetowrite'] += "PropFail"
                else:
                    matrix_format['linetowrite'] += "PASS"
            matrix_format['linetowrite'] += "\n"
        # Determine if a line should be present in a particular matrix
        # bestsnp outputs
        if call_data['snpcall'] == 0 or call_data['indelcall'] > 0 or call_data['snpcall'] + call_data[
            'refcall'] < genome_count or dups_call or not consensus_check:
            for matrix_format in matrix_formats:
                if matrix_format['filter'] == 'bestsnp':
                    if matrix_format['dataformat'] == 'matrix' or matrix_format['dataformat'] == 'vcf':
                        matrix_format['linetowrite'] = None
        else:
            for matrix_format in matrix_formats:
                if matrix_format['dataformat'] == 'fasta' and matrix_format['filter'] == 'bestsnp':
                    for genome_identifier in fasta_pending_data:
                        matrix_format['fastadata'].append_contig(fasta_pending_data[genome_identifier],
                                                                 genome_identifier)
        # missing data outputs
        if call_data['snpcall'] == 0 or call_data['indelcall'] > 0 or dups_call:
            for matrix_format in matrix_formats:
                if matrix_format['filter'] == "missingdata":
                    if matrix_format['dataformat'] == 'matrix' or matrix_format['dataformat'] == 'vcf':
                        matrix_format['linetowrite'] = None
        else:
            for matrix_format in matrix_formats:
                if matrix_format['dataformat'] == 'fasta' and matrix_format['filter'] == 'missingdata':
                    for genome_identifier in fasta_pending_data:
                        matrix_format['fastadata'].append_contig(fasta_pending_data[genome_identifier],
                                                                 genome_identifier)
        # including ref outputs
        if call_data['indelcall'] > 0 or call_data['snpcall'] + call_data[
            'refcall'] < genome_count or dups_call or not consensus_check:
            for matrix_format in matrix_formats:
                if matrix_format['filter'] == 'includeref':
                    if matrix_format['dataformat'] == 'matrix' or matrix_format['dataformat'] == 'vcf':
                        matrix_format['linetowrite'] = None
        else:
            for matrix_format in matrix_formats:
                if matrix_format['dataformat'] == 'fasta' and matrix_format['filter'] == 'includeref':
                    for genome_identifier in fasta_pending_data:
                        matrix_format['fastadata'].append_contig(fasta_pending_data[genome_identifier],
                                                                 genome_identifier)
        self.flush_cumulative_stat_cache()

    # NOTE(jtravis): Fix documentation _write_matrix_line() does not exist
    def send_to_matrix_handles(self, matrix_formats):
        """
        Writes headers and handles per-matrix logic.  Calls _write_matrix_line
        to handle the per-line computation and analysis.
        """
        for matrix_format in matrix_formats:
            if matrix_format['dataformat'] == 'matrix':
                matrix_format['handle'].write("LocusID\tReference\t")
                for genome in self._genomes:
                    matrix_format['handle'].write("{0}\t".format(genome.identifier()))
                for genome_path in self._failed_genomes:
                    matrix_format['handle'].write("{0}\t".format(genome_path))
                if matrix_format['filter'] == 'bestsnp' or matrix_format['filter'] == 'includeref':
                    matrix_format['handle'].write(
                        "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tPattern\tPattern#\n")
                else:
                    # must be all callable or missing data
                    matrix_format['handle'].write(
                        "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\t" +
                            "Contig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\tPattern#\n")
                matrix_format['linetowrite'] = ''
            elif matrix_format['dataformat'] == 'fasta':
                matrix_format['fastadata'] = GenomeStatus()
                for genome in self._genomes:
                    matrix_format['fastadata'].add_contig(genome.identifier())
            elif matrix_format['dataformat'] == 'vcf':
                matrix_format['handle'].write("##fileFormat=VCFv4.2\n##source=NASPv{0}\n".format(__version__))
                for current_contig in self._reference.get_contigs():
                    matrix_format['handle'].write("##contig=<ID=\"{0}\",length={1}>\n".format(current_contig, self._reference.get_contig_length(current_contig)))
                for genome in self._genomes:
                    matrix_format['handle'].write("##SAMPLE=<ID=\"{0}\",Genomes=\"{0}\",Mixture=1.0>\n".format(genome.identifier()))
                matrix_format['handle'].write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n")
                matrix_format['handle'].write("##FILTER=<ID=NoCall,Description=\"No call for this sample at this position\">\n")
                matrix_format['handle'].write("##FILTER=<ID=CovFail,Description=\"Insufficient depth of coverage for this sample at this position\">\n")
                matrix_format['handle'].write("##FILTER=<ID=PropFail,Description=\"Insufficient proportion of reads were variant for this sample at this position\">\n")
                matrix_format['handle'].write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                matrix_format['handle'].write("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filters that failed for this sample at this position\">\n")
                matrix_format['handle'].write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
                for genome in self._genomes:
                    matrix_format['handle'].write("\t" + genome.identifier())
                matrix_format['handle'].write("\n")
        # Key None stores next unused
        pattern_data = {None: 1}
        for current_contig in self.get_contigs():
            for current_pos in range(1, self._reference.get_contig_length(current_contig) + 1):
                self._format_matrix_line(current_contig, current_pos, matrix_formats, pattern_data)
                for matrix_format in matrix_formats:
                    if matrix_format['dataformat'] == 'matrix' and matrix_format['linetowrite'] is not None:
                        matrix_format['handle'].write(matrix_format['linetowrite'])
                        matrix_format['linetowrite'] = ''
                    if matrix_format['dataformat'] == 'vcf' and matrix_format['linetowrite'] is not None:
                        matrix_format['handle'].write(matrix_format['linetowrite'])
                        matrix_format['linetowrite'] = ''
        for matrix_format in matrix_formats:
            if matrix_format['dataformat'] == 'fasta':
                matrix_format['fastadata'].send_to_fasta_handle(matrix_format['handle'])

    def write_to_matrices(self, matrix_formats):
        """ Opens files for writing; abstracted for flexibility/testing. """
        for matrix_format in matrix_formats:
            matrix_format['handle'] = open(matrix_format['filename'], 'w')
        self.send_to_matrix_handles(matrix_formats)
        for matrix_format in matrix_formats:
            matrix_format['handle'].close()

    def _write_general_stats(self, general_handle):
        """ Writes finalized general stats to file. """
        general_stat_array = ['reference_length', 'reference_clean', 'reference_duplicated', 'all_called',
                              'all_passed_coverage', 'all_passed_proportion', 'all_passed_consensus', 'quality_breadth',
                              'any_snps', 'best_snps']
        denominator_stat = 'reference_length'
        general_handle.write("Contig\t")
        for current_stat in general_stat_array:
            general_handle.write("{0}\t".format(current_stat))
            if current_stat != denominator_stat:
                general_handle.write("{0} (%)\t".format(current_stat))
        general_handle.write("\n")
        general_handle.write("\tstat descriptions go here\n") # FIXME
        for current_contig in ( [None] + self.get_contigs() ):
            denominator_value = self.get_contig_stat(denominator_stat, current_contig)
            if current_contig is None:
                general_handle.write("Whole Genome\t")
            else:
                general_handle.write("{0}\t".format(current_contig))
            for current_stat in general_stat_array:
                general_handle.write("{0}\t".format(str(self.get_contig_stat(current_stat, current_contig))))
                if current_stat != denominator_stat:
                    general_handle.write(
                        "%.2f%%\t" % ( self.get_contig_stat(current_stat, current_contig) / denominator_value * 100 ))
            general_handle.write("\n")

    def _write_sample_stats(self, sample_handle):
        """
        Writes finalized sample stats to file.
        Includes per-sample counts, any/all counts across each set of
        sample-analyses, and any/all counts overall.
        """
        sample_stat_array = ['was_called', 'passed_coverage_filter', 'passed_proportion_filter', 'quality_breadth',
                             'called_reference', 'called_snp', 'called_degen']
        denominator_stat = 'reference_length'
        denominator_value = self.get_contig_stat(denominator_stat)
        sample_handle.write("Sample\tSample::Analysis\t")
        for current_stat in sample_stat_array:
            sample_handle.write("{0}\t".format(current_stat))
            sample_handle.write("{0} (%)\t".format(current_stat))
        sample_handle.write("\n")
        sample_handle.write("\tstat descriptions go here\n\n") # FIXME
        for current_sample in ( [None] + sorted(self._genome_identifiers.keys()) ):
            for current_analysis in ['any', 'all']:
                if current_sample is not None:
                    sample_handle.write("{0}\t".format(current_sample))
                sample_handle.write("[{0}]\t".format(current_analysis))
                if current_sample is None:
                    sample_handle.write("\t")
                for current_stat in sample_stat_array:
                    sample_handle.write("{0}\t".format(str(self.get_cumulative_stat(current_stat, current_analysis, current_sample))))
                    if current_stat != denominator_stat:
                        sample_handle.write("%.2f%%\t" % (self.get_cumulative_stat(current_stat, current_analysis, current_sample) / denominator_value * 100))
                sample_handle.write("\n")
            if current_sample is not None:
                for current_analysis in sorted(self._genome_identifiers[current_sample].keys()):
                    ( sample_identifier, sample_path ) = current_analysis
                    sample_handle.write("{0}\t{1}\t".format(current_sample, sample_identifier))
                    for current_stat in sample_stat_array:
                        sample_handle.write("{0}\t".format(str(self.get_sample_stat(current_stat, current_sample, sample_identifier, sample_path))))
                        if current_stat != denominator_stat:
                            sample_handle.write("%.2f%%\t" % (self.get_sample_stat(current_stat, current_sample, sample_identifier, sample_path) / denominator_value * 100))
                    sample_handle.write("\n")
            sample_handle.write("\n")

    def write_to_stats_files(self, general_filename, sample_filename):
        """ Opens files for writing; abstracted for flexibility/testing. """
        with open(general_filename, 'w') as general_handle, open(sample_filename, 'w') as sample_handle:
            self._write_general_stats(general_handle)
            self._write_sample_stats(sample_handle)
        # print( self._stats._contig_stats )
        # print( self._stats._sample_stats )


class VCFRecord(object):
    """ VCF parser, object representing an input VCF being read. """

    def __init__(self, file_path):
        self._file_path = file_path
        self._file_handle = open(self._file_path, 'r')
        self._header_list = []
        self._sample_list = []
        self._current_record = {}
        self._get_header_map()

    def _get_header_map(self):
        """
        Pulls the headers from the file, so that we know which columns mean
        what as we read the file in line-by-line.
        The header map is expected to start with "#CHROM".
        """
        current_line = ''
        while current_line[0:6] != "#CHROM":
            current_line = self._file_handle.readline()
            if current_line == '':
                raise MalformedInputFile(self._file_path, "mandatory VCF header not found or not recognized")
        header_list = current_line.rstrip()[1:].split("\t")
        sample_headers_started = False
        for current_header in header_list:
            self._header_list.append(current_header)
            if sample_headers_started:
                self._sample_list.append(current_header)
            elif current_header == "FORMAT":
                sample_headers_started = True
        if not sample_headers_started:
            self._sample_list.append('vcf_sample')
            # print( self._header_list )
            # print( self._sample_list )

    def _get_record_alts(self):
        """
        Get and parse the list of non-reference calls that might appear in
        one or more samples.
        """
        self._current_record['alts'] = [self._current_record['global']['REF']]
        if self._current_record['global']['ALT'] != '.':
            self._current_record['alts'] += self._current_record['global']['ALT'].split(',')

    def _get_record_info(self):
        self._current_record['info'] = {}
        info_strings = self._current_record['global']['INFO'].split(';')
        for info_string in info_strings:
            info_keyval = info_string.split('=', 1)
            if len(info_keyval) > 1:
                self._current_record['info'][info_keyval[0]] = info_keyval[1]
            else:
                self._current_record['info'][info_keyval[0]] = None

    def _get_record_samples(self):
        if 'FORMAT' in self._current_record['global']:
            self._current_record['samples'] = {}
            sample_header = self._current_record['global']['FORMAT'].split(':')
            for current_sample in self._sample_list:
                self._current_record['samples'][current_sample] = dict(
                    zip(sample_header, self._current_record['global'][current_sample].split(':')))

    def fetch_next_record(self):
        """
        We're done with the current line, let's move to the next.
        Populates all the information at the position we'll need for parsing
        parts of the next line.
        """
        current_line = "#"
        while current_line[0:1] == "#":
            current_line = self._file_handle.readline()
        return_value = False
        if current_line != '':
            record_list = current_line.rstrip().split("\t")
            self._current_record['global'] = dict(zip(self._header_list, record_list))
            self._get_record_alts()
            self._get_record_info()
            self._get_record_samples()
            return_value = True
        # print( self._current_record )
        return return_value

    def get_samples(self):
        return self._sample_list

    def get_contig(self):
        return self._current_record['global']['CHROM']

    def get_position(self):
        return int(self._current_record['global']['POS'])

    def get_reference_call(self):
        return self._current_record['global']['REF']

    def get_sample_call(self, current_sample):
        # FIXME indels
        return_value = None
        if len(self._current_record['alts']) == 1:
            return_value = self._current_record['alts'][0]
        elif 'FORMAT' in self._current_record['global']:
            if 'GT' in self._current_record['samples'][current_sample]:
                alt_number = self._current_record['samples'][current_sample]['GT'].split('/', 1)[0].split('|', 1)[0]
                if alt_number.isdigit():
                    return_value = self._current_record['alts'][int(alt_number)]
                    # OMG varscan
                    if len(self._current_record['global']['REF']) > 1 and (
                                len(self._current_record['global']['REF']) - 1 ) == len(return_value) and \
                                    self._current_record['global']['REF'][:len(return_value)] != return_value and \
                                    self._current_record['global']['REF'][-len(return_value):] == return_value:
                        return_value = self._current_record['alts'][0]
        #elif len(self._current_record['alts']) > 1:
        #    return_value = self._current_record['alts'][1]
        return return_value

    def get_coverage(self, current_sample):
        sample_coverage = None
        if 'FORMAT' in self._current_record['global'] and 'DP' in self._current_record['samples'][current_sample] and \
                        self._current_record['samples'][current_sample]['DP'] is not None and \
                self._current_record['samples'][current_sample]['DP'].isdigit():
            sample_coverage = int(self._current_record['samples'][current_sample]['DP'])
        elif 'DP' in self._current_record['info'] and self._current_record['info']['DP'] is not None and \
                self._current_record['info']['DP'].isdigit():
            sample_coverage = int(self._current_record['info']['DP']) / len(self._sample_list)
        elif 'ADP' in self._current_record['info'] and self._current_record['info']['ADP'] is not None and \
                self._current_record['info']['ADP'].isdigit():
            sample_coverage = int(self._current_record['info']['ADP']) / len(self._sample_list)
        # NASP output
        elif 'FORMAT' in self._current_record['global'] and 'FT' in self._current_record['samples'][current_sample]:
            failed_filters = self._current_record['samples'][current_sample]['FT'].split(',')
            if 'CovFail' in failed_filters:
                sample_coverage = -1
            elif 'PASS' in failed_filters or 'PropFail' in failed_filters:
                sample_coverage = 'PASS'
        return sample_coverage

    def get_proportion(self, current_sample, sample_coverage, is_a_snp):
        sample_proportion = None
        if 'FORMAT' in self._current_record['global'] and 'AD' in self._current_record['samples'][current_sample]:
            call_depths = self._current_record['samples'][current_sample]['AD'].split(',')
            # gatk, reliable and documented
            if len(call_depths) > 1:
                alt_number = self._current_record['samples'][current_sample]['GT'].split('/', 1)[0].split('|', 1)[0]
                if alt_number.isdigit():
                    sample_proportion = int(call_depths[int(alt_number)]) / sample_coverage
                    # varscan, reliable and documented
            elif is_a_snp:
                sample_proportion = int(call_depths[0]) / sample_coverage
            elif not is_a_snp and 'RD' in self._current_record['samples'][current_sample]:
                sample_proportion = int(self._current_record['samples'][current_sample]['RD']) / sample_coverage
        # solsnp, undocumented, no multi-sample support
        elif 'AR' in self._current_record['info']:
            sample_proportion = float(self._current_record['info']['AR'])
            if not is_a_snp:
                sample_proportion = 1 - sample_proportion
        # samtools, estimate, dubious accuracy
        elif 'DP4' in self._current_record['info']:
            call_depths = self._current_record['info']['DP4'].split(',')
            if is_a_snp:
                sample_proportion = ( int(call_depths[2]) + int(call_depths[3]) ) / (
                    sample_coverage * len(self._sample_list) )
            else:
                sample_proportion = ( int(call_depths[0]) + int(call_depths[1]) ) / (
                    sample_coverage * len(self._sample_list) )
        # NASP output
        elif 'FORMAT' in self._current_record['global'] and 'FT' in self._current_record['samples'][current_sample]:
            failed_filters = self._current_record['samples'][current_sample]['FT'].split(',')
            if 'PropFail' in failed_filters:
                sample_proportion = -1
            elif 'PASS' in failed_filters:
                sample_proportion = 'PASS'
        return sample_proportion

    def get_sample_info(self, current_sample):
        sample_info = {'call': self.get_sample_call(current_sample)}
        # FIXME indels
        if sample_info['call'] is not None and len(sample_info['call']) > 1:
            sample_info['call'] = sample_info['call'][0]
        sample_info['was_called'] = False
        sample_info['is_a_snp'] = False
        if sample_info['call'] is not None and sample_info['call'] != 'N':
            sample_info['was_called'] = True
            if Genome.simple_call(sample_info['call']) != Genome.simple_call(self._current_record['global']['REF']):
                sample_info['is_a_snp'] = True
        # FIXME indels
        sample_info['is_an_insert'] = None
        sample_info['is_a_delete'] = None
        sample_info['coverage'] = self.get_coverage(current_sample)
        sample_info['proportion'] = self.get_proportion(current_sample, sample_info['coverage'],
                                                        sample_info['is_a_snp'])
        return sample_info


class InvalidContigName(Exception):
    """
    A contig was referenced that is not known to exist in this run.
    For the most part, this exception isn't used because encountering
    an unknown contig either means we should create it or we should
    discard the data that follows.  However, there are potential times
    where we need to take more extreme action.
    """

    def __init__(self, invalid_contig, contig_list):
        self._invalid_contig = invalid_contig
        self._contig_list = contig_list

    def __str__(self):
        return "Contig '{0}' not in {1}!".format(self._invalid_contig, str(self._contig_list))


class ReferenceCallMismatch(Exception):
    """
    At two different points in the run, the reference call appeared, and the
    two calls disagreed.  This probably means the user is analyzing two
    different sets of data as if they belong together, in error.  This is
    a very bad thing.
    """

    def __init__(self, old_call, new_call, data_file=None, contig_name=None, position_number=None):
        self._old_call = old_call
        self._new_call = new_call
        self._data_file = data_file
        self._contig_name = contig_name
        self._position_number = position_number

    def __str__(self):
        """
        The file, contig, and position are optional.  When given in any
        combination, what's given will be in the error message.
        """
        return_value = "Expected reference call of '{0}' but got '{1}'".format(self._old_call, self._new_call)
        if self._data_file is not None:
            return_value += " while reading '{0}'".format(self._data_file)
        if self._position_number is not None:
            return_value += " at position {0}".format(str(self._position_number))
        if self._contig_name is not None:
            return_value += " on contig '{0}'".format(self._contig_name)
        return_value += "!"
        return return_value


class MalformedInputFile(Exception):
    """
    There was an error parsing the input file.
    Give the information we can to the user in the stderr stream.
    Let's hope they read it.
    """

    def __init__(self, data_file, error_message=None):
        self._data_file = data_file
        self._error_message = error_message

    def __str__(self):
        return_value = "Input file '{0}' seems to be malformed".format(self._data_file)
        if self._error_message is not None:
            return_value += ": {0}".format(self._error_message)
        return_value += "!"
        return return_value


def main():
    print("This is just a class definition, to be included from other python programs.")
    print("  ...You sillyface.")


if __name__ == "__main__":
    main()


