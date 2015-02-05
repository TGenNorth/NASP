__author__ = 'jtravis'

from abc import ABCMeta, abstractmethod
from collections import namedtuple
import csv

PositionInfo = namedtuple('PositionInfo', [
    # List of SampleInfo
    'samples',
    #
    'reference',
    # True if all samples called A/C/G/T
    'is_all_called',
    # True if reference call is in a duplicated region.
    'is_reference_duplicated',
    # True if all sample calls match.
    'is_consensus',
    # Count of sample calls that differ from the reference.
    'any_snps',
    # True if all sample calls differ from the reference and passed all the filters (coverage, proportion, consensus).
    'is_best_snp',
    'called_reference',
    'num_A',
    'num_C',
    'num_G',
    'num_T',
    'num_N',
    'CallWasMade',
    'PassedDepthFilter',
    'PassedProportionFilter',
    'Pattern'
])


# SampleInfo = namedtuple('SampleInfo', [
# # The base identified at the current contig position
#     'call',
#     'simple_call'
#     #
#     'coverage',
#     #
#     'proportion',
#     # TODO: add quality_breadth
# ])
class SampleInfo(namedtuple('SampleInfo', ['call', 'simple_call', 'coverage', 'proportion'])):
    """
    Attributes:
        call (str): Nucleotide call at the current position
        is_called (bool): True if call != X (no call) or N (any call)
        is_snp (bool): True if call != reference call
        coverage (float):
        proportion (float):
    """
    __slots__ = ()

    def __new__(cls, call='X', simple_call='X', **sample_info):
        """
        Set is_called based on the value of call
        """
        return super().__new__(cls, call=call, simple_call=cls.simple_call(call), **sample_info)

    # NOTE: There are a lot of samples, does attaching this method add a lot of accumulative weight?
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


class Sample(metaclass=ABCMeta):
    def __init__(self, filepath, name, aligner, snpcaller):
        self._filepath = filepath
        self._name = name
        self._aligner = aligner
        self._snpcaller = snpcaller
        self._index = self._index_contigs()

    @property
    def identifier(self):
        """
        Returns:
            str: sample_name::aligner,snpcaller
        """
        return "{0}::{1},{2}".format(self._name, self._aligner, self._snpcaller)

    @abstractmethod
    def get_contig(self, contig_name):
        """
        Returns:
            generator: Yields SampleInfo at every position.
        """
        pass

    @abstractmethod
    def _index_contigs(self):
        pass


class Vcf(Sample):
    def __init__(self, filepath, name, aligner, snpcaller):
        super(Vcf, self).__init__(filepath, name, aligner, snpcaller)

    def get_contig(self, contig_name):
        """
        Return a Contig instance for a contig with the given name.

        Args:
            contig_name (str): Name of contig to retrieve.

        Returns:
            VcfContig or EmptyContig: VcfContig if the contig exists, otherwise an EmptyContig.
        """
        file_position = self._index.get(contig_name)
        return VcfContig(contig_name, self._filepath, file_position) if file_position is not None else EmptyContig(contig_name)

    def _index_contigs(self):
        """
        Return:
            dict: A dictionary of contig names and their starting file position.
        """
        index = {}
        file_position = 0
        with open(self._filepath) as handle:
            # Skip metadata and header
            for line in handle:
                if line.startswith('#CHROM'):
                    break
            else:
                # FIXME: use appropriate error
                raise Exception('Mandatory header not found in ' + self._filepath)
            for line in handle:
                contig_name = line[:line.index('\t')]
                if contig_name not in index:
                    index[contig_name] = file_position
                file_position += len(line)
        return index


class Fasta(Sample):

    def __init__(self, filepath, name, aligner):
        super(Fasta, self).__init__(filepath, name, aligner, None)

    def get_contig(self, contig_name):
        """
        Return a Contig instance for a contig with the given name.

        Args:
            contig_name (str): Name of contig to retrieve.

        Returns:
            VcfContig or EmptyContig: VcfContig if the contig exists, otherwise an EmptyContig.
        """
        file_position = self._index.get(contig_name)
        return FastaContig(contig_name, file_position) if file_position else EmptyContig(contig_name)

    def _index_contigs(self):
        """
        Return:
            dict: A dictionary of contig names and their starting file position.
        """
        with open(self._filepath) as handle:
            index = {}
            file_position = 0
            for line in handle:
                if line.startswith('>'):
                    # It is assumed the fasta has the following format:
                    # >contig_name [description]
                    # GATC...
                    # >contig_name [description]
                    # GATC...
                    contig_name = line.partition(' ')[0][1:]
                    if contig_name not in index:
                        index[contig_name] = file_position
                file_position += len(line)
            return index


class Contig(metaclass=ABCMeta):
    @property
    @abstractmethod
    def name(self):
        """
        Returns:
            str: The contig name
        """
        pass

    @property
    @abstractmethod
    def positions(self):
        """
        Returns:
            generator: Yields SampleInfo at every position.
        """
        pass


class EmptyContig(Contig):
    """
    EmptyContig is a placeholder that yields empty positions for a contig that does not exist in a Sample.
    """
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def positions(self):
        while True:
            yield SampleInfo(
                call='X',
                coverage=None,
                proportion=None,
            )


class FastaContig(Contig):
    """
    TODO: Add support for the .fai index format:
    http://manpages.ubuntu.com/manpages/trusty/man5/faidx.5.html

    FastaContig can be used to read frankenfastas, duplicates.txt, and the reference fasta.
    """

    def __init__(self, name, file_path, file_position=0):
        """
        Attributes:
            _file_path (str): Path to fasta file.
            _file_position (str): File position of the contig description as distinguished by a line starting with '>'
            _name (str): The contig name defined in the fasta file.
        """
        self._file_path = file_path
        self._name = name
        self._file_position = file_position

    @property
    def name(self):
        # """
        # Returns:
        #     str: The contig name.
        # """
        # if self._name is None:
        #     # If undefined, read the contig name from the file.
        #     with open(self._file_path) as handle:
        #         handle.seek(self._file_position)
        #         # TODO: Strip the description after the first whitespace character.
        #         self._name = handle.readline().rstrip()[1:]
        return self._name

    @property
    def positions(self):
        """
        Note:
            The generator assumes the fasta has the following format:

            >name [description]
            GATC...
            >name [description]
            GATC...

        Returns:
            generator: Yields SampleInfo.
        """
        with open(self._file_path) as handle:
            # Seek to the contig file position, read the name, and move to the first contig position.
            handle.seek(self._file_position)
            # TODO: Strip description following the first whitespace?
            self.contig_name = handle.readline().rstrip()[1:]
            for line in handle:
                # Discard newline character
                line = line.rstrip()
                # Stop iterating at EOF, blank line, or the start of the next contig
                if line == "" or line.startswith('>'):
                    return
                for call in line:
                    yield SampleInfo(
                        call=call,
                        coverage='-',
                        proportion='-',
                    )


# FIXME: Support for multiple samples in a vcf is disabled.
class VcfContig(Contig):
    """
    See http://samtools.github.io/hts-specs/VCFv4.2.pdf
    """

    def __init__(self, name, num_samples, file_path, file_position=0):
        """
        Args:
            file_path (str): Path to the VCF file.
            file_positions (int): Seek offset from the beginning of the file to the first position of this contig.
            sample_names (tuple):
            index (int):
        """
        self._name = name
        self._num_samples = num_samples
        self._file_path = file_path
        self._file_position = file_position

    def _get_coverage(self, record):
        """
        Args:
            record (dict):

        Returns:
            int or 'PASS' or None:
        """
        sample_coverage = None
        if record[self._name].get('DP') and record[self._name]['DP'].isdigit():
            sample_coverage = int(record[self._name]['DP'])
        elif record['INFO'].get('DP') and record['INFO']['DP'].isdigit():
            sample_coverage = int(record['INFO']['DP']) / self._num_samples
        elif record['INFO'].get('ADP') and record['INFO']['ADP'].isdigit():
            sample_coverage = int(record['INFO']['ADP']) / self._num_samples
        # NASP output
        elif 'FT' in record[self._name]:
            failed_filters = record[self._name]['FT'].split(',')
            if 'CovFail' in failed_filters:
                sample_coverage = -1
            elif 'PASS' in failed_filters or 'PropFail' in failed_filters:
                sample_coverage = 'PASS'
        return sample_coverage

    def _get_proportion(self, record, sample_coverage, is_snp):
        """
        Args:
            record (dict):
            sample_name (str):
            sample_coverage (int):
            is_snp (bool):

        Returns:
            int or float or 'PASS' or None:
        """
        sample_proportion = None
        if 'AD' in record[self._name]:
            call_depths = record[self._name]['AD'].split(',')
            # gatk, reliable and documented
            if len(call_depths) > 1:
                # The genotype (GT) may be phased (|), unphased (/).
                alt_number = record[self._name]['GT'].split('|', 1)[0].split('/', 1)[0]
                if alt_number.isdigit():
                    sample_proportion = int(call_depths[int(alt_number)]) / sample_coverage
                    # varscan, reliable and documented
            elif is_snp:
                sample_proportion = int(call_depths[0]) / sample_coverage
            elif not is_snp and 'RD' in record[self._name]:
                sample_proportion = int(record[self._name]['RD']) / sample_coverage
        # solsnp, undocumented, no multi-sample support
        elif 'AR' in record['INFO']:
            sample_proportion = float(record['INFO']['AR'])
            if not is_snp:
                sample_proportion = 1 - sample_proportion
        # samtools, estimate, dubious accuracy
        elif 'DP4' in record['INFO']:
            call_depths = record['INFO']['DP4'].split(',')
            if is_snp:
                sample_proportion = ( int(call_depths[2]) + int(call_depths[3]) ) / (
                    sample_coverage * self._num_samples)
            else:
                sample_proportion = ( int(call_depths[0]) + int(call_depths[1]) ) / (
                    sample_coverage * self._num_samples )
        # NASP output
        elif 'FT' in record[self._name]:
            failed_filters = record[self._name]['FT'].split(',')
            if 'PropFail' in failed_filters:
                sample_proportion = -1
            elif 'PASS' in failed_filters:
                sample_proportion = 'PASS'
        return sample_proportion

    @property
    def name(self):
        """
        Returns:
            str: Contig name read from the header column this instance represents.
        """
        return self._name

    @property
    def positions(self):
        """
        Note:
            It is assumed the reference column in the VCF matches the analysis reference.
            The FORMAT column is defined with a GT attribute.

        Returns:
            generator: Yields SampleInfo.

        Examples:
            If there are gaps between positions, the generator returns a result equivalent to an empty/uncalled position.

            .. csv-table:: Position Gaps
                :header: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,example_1_L001

                "500WT1_test","320",".","GTT","G","900848.97",".","AC=1;AF=1.00;AN=1;BaseQRankSum=1.732;DP=18490;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=59.72;MQ0=0;MQRankSum=0.013;QD=24.36;ReadPosRankSum=0.538","GT:AD:DP:GQ:MLPSAC:MLPSAF:PL","1:1,16410:18355:99:1:1.00:900888,0"
                "500WT1_test","323",".","G",".","676339",".","AN=1;DP=18212;MQ=59.72;MQ0=9","GT:DP:MLPSAC:MLPSAF","0:18182"

            >>> SampleInfo(
            >>>    call='X',
            >>>    coverage=None,
            >>>    proportion=None,
            >>> )

            If there are duplicate positions from indels, the generators skips past them.

            .. csv-table:: Duplicate Position
                :header: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,example_1_L001

                "500WT1_test","320",".","G",".","714059",".","AN=1;DP=18490;MQ=59.72;MQ0=10","GT:DP:MLPSAC:MLPSAF","0:18423"
                "500WT1_test","320",".","GTT","G","900848.97",".","AC=1;AF=1.00;AN=1;BaseQRankSum=1.732;DP=18490;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=59.72;MQ0=0;MQRankSum=0.013;QD=24.36;ReadPosRankSum=0.538","GT:AD:DP:GQ:MLPSAC:MLPSAF:PL","1:1,16410:18355:99:1:1.00:900888,0"
        """
        with open(self._file_path, 'r') as handle:
            position = 1
            fieldnames = []
            # Skip past the metadata to the header
            for line in handle:
                if line.startswith('#CHROM'):
                    # Fieldnames are parsed here instead of by the DictReader because the file pointer has already moved
                    # passed them.
                    fieldnames = line.strip().split('\t')
                    break

            for row in csv.DictReader(handle, fieldnames=fieldnames, delimiter='\t'):
                # Stop iterating if a new contig has started.
                if row['#CHROM'] != self._name:
                    return

                row['POS'] = int(row['POS'])

                # TODO: Collect indel statistics instead of skipping them.
                # Skip past duplicate positions from indels.
                if row['POS'] < position:
                    continue

                # TODO: what should the value of coverage and proportion be? Should it be a question mark? A negative number?
                # Yield empty rows to bridge position gaps.
                while row['POS'] > position:
                    position += 1
                    yield SampleInfo(
                        call='X',
                        coverage=None,
                        proportion=None,
                    )

                # REF is concatenated with ALT because the genotype (GT) allele values are indices for which ALT belongs
                # to which sample with REF at index 0.
                if row['ALT'] == '.':
                    row['ALT'] = [row['REF']]
                else:
                    row['ALT'] = [row['REF']] + row['ALT'].split(',')

                # Parse INFO column such as NS=3;DP=14;AF=0.5;DB;H2 to a dictionary:
                # { 'NS': '3', 'DP': '4', 'AF': '0.5', 'DB': None, 'H2': None }
                # The list comprehension splits the string by the ';' delimiter and each key=value pair by '='.
                # If the element is not a key=value pair, the if condition assigns it the value None to guard against
                # IndexErrors. Finally, the dictionary comprehension assembles the final result.
                row['INFO'] = {elem[0]: elem[1] if len(elem) == 2 else None for elem in
                               [key_value.split('=', 1) for key_value in row['INFO'].split(';')]}

                # Join the FORMAT column vector with the sample column value vector as a dictionary.
                # If the FORMAT column and a sample column named 'NA00001' contained the following vectors:
                #
                # +--------------+-----------------+
                # | FORMAT       |  NA00001        |
                # | ===========  +  ============== |
                # | GT:GQ:DP:HQ  |  0|0:48:1:51,51 |
                # +--------------+-----------------+
                #
                # The row key 'NA00001' would have the following value:
                # { 'GT': 0|0, 'GQ': 48, 'DP': 1, 'HQ': 51,51 }
                row[self._name] = dict(zip(row['FORMAT'].split(':'), row[self._name].split(':')))

                # TODO: Handle '.' and VarScan special case. See VCFRecord._get_sample_call()
                # Use the GT index to lookup the ALT allele for the sample. Strip all but the first character.
                call = row['ALT'][int(row[self._name]['GT'].split('/', 1)[0].split('|', 1)[0])][:1]
                coverage = self._get_coverage(row)
                proportion = self._get_proportion(row, coverage, call != row['REF'])

                # Yield the position.
                position += 1
                yield SampleInfo(
                    # FIXME: call is a raw value from the ALT column
                    call=call,
                    coverage=coverage,
                    proportion=proportion,
                )

                # TODO: Is there a situation where the sample contig would be shorter than the reference?
                # The generator could use a while loop to yield empty rows after the file is exhausted.


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


        # ReferenceCalls = namedtuple("ReferenceCalls", ['reference_call', 'duplicate_call'])
        # {
        #     contig_name: (file_position, length)
        # }


        # class Sample(metaclass=ABCMeta):
        #     @property
        #     @abstractmethod
        #     def records(self):
        #         """
        #         Returns:
        #             generator: Yields contig position information.
        #         """
        #         pass

        # def index_fasta(filepath):
        #     """
        #     Indexes the file positions of all the contigs in the fasta file and returns a Fasta initialized with the indices.
        #
        #     Args:
        #         filepath (str):
        #
        #     Returns:
        #         Fasta:
        #     """
        #     with open(filepath, 'r') as handle:
        #         index = {}
        #         file_position = 0
        #         for line in handle:
        #             file_position += len(line)
        #             if line.startswith('>'):
        #                 index[line[1:].rstrip()] = file_position
        #
        #     print(index)
        #
        #     return Fasta(filepath, index)
        #
        # filepaths = ['/Users/jtravis/results/Washington_CQ3.fasta']
        #
        # cache_size = 80
        #
        # with ProcessPoolExecutor() as executor:
        #     fastas = [fasta for fasta in executor.map(index_fasta, filepaths)]
        #
        #     for record in fastas[0].records:
        #         print(record)

        # class ReferenceContig(FastaContig):
        #
        #     def __init__(self, ref_path, dups_path=None, file_position=0):
        #         """
        #         Args:
        #             ref_path (str): Path to the reference genome fasta.
        #             dups_path (str): Path to the duplicates file.
        #             file_position (int): File position of the first contig position after the contig description
        #         """
        #         super().__init__(ref_path, file_position)
        #
        #         # with open(self.filepath) as reference_handle, open(self.dups_path) if self.dups_path else none_context as dups_handle:
        #         #     # TODO: handle (multiple) contigs
        #         #     for line in zip(reference_handle, dups_handle):
        #         #         # TODO: discard newline character
        #         #         # Stop iterating at EOF or the start of the next contig
        #         #         if line[0].startswith('>'):
        #         #             return
        #         #         for c in line:
        #         #             yield ReferenceCalls(line[0], line[1])


        # class MatrixContig(Contig):
        #     # TODO: Raise exception if matrix reference != analysis reference
        #     # "The reference in Matrix <filename> does not match the analysis reference"
        #
        #     def __init__(self, file_path, file_position, contig_name):
        #         """
        #         Attributes:
        #             _file_path (str): Path to the NASP matrix file.
        #             _file_position (int): File position of the first contig record.
        #             _filter (list): Samples to filter from the matrix.
        #             _name (str): Contig description read from the #CHROM column.
        #         """
        #         self._file_path = file_path
        #         self._file_position = file_position
        #         self._name = None
        #
        #     @property
        #     def name(self):
        #         """
        #         Returns:
        #             str: Contig description read from the #CHROM column.
        #         """
        #         if self._name is None:
        #             with open(self._file_path) as handle:
        #                 handle.seek(self._file_position)
        #                 # TODO: Can the read stop at the first tab instead of the entire line?
        #                 self._name, _, _ = handle.readline().partition('\t')
        #
        #         return self._name
        #
        #     @property
        #     def positions(self):
        #         """
        #         Returns:
        #             dict generator: Yields contig records.
        #         """
        #         with open(self._file_path) as handle:
        #             # Parse the headers from the first line.
        #             headers = handle.readline().rstrip().split('\t')
        #             try:
        #                 # Samples are between the Reference and #SNPcall columns in the Master Matrix
        #                 samples = headers[2:headers.index('#SNPcall')]
        #             except ValueError:
        #                 # TODO: Raise exception, is this a master matrix?
        #                 pass
        #
        #             # Move file pointer to the contig
        #             handle.seek(self._file_position)
        #             for row in csv.DictReader(handle, fieldnames=headers, delimiter='\t'):
        #                 # Stop iterating at EOF or the start of the next contig
        #                 if row['#CHROM'] != self._name:
        #                     return
        #                 # TODO: Create a selectors table representing the filtered samples and use the itertools.compress()
        #                 # function to filter the position strings: CallWasMade, PassedDepthFilter, PassedProportionFilter,
        #                 # Pattern. How can the other fields, such as the passed filter columns be updated?
        #                 yield SampleInfo(
        #                     call=row[self.sample_name],
        #                     coverage=row['PassedDepthFilter'][self.i],
        #                     proportion=row['PassedProportionFilter'][self.i],
        #                 )