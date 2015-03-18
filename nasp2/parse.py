"""
The parse module reads Fasta and Vcf files created using combinations of aligners and snpcallers.

To minimize memory usage the files are read line by line.

To compare the files, the formats are abstracted with an interface that creates the illusion that all files have all
have corresponding contigs, positions, and data formats.

The key to the interface is that the iterators yield information for each position even when there is no data and all
the positions have the same format.


"""

__author__ = 'jtravis'

from abc import ABCMeta, abstractmethod
from collections import namedtuple, OrderedDict
import csv
import logging
logging.basicConfig(filename='parse.log', level=logging.DEBUG)


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
class Position(namedtuple('SampleInfo', ['call', 'simple_call', 'coverage', 'proportion'])):
    """
    SampleInfo is all the data collected for a single Sample position.

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
        Set simple_call based on the value of call
        """
        return super().__new__(cls, call=call, simple_call=cls._simple_call(call), **sample_info)

    # TODO: Remove allow_x and allow_del parameters
    # NOTE: There are a lot of samples, does attaching this method add a lot of accumulative weight?
    @staticmethod
    def _simple_call(dna_string, allow_x=False, allow_del=False):
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


class SampleAnalysis(metaclass=ABCMeta):

    def __init__(self, filepath, name, aligner, snpcaller):
        """
        A SampleAnalysis represents a file composed of DNA Contigs for a specific aligner and snpcaller combination.
        Multiple SampleAnalyses may be run on the same sample.

        Args:
            filepath (str): Path to the SampleAnalysis file.
            name (str): Name of the sample.
            aligner (str): Name of the aligner program used to create this analysis.
            snpcaller (str): Name of the snpcaller program used to create this analysis.
        """
        self._filepath = filepath
        # TODO: Remove rsplit if the name property in the matrix_dto.xml no longer includes the suffix.
        # Strip the '-aligner-snpcaller' suffix if it exists.
        self._name = name.rsplit('-', 2)[0]
        self._aligner = aligner
        self._snpcaller = snpcaller
        self._index = self._index_contigs()
        logging.debug("{0}._index => {1!r}".format(self.__repr__(), self._index))

    def __lt__(self, other):
        """
        __lt__ is used to sort SampleAnalyses before grouping by sample name.

        Args:
            other (SampleAnalysis):

        Return
            bool: True if the name of the sample this analysis represents is alphanumerically less than the `other`.
        """
        # return self._name < other.name
        return self.identifier < other.identifier

    @property
    def name(self):
        """
        Return
            str: The sample name which may be shared by many SampleAnalyses.
        """
        return self._name

    @property
    def identifier(self):
        """
        A label for the user to identify a SampleAnalysis amidst others for the same sample.

        Note:
            While it should not happen, there are no guards against two SampleAnalyses having the same identifier.

        Returns:
            str: sample_name::aligner,snpcaller
        """
        return "{0}::{1},{2}".format(self._name, self._aligner, self._snpcaller)

    # @property
    # @abstractmethod
    # def contigs(self):
    #     """
    #     The contigs property is used to iterate over the reference contigs.
    #     It is not needed for any other SampleAnalysis.
    #
    #     Return:
    #         Contig generator: All the contigs in the SampleAnalysis
    #     """
    #     pass

    @abstractmethod
    def get_contig(self, contig_name):
        """
        Return a Contig instance for a contig with the given name.

        Args:
            name (str): Name of contig to retrieve.

        Returns:
            Contig or EmptyContig: Contig if the contig exists, otherwise an EmptyContig.
        """
        pass

    @abstractmethod
    def _index_contigs(self):
        """
        Build an index of file positions for all the contigs in the SampleAnalysis. The file is scanned once so that the
        Contigs returned by `get_contig` can seek straight to the file position of their base.
        """
        pass

    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self._filepath, self._name, self._aligner, self._snpcaller)


class Vcf(SampleAnalysis):

    def __init__(self, filepath, name, aligner, snpcaller):
        """
        A Vcf represents a Variant Call Format file.

        TODO: is name and _sample_name the same thing?

        Args:
            filepath (str): Path to Vcf file.
            name (str): Name of the sample.
            aligner (str): Name of the aligner used to create this file.
            snpcaller (str): Name of the snpcaller used to create this file.

        Attributes:
            _sample_name: Name of the sample as read from the header following the FORMAT column.
        """
        # _sample_name is declared before super because it is initialized when super calls _index_contigs.
        # It is declared here to avoid a hidden attribute.
        self._sample_name = None
        super(Vcf, self).__init__(filepath, name, aligner, snpcaller)

    # TODO: Disable/remove?
    @property
    def contigs(self):
        """
        Return:
            Contig generator: All the contigs in the SampleAnalysis
        """
        return (VcfContig(contig_name, self._sample_name, self._filepath, file_position) for contig_name, file_position in self._index.items())

    def get_contig(self, contig_name):
        """
        Return a Contig instance for a contig with the given name.

        Args:
            contig_name (str): Name of contig to retrieve.

        Returns:
            VcfContig or EmptyContig: VcfContig if the contig exists, otherwise an EmptyContig.
        """
        file_position = self._index.get(contig_name)
        if file_position is None:
            logging.info("{0}.get_contig({1!r}) => EmptyContig({1!r})".format(self.__repr__(), contig_name))
            return EmptyContig(contig_name)
        return VcfContig(contig_name, self._sample_name, self._filepath, file_position)

    def _index_contigs(self):
        """
        Note:
            If the file contains multiple contigs with the same name, all but the first will be ignored.
            A multiple contigs error is not detected.

        Return:
            dict: A dictionary of contig names and their starting file position.
        """
        # An OrderedDict is not required, it is used here to be consistent with the Fasta class.
        index = OrderedDict()
        file_position = 0
        with open(self._filepath) as handle:
            # Skip metadata and header
            for line in handle:
                file_position += len(line)
                if line.startswith('#CHROM'):
                    # FIXME: Add support for multiple samples
                    # The sample columns follow the FORMAT column. Here is is assumed there is only one sample.
                    self._sample_name = line.rstrip()[line.index('FORMAT')+len('FORMAT\t'):].split('\t', 1)[0]
                    break
            else:
                # FIXME: use appropriate error
                raise Exception('Mandatory header not found in Vcf file:' + self._filepath)
            for line in handle:
                # The first column, #CHROM, is the contig name.
                contig_name = line[:line.index('\t')]
                if contig_name not in index:
                    index[contig_name] = file_position
                file_position += len(line)
        return index


ContigIndex = namedtuple('ContigIndex', ['length', 'file_position'])


class Fasta(SampleAnalysis):

    def __init__(self, filepath, name, aligner, is_reference=False):
        """
        A Fasta represents a Fasta file.

        Args:
            filepath (str): Path to Vcf file.
            name (str): Name of the sample.
            aligner (str): Name of the aligner used to create this file.
            is_reference (bool): True if the file is the reference genome. All other files will yield empty when the contig is exhausted.
        """
        super(Fasta, self).__init__(filepath, name, aligner, None)
        self._is_reference = is_reference

    @property
    def contigs(self):
        """
        Return
            FastaContig generator: Yields contigs sorted by length, longest to shortest.
        """
        # return (FastaContig(name, contig_index.length, contig_index.file_position, self._filepath, self._is_reference) for name, contig_index in self._index.items())
        return (FastaContig(name, contig_index.length, contig_index.file_position, self._filepath, self._is_reference) for name, contig_index in sorted(self._index.items(), key=lambda x: x[1].length, reverse=True))

    def get_contig(self, contig_name):
        """
        Return a Contig instance that can iterate over the positions for a contig with the given name.

        Args:
            contig_name (str): Name of contig to retrieve.

        Returns:
            FastaContig or EmptyContig: FastaContig if the contig exists, otherwise an EmptyContig.
        """
        contig_index = self._index.get(contig_name)
        if contig_index is None:
            logging.info("{0}.get_contig({1!r}) => EmptyContig({1!r})".format(self.__repr__(), contig_name))
            return EmptyContig(contig_name, is_fasta=True)
        return FastaContig(contig_name, contig_index.length, contig_index.file_position, self._filepath, self._is_reference)

    def _index_contigs(self):
        """
        It is assumed the fasta has the following format:

            >contig_name [description]\n
            GATC...
            >contig_name [description]\n
            GATC...

        There may be blank lines between contigs and a contig sequence may only be broken by newline characters.

        Note:
            If the file contains multiple contigs with the same name, all but the last will be ignored.
            A multiple contigs error is not detected.
            TODO: Raise DuplicateContigNameError

        Return:
            dict: A dictionary of contig names and their starting file position.
        """
        with open(self._filepath) as handle:
            # OrderedDict is used so that the reference can iterate over its Contigs in the encountered file order.
            index = OrderedDict()

            # The entire contig block must be scanned in order to get both the file position and length.
            # A contig name of None indicates this is the first contig block. There is no previous block to index.
            contig_name = None
            contig_len = 0
            # contig_file_position marks the first base of the contig_name.
            contig_file_position = 0
            # file_position is a running counter for the file pointer since ftell() is disabled by the file iterator.
            file_position = 0

            for line in handle:
                file_position += len(line)
                if line.startswith('>'):
                    # TODO: Raise exception. There should never be two contigs with the same name.
                    # if contig_name in index:
                    #     raise Exception()

                    # Starting a new Contig block. If there was a previous block we should be able to index it having
                    # both the file position and length.
                    if contig_name is not None:
                        index[contig_name] = ContigIndex(length=contig_len, file_position=contig_file_position)

                    # contig_name is everything between the '>' and the first whitespace character.
                    # As a result it excludes the optional description.
                    contig_name = line.partition(' ')[0][1:].rstrip()
                    # TODO: This if statement removes a prefix added to the franken fasta contigs. Why is this prefix needed? Remove both.
                    if contig_name.startswith('franken::'):
                        contig_name = contig_name[len('franken::'):]
                    contig_file_position = file_position
                    contig_len = 0
                else:
                    contig_len += len(line.rstrip())

            # The end of the last contig block is marked by EOF instead of a new contig block so it is added to
            # the index explicitly.
            # TODO: Raise exception. There should never be two contigs with the same name.
            # if contig_name in index:
            #     raise Exception()
            if contig_name is not None:
                index[contig_name] = ContigIndex(length=contig_len, file_position=contig_file_position)
            return index


class Contig(metaclass=ABCMeta):
    # An empty position represents gaps between positions or the sample contig is shorter than the
    # reference contig.
    FASTA_EMPTY_POSITION = Position(
        call='X',
        coverage='-',
        proportion='-'
    )

    VCF_EMPTY_POSITION = Position(call='X', coverage='?', proportion='?')

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

    def __repr__(self):
        return "{0}({1!r})".format(self.__class__.__name__, self.name)


class EmptyContig(Contig):
    def __init__(self, name, is_fasta=False):
        """
        EmptyContig is a placeholder that yields empty positions for a contig that does not exist in a SampleAnalysis.

        Args:
            name (str): Name of the contig.
            is_fasta (bool): A Fasta position always passes coverage and proportion filters.
        """
        self._name = name
        self._is_fasta = is_fasta

    @property
    def name(self):
        return self._name

    @property
    def positions(self):
        if self._is_fasta:
            while True:
                yield self.FASTA_EMPTY_POSITION
        else:
            while True:
                yield self.VCF_EMPTY_POSITION

    def __repr__(self):
        return "{0}(name={1!r}, is_fasta={2!r})".format(self.__class__.__name__, self.name, self.is_fasta)


class FastaContig(Contig):

    def __init__(self, name, length, file_position, file_path, is_reference):
        """
        FastaContig can be used to read frankenfastas, duplicates.txt, and the reference fasta.

        Attributes:
            _name (str): The contig name defined in the fasta file.
            _length (int): Length of the contig sequence. It is used to
            _file_path (str): Path to fasta file.
            _file_position (str): File position of the first nucleotide base following the contig description.
        """
        self._name = name
        self._length = length
        self._file_path = file_path
        self._file_position = file_position
        self._is_reference = is_reference

    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self.name, self._length, self._file_position, self._file_path)

    def __len__(self):
        return self._length

    @property
    def name(self):
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
            generator: Yields SampleInfo. Since a Fasta file has no coverage or proportion information it is assumed
             all positions pass coverage and proportion.
        """
        with open(self._file_path) as handle:
            # Seek to the file position of the first base of this contig.
            handle.seek(self._file_position)
            for line in handle:
                # Discard newline character
                line = line.rstrip()
                # Stop iterating at EOF, blank line, or the start of the next contig
                if line == "" or line.startswith('>'):
                    return
                for call in line:
                    # FIXME: Will having "infinite" values create problems later on?
                    yield Position(
                        call=call,
                        #coverage=float('Infinity'),
                        #proportion=float('Infinity'),
                        coverage='-',
                        proportion='-',
                    )

        # FIXME: If the FastaContig is the reference, this while loop creates an infinite loop.
        # Yield empty positions when the contig is exhausted.
        while not self._is_reference:
            yield self.FASTA_EMPTY_POSITION


# FIXME: Support for multiple samples in a vcf is disabled.
class VcfContig(Contig):
    """
    See http://samtools.github.io/hts-specs/VCFv4.2.pdf
    """

    def __init__(self, contig_name, sample_name, filepath, file_position):
        """
        Args:
            contig_name (str): Contig name from the #CHROM column
            sample_name (str): Sample name header for the sample this contig represents.
            filepath (str): Path to the VCF file.
            file_position (int): Seek offset from the beginning of the file to the first position of this contig.
        """
        self._name = contig_name
        self._sample_name = sample_name
        # FIXME: Support for multiple samples in a vcf is disabled. _num_samples should be the number of samples in the VCF.
        self._num_samples = 1
        self._filepath = filepath
        self._file_position = file_position

    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self._name, self._sample_name, self._filepath, self._file_position)

    def _get_sample_call(self, record):
        """
        Args:
            record (dict):

        Note:
            Has the side effect of modifying the `record` ALT key.
        """

        # REF is concatenated with ALT because the genotype (GT) allele values are indices for which ALT belongs
        # to which sample with REF at index 0.
        if record['ALT'] == '.':
            record['ALT'] = [record['REF']]
        else:
            record['ALT'] = [record['REF']] + record['ALT'].split(',')

        # FIXME indels
        return_value = None
        if len(record['ALT']) == 1:
            return_value = record['ALT'][0]
        if 'GT' in record[self._sample_name] and record[self._sample_name]['GT'] != '.':
            return_value = record['ALT'][int(record[self._sample_name]['GT'])]
            # OMG varscan
            # Handles a scenerio such as this:
            # +-----+-----+
            # | REF | ALT |
            # +-----+-----+
            # | GTT | TT  |
            # +-----+-----+
            # TODO: clarify what the if statement is doing for parsing varscan output.
            if len(record['REF']) > 1 and (len(record['REF']) - 1 ) == len(return_value) and \
                            record['REF'][:len(return_value)] != return_value and \
                            record['REF'][-len(return_value):] == return_value:
                return_value = record['ALT'][0]

        if return_value is not None:
            return return_value[:1]

        return return_value

    def _get_coverage(self, record):
        """
        Args:
            record (dict):

        Returns:
            int or 'PASS' or None:
        """
        sample_coverage = '-'
        if record[self._sample_name].get('DP') and record[self._sample_name]['DP'].isdigit():
            sample_coverage = int(record[self._sample_name]['DP'])
        elif record['INFO'].get('DP') and record['INFO']['DP'].isdigit():
            sample_coverage = int(record['INFO']['DP']) / self._num_samples
        elif record['INFO'].get('ADP') and record['INFO']['ADP'].isdigit():
            sample_coverage = int(record['INFO']['ADP']) / self._num_samples
        # TODO: Deprecated? The coverage and proportion number was not included in the VCF's produced by NASP because it
        # was thrown away after the file was read. The information is still available with the new parser so the actual
        # value may be used instead.
        # NASP output
        elif 'FT' in record[self._sample_name]:
            failed_filters = record[self._sample_name]['FT'].split(',')
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
            float or 'PASS' or None:
        """
        sample_proportion = '-'
        if 'AD' in record[self._sample_name]:
            call_depths = record[self._sample_name]['AD'].split(',')
            # gatk, reliable and documented
            if len(call_depths) > 1 and record[self._sample_name]['GT'] != '.':
                sample_proportion = int(call_depths[int(record[self._sample_name]['GT'])]) / sample_coverage
                # varscan, reliable and documented
            elif is_snp:
                sample_proportion = int(call_depths[0]) / sample_coverage
            elif not is_snp and 'RD' in record[self._sample_name]:
                sample_proportion = int(record[self._sample_name]['RD']) / sample_coverage
        # solsnp, undocumented, no multi-sample support
        elif 'AR' in record['INFO']:
            sample_proportion = float(record['INFO']['AR'])
            if not is_snp:
                sample_proportion = 1 - sample_proportion
        # samtools, estimate, dubious accuracy
        elif 'DP4' in record['INFO']:
            call_depths = record['INFO']['DP4'].split(',')
            if is_snp:
                sample_proportion = (int(call_depths[2]) + int(call_depths[3])) / (sample_coverage * self._num_samples)
            else:
                sample_proportion = (int(call_depths[0]) + int(call_depths[1])) / (sample_coverage * self._num_samples)
        # TODO: Deprecated? The coverage and proportion number was not included in the VCF's produced by NASP because it
        # was thrown away after the file was read. The information is still available with the new parser so the actual
        # value may be used instead.
        # NASP output
        elif 'FT' in record[self._sample_name]:
            failed_filters = record[self._sample_name]['FT'].split(',')
            if 'PropFail' in failed_filters:
                sample_proportion = -1
            elif 'PASS' in failed_filters:
                sample_proportion = 'PASS'

        # Some big SNP callers, like GATK, do not provide proportion information when
        # the position is called reference.  We cannot filter these positions.
        # TODO: Remove - unnecessary
        if sample_proportion is None and is_snp:
            sample_proportion = '-'

        return sample_proportion

    @property
    def name(self):
        """
        Returns:
            str: Contig name read from the #CHROM column.
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
            If there are gaps between positions or the sample contig is shorter than the reference contig, the generator
            yields values representing empty positions.

            .. csv-table:: Position Gaps
                :header: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,example_1_L001

                "500WT1_test","320",".","GTT","G","900848.97",".","AC=1;AF=1.00;AN=1;BaseQRankSum=1.732;DP=18490;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=59.72;MQ0=0;MQRankSum=0.013;QD=24.36;ReadPosRankSum=0.538","GT:AD:DP:GQ:MLPSAC:MLPSAF:PL","1:1,16410:18355:99:1:1.00:900888,0"
                "500WT1_test","323",".","G",".","676339",".","AN=1;DP=18212;MQ=59.72;MQ0=9","GT:DP:MLPSAC:MLPSAF","0:18182"

            >>> Position(
            >>>    call='X',
            >>>    coverage='?',
            >>>    proportion='?',
            >>> )

            If there are duplicate positions from indels, the generators skips past them.

            .. csv-table:: Duplicate Position
                :header: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,example_1_L001

                "500WT1_test","320",".","G",".","714059",".","AN=1;DP=18490;MQ=59.72;MQ0=10","GT:DP:MLPSAC:MLPSAF","0:18423"
                "500WT1_test","320",".","GTT","G","900848.97",".","AC=1;AF=1.00;AN=1;BaseQRankSum=1.732;DP=18490;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=59.72;MQ0=0;MQRankSum=0.013;QD=24.36;ReadPosRankSum=0.538","GT:AD:DP:GQ:MLPSAC:MLPSAF:PL","1:1,16410:18355:99:1:1.00:900888,0"
        """
        with open(self._filepath, 'r') as handle:
            position = 1
            fieldnames = []
            # Skip past the metadata to the header
            for line in handle:
                if line.startswith('#CHROM'):
                    # Fieldnames are parsed here instead of by the DictReader because the file pointer has already moved
                    # passed them.
                    fieldnames = line.strip().split('\t')
                    break
            else:
                # TODO: Raise an appropriate exception for a missing header.
                raise Exception('VCF missing required header.')

            handle.seek(self._file_position)

            for row in csv.DictReader(handle, fieldnames=fieldnames, delimiter='\t'):

                # Stop iterating if a new contig has started.
                if row['#CHROM'] != self._name:
                    # print('\t', position, row['POS'], 'NEWCONTIG', self._filepath, row)
                    break

                row['POS'] = int(row['POS'])

                # TODO: Collect indel statistics instead of skipping them.
                # Skip past duplicate positions from indels.
                if row['POS'] < position:
                    # print('\t', position, row['POS'], 'INDEL', self._filepath, row)
                    continue

                # Yield empty rows to bridge position gaps.
                while row['POS'] > position:
                    position += 1
                    # TODO: Remove
                    # print('\t', position, row['POS'], 'GAP', self._filepath, row)
                    yield self.VCF_EMPTY_POSITION

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
                row[self._sample_name] = dict(zip(row['FORMAT'].split(':'), row[self._sample_name].split(':')))

                # The genotype (GT) is an index to lookup the ALT allele for the sample.
                # It may be phased (|), unphased (/), or neither. If the value is unavailable it is a '.'
                if 'GT' in row[self._sample_name]:
                    row[self._sample_name]['GT'] = row[self._sample_name]['GT'].split('|', 1)[0].split('/', 1)[0]

                call = self._get_sample_call(row)
                coverage = self._get_coverage(row)
                # call != ref checks if the sample is a SNP. The ref may have more than one position if it's an indel.
                proportion = self._get_proportion(row, coverage, call != row['REF'][:1])

                position += 1
                # print('\t', position, row['POS'], 'NORMAL', self._filepath, row)
                # Yield the position.
                yield Position(
                    call=call,
                    coverage=coverage,
                    proportion=proportion,
                )


            # print('\t', position, row['POS'], 'EXHAUSTED', self._filepath, row)
            # If the sample contig is shorter than the reference, yield empty positions after the file is exhausted.
            while True:
                position += 1
                yield self.VCF_EMPTY_POSITION


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
