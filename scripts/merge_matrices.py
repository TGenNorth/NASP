__author__ = 'jtravis'

from contextlib import ExitStack
import csv
from collections import Counter
import itertools


# TODO: Add a commandline interface.
def merge(outfile, filenames):
    """
    Merge NASP master matrices into a single matrix.

    Note:
        It is assumed the matrices have corresponding rows and use the same reference.
        Each matrix must contain at least one sample column.

    Args:
        outfilename (str): Name of the file where the merged matrix will be written.
        filenames (list): List of NASP Master Matrix files to merge.

    Raises:
        ValueError: The #SNPcall header was not found in one of the matrices.
        DuplicateSampleError: The same sample name was found in the header of one or more matrices. A merged matrix
        would count the values twice.
    """
    # All opened files will automatically be closed at the end of the with statement.
    with ExitStack() as stack:
        # Open all the files as TSV DictReaders.
        matrices = list(map(lambda x: csv.DictReader(x, delimiter='\t'), (stack.enter_context(open(fname))
                                                                          for fname in filenames)))

        # Get a list of sample columns from the header of each matrix so we know which sample keys exist when iterating
        # through each matrix row.
        sample_names = []
        for matrix in matrices:
            fieldnames = matrix.fieldnames
            # FIXME: A ValueError will be raised if #SNPcall is not found which implies an invalid matrix.
            # The sample columns are between the Reference and #SNPcall columns
            sample_names.append(fieldnames[2:fieldnames.index('#SNPcall')])
        # The flattened version is used to get a sample count and build the merged matrix header.
        flattened_sample_names = tuple(itertools.chain.from_iterable(sample_names))
        num_samples = len(flattened_sample_names)

        if len(set(flattened_sample_names)) != num_samples:
            raise DuplicateSampleError('Filtering duplicate samples is not implemented.')

        # Merged Master Matrix headers
        header = ['LocusID', 'Reference']
        header.extend(flattened_sample_names)
        header.extend(('#SNPcall', '#Indelcall', '#Refcall', '#CallWasMade', '#PassedDepthFilter',
                       '#PassedProportionFilter', '#A', '#C', '#G', '#T', '#Indel', '#NXdegen', 'Contig',
                       'Position', 'InDupRegion', 'SampleConsensus', 'CallWasMade', 'PassedDepthFilter',
                       'PassedProportionFilter', 'Pattern', 'Pattern#'))

        writer = csv.DictWriter(stack.enter_context(open(outfile, 'w')), fieldnames=header, delimiter='\t')
        writer.writeheader()

        # int_columns are all the counter columns that will be summed.
        int_columns = ('#SNPcall', '#Indelcall', '#Refcall', '#A', '#C', '#G', '#T', '#Indel', '#NXdegen')
        # str_columns are all the string pattern columns that will be concatenated.
        str_columns = ('CallWasMade', 'PassedDepthFilter', 'PassedProportionFilter', 'Pattern')

        # count_template is a counter initialized to sum and concatenate columns. It is copied at the start of each row
        # instead of re-initializing a new one or resetting the old one which would clear the str_columns.
        # str_columns must be initialized with the empty string in order for Counter to concatenate them as strings.
        # int_columns will default to zero.
        count_template = Counter({k: '' for k in str_columns})

        # pattern_index is an incrementing index of each unique pattern encountered. Since None is an invalid pattern,
        # it is used to keep track of the next available index number in the series.
        pattern_index = {None: 1}

        try:
            # Each iteration will write a merged row until a StopIteration exception is raised at EOF.
            while True:
                # row will contain all the data necessary to write a merged matrix row as key-value pairs.
                row = None
                count = count_template.copy()
                # Merge the columns from each matrix.
                for index, matrix in enumerate(matrices):
                    # Read a dictionary record from each matrix. The next() function is used instead of a loop
                    # in order to read across all the matrices before advancing to the next row.
                    record = next(matrix)
                    # On the first iteration, initialize the row with columns that should be the same for all matrices.
                    if row is None:
                        row = {k: record[k] for k in
                               ('LocusID', 'Reference', 'Contig', 'Position', 'InDupRegion', 'SampleConsensus')}
                    # On all following iterations, if any matrix does not have consensus or the first sample of this
                    # matrix does not match the first sample of the previous matrix, then the merged matrix does not
                    # have consensus. It is assumed all matrices have at least one sample column.
                    elif not record['SampleConsensus'] or \
                            record[sample_names[index][0]] != record[sample_names[index - 1][0]]:
                        row['SampleConsensus'] = False

                    # On all following iterations, strip the '1' representing the reference from the Pattern string.
                    if row is not None:
                        record['Pattern'] = record['Pattern'][1:]
                    # Sum the int_columns and concatenate the str_columns.
                    count.update({k: int(record[k]) for k in int_columns})
                    count.update({k: record[k] for k in str_columns})

                    # Add sample calls to the row.
                    # FIXME: What if two matrices have the same sample? The last one will overwrite the previous,
                    # but the counts and string patterns will be off. Use itertools.compress to filter duplicates?
                    row.update({k: record[k] for k in sample_names[index]})

                # Add the str and int columns to the row.
                row.update(count)

                # Assign an incrementing index value to each unique Pattern encountered.
                if count['Pattern'] not in pattern_index:
                    pattern_index[count['Pattern']] = pattern_index[None]
                    pattern_index[None] += 1
                row['Pattern#'] = pattern_index[count['Pattern']]

                # Each of these is '#/#' as the number of occurrences out of the total number of samples.
                row['#CallWasMade'] = "{0} / {1}".format(num_samples - count['CallWasMade'].count('N'), num_samples)
                row['#PassedDepthFilter'] = "{0} / {1}".format(num_samples - count['PassedDepthFilter'].count('N'),
                                                               num_samples)
                row['#PassedProportionFilter'] = "{0} / {1}".format(
                    num_samples - count['PassedProportionFilter'].count('N'), num_samples)

                writer.writerow(row)
        except StopIteration:
            # Fallthrough: EOF
            pass


class DuplicateSampleError(Exception):
    pass