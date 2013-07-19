#!/usr/bin/env python
"""convert matrix to fasta and fasta filtered for "N" and ".", version 2"""

from optparse import OptionParser
import os
import sys

from nasp.util import get_field_index
from nasp.util import raw_matrix
from nasp.util import matrix_to_fasta
from nasp.util import filter_matrix
from nasp.util import filter_singletons
	
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix_in, prefix):
    last=get_field_index(matrix_in)
    raw_matrix(matrix_in)
    matrix_to_fasta(matrix_in, prefix, "raw", last)
    filter_matrix(matrix_in, last)
    matrix_to_fasta("tmp.matrix", prefix, "filtered", last)
    filter_singletons(last, "tmp.matrix")
    matrix_to_fasta("tmp2.matrix", prefix, "filtered_PI_snps_only", last)
    os.system("mv tmp.matrix clean_matrix.txt")
    os.system("mv tmp2.matrix clean_PI_matrix.txt")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-m", "--input_matrix", dest="matrix_in",
                      help="/path/to/isg_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--output_prefix", dest="prefix",
                      help="prefix for outfiles [REQUIRED]",
                      action="store", type="string")    

    options, args = parser.parse_args()
    
    mandatories = ["matrix_in", "prefix"]
    for m in mandatories:
	    if not options.__dict__[m]:
		print "\nMust provide %s.\n" %m
		parser.print_help()
		exit(-1)

    main(options.matrix_in, options.prefix)
