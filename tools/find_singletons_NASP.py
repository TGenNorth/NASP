#!/usr/bin/env python

"""find autapomorphic SNPs for a single isolate
from a SNP matrix created by NASP"""

from optparse import OptionParser
from nasp.util import get_field_index
from nasp.util import get_genome_index
from nasp.util import get_singleton_snps
from nasp.util import inverse_filter_matrix

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()
        
def main(matrix_in, genome_name):
    last=get_field_index(matrix_in)
    my_genome=get_genome_index(matrix_in, genome_name)
    get_singleton_snps(matrix_in, my_genome, last)
    """this is used to filter unique SNPs out of a NASP matrix"""
    inverse_filter_matrix(matrix_in, "filtered", "unique_snps.txt")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-m", "--input_matrix", dest="matrix_in",
                    help="/path/to/isg_matrix [REQUIRED]",
                    action="callback", callback=test_file, type="string")
    parser.add_option("-g", "--genome_name", dest="genome_name",
                    help="genome to look for singleton SNPs [REQUIRED]",
                    action="store", type="string")

    options, args = parser.parse_args()
    
    mandatories = ["matrix_in", "genome_name"]
    for m in mandatories:
	if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
	    parser.print_help()
	    exit(-1)
            
    main(options.matrix_in, options.genome_name)
