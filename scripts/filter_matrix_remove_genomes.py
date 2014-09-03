#!/usr/bin/env python

"""filter a raw ISG matrix for a 
list of genomes"""

from __future__ import print_function
from optparse import OptionParser
import sys
try:
    from collections import deque
except:
    print("must have Python 2.7+ installed or Python 3+")
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()
    
def filter_genomes(genomes, in_matrix):
    in_matrix = open(in_matrix, "rU")
    firstLine = in_matrix.readline()
    first_fields = firstLine.split()
    last=first_fields.index("#SNPcall")
    all_genomes=first_fields[:last]
    genomes_file = open(genomes, "r").read().splitlines()
    to_keep = [ ]
    for x in all_genomes[1:]:
        if x in genomes_file:
            to_keep.append(all_genomes.index(x))
    return to_keep
    in_matrix.close()

def filter_matrix(to_keep, in_matrix, prefix):
    matrix = open(in_matrix, "rU")
    outfile = open("%s_genomes.matrix" % prefix, "w")
    for line in matrix:
        fields = line.split()
        deque((list.pop(fields, i) for i in sorted(to_keep, reverse=True)), maxlen=0)
        outfile.write("\t".join(fields))
    outfile.close()
    matrix.close()
    
def main(in_matrix, prefix, genomes):
    to_keep=filter_genomes(genomes, in_matrix)
    filter_matrix(to_keep, in_matrix, prefix)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-m", "--in_matrix", dest="in_matrix",
                      help="/path/to/file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--out_prefix", dest="prefix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-g", "--genomes", dest="genomes",
                      help="/path/to/genomes_file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    
    mandatories = ["in_matrix", "prefix", "genomes"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.prefix,options.genomes)
