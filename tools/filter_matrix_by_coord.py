#!/usr/bin/env python

"""filter matrix based on a list of known SNPs"""

from optparse import OptionParser
import itertools
import collections

from nasp.util import get_field_index
from nasp.util import filter_matrix_by_coord
from nasp.util import matrix_to_fasta


def matrix_to_fasta(matrix_in, prefix, type, last):
def matrix_to_fasta(prefix, kind, type, last):
    """converts a NASP matrix to fasta format"""
    reduced = [ ]
    in_matrix = open("%s.%s.matrix" % (prefix,kind), "rU")
    out_fasta = open("%s.%s.fasta" % (prefix, type), "w")
    for line in in_matrix:
        fields = line.split("\t")
        reduced.append(fields[2:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
    return last

def filter_missing_matrix(prefix, last):
    """filter a NASP matrix position if it contains either missing
    or ambiguous data"""
    matrix=open("%s.coord_filtered.matrix" % prefix, "rU")
    file_out=open("%s.coord_missing_filtered.matrix" % prefix, "w")
    firstLine = matrix.readline()
    print >> file_out, firstLine,
    lines = [ ]
    for line in matrix:
	fields = line.split("\t")
        if "X" not in fields[1:last] and "N" not in fields[1:last]: print >> file_out, line,
	if "X" not in fields[1:last] and "N" not in fields[1:last]: lines.append(line)
    print "number of filtered SNPs:", len(lines)
	
def filter_singletons(prefix, matrix_prefix, last):
    """identify the number of autapomorphic SNPs in an ISG matrix"""
    matrix=open("%s.%s" % (prefix, matrix_prefix))
    firstLine = matrix.readline()
    next(matrix)
    lines = [ ]
    for line in matrix:
	fields=line.split("\t")
	counter=collections.Counter(fields[2:last])
        values=counter.values()
        values.sort(key=int)
	for i in range(2,5):
		if len(values)==int(i) and values[int(i-2)]!=1: lines.append(line)
    print "number of parsimony-informative SNPs:", len(lines)

def main(in_matrix, prefix, coords):
    last=get_field_index(in_matrix)
    filter_matrix(in_matrix, prefix, coords)
    matrix_to_fasta(prefix, "coord_filtered", "raw", last)
    filter_missing_matrix(prefix, last)
    matrix_to_fasta(prefix, "coord_missing_filtered", "filtered", last)
    filter_singletons(prefix, "coord_missing_filtered.matrix", last)
	
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-i", "--in_matrix", dest="in_matrix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-p", "--out_prefix", dest="prefix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-c", "--in_coords", dest="coords",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")

    options, args = parser.parse_args()
    
    mandatories = ["in_matrix", "prefix", "coords"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.prefix,options.coords)
