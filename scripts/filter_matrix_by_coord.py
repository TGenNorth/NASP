#!/usr/bin/env python

"""filter matrix to only include a list of known SNPs"""

from __future__ import print_function
from optparse import OptionParser
import fileinput
import itertools
import sys
try:
    import collections
except:
    print("you need Python 2.7+ or Python 3+ to run")
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()
   
def get_index_range(in_matrix):
    matrix=open(in_matrix, "rU")
    firstLine = matrix.readline()
    first_fields = firstLine.split("\t")
    try:
        last=first_fields.index("#SNPcall")
    except:
        print("matrix is not readable, is this a valid NASP matrix?")
        sys.exit()
    matrix.close()
    return last

def filter_matrix(in_matrix, prefix, coords):
    """filter NASP matrix for a user provided list of genome
    coordinates"""
    infile=iter(fileinput.input([in_matrix]))
    coords_file = open(coords, "r").read().splitlines()
    out_file=open("%s.coord_filtered.matrix" % prefix, "w")
    firstLine = infile.readline()
    sorted(coords_file)
    out_file.write(firstLine,)
    lines = [ ]
    for line in infile:
        fields=line.split("\t")
        if fields[0] in coords_file:
            lines.append(line)
            out_file.write(line,)
    print("Number of provided SNPs:", len(coords_file))
    print("Number of matching SNPs:", len(lines))
    out_file.close()

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
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:])+"\n")
    return last
    in_matrix.close()
    out_fasta.close()

def filter_missing_matrix(prefix, last):
    """filter a NASP matrix position if it contains either missing
    or ambiguous data"""
    matrix=open("%s.coord_filtered.matrix" % prefix, "rU")
    file_out=open("%s.coord_missing_filtered.matrix" % prefix, "w")
    firstLine = matrix.readline()
    file_out.write(firstLine,)
    lines = [ ]
    for line in matrix:
        fields = line.split("\t")
        counter=collections.Counter(fields[1:last])
        values=counter.values()
        sorted(values,key=int)
    if "X" not in fields[1:last] and "N" not in fields[1:last] and len(values)>1: file_out.write(line,)
    if "X" not in fields[1:last] and "N" not in fields[1:last] and len(values)>1: lines.append(line)
    print("number of filtered SNPs:", len(lines))
    matrix.close()
    file_out.close()
	
def filter_singletons(prefix, matrix_prefix, last):
    """identify the number of autapomorphic SNPs in an ISG matrix"""
    matrix=open("%s.%s" % (prefix, matrix_prefix))
    firstLine = matrix.readline()
    lines = [ ]
    for line in matrix:
        fields=line.split("\t")
        counter=collections.Counter(fields[2:last])
        values=counter.values()
        new_values=list(sorted(values, key=int))
    for i in range(2,5):
	    if len(new_values)==int(i) and new_values[int(i-2)]!=1: lines.append(line)
    print("number of parsimony-informative SNPs:", len(lines))
    matrix.close()

def main(in_matrix, prefix, coords):
    last=get_index_range(in_matrix)
    filter_matrix(in_matrix, prefix, coords)
    matrix_to_fasta(prefix, "coord_filtered", "raw", last)
    filter_missing_matrix(prefix, last)
    matrix_to_fasta(prefix, "coord_missing_filtered", "filtered", last)
    filter_singletons(prefix, "coord_missing_filtered.matrix", last)
	
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-m", "--in_matrix", dest="in_matrix",
                      help="/path/to/file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--out_prefix", dest="prefix",
                      help="/path/to/file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-c", "--in_coords", dest="coords",
                      help="/path/to/file [REQUIRED]",
                      action="callback", callback=test_file, type="string")

    options, args = parser.parse_args()
    
    mandatories = ["in_matrix", "prefix", "coords"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.prefix,options.coords)
