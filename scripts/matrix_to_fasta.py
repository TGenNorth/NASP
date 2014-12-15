#!/usr/bin/env python
"""convert matrix to fasta and fasta filtered for "N" and ".", version 2"""

from __future__ import print_function
from optparse import OptionParser
import fileinput
# TODO(jtravis): remove unused subprocess import
import subprocess
import sys
try:
    import collections
except:
    print("python 2.7+ must be installed")
    sys.exit()
import os

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()

def get_field_index(matrix_in):
    """index to find where the SNP calls end"""
    matrix=open(matrix_in, "rU")
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    matrix.close()
    return last

def raw_matrix(matrix_in):
    """bring in ISG matrix and report the total number of SNPs"""
    matrix=iter(fileinput.input([matrix_in]))
    next(matrix)
    lines=[ ]
    for line in matrix:
        lines.append(line)
    print("Total SNPs:",len(lines))
	
def matrix_to_fasta(matrix_in, prefix, type, last):
    """converts an ISG matrix to fasta format"""
    reduced = [ ]
    out_fasta = open("%s.%s.fasta" % (prefix, type), "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:])+"\n")
    out_fasta.close()
   
def filter_matrix(matrix_in, last):
    """filter an ISG matrix position if it contains either missing
    or ambiguous data.  Also removes positions if they are monomorphic"""
    matrix = open(matrix_in, "rU")
    file_out=open("tmp.matrix", "w")
    firstLine = matrix.readline()
    file_out.write(firstLine,)
    lines = [ ]
    for line in matrix:
        doubles=[ ]
        fields = line.split("\t")
        for field in fields[1:last]:
            if len(field)>=2:
	            doubles.append("1")
        if len(doubles)>=1:
            pass
        else:
            counter=collections.Counter(fields[1:last])
            values=counter.values()
            sorted(values, key=int)
            """print line if the calls for a given SNP do not include N or X or is not polymorphic"""
            if "X" not in fields[1:last] and "N" not in fields[1:last] and "." not in fields[1:last] and len(values)>1: file_out.write(line,)
            if "X" not in fields[1:last] and "N" not in fields[1:last] and "." not in fields[1:last] and len(values)>1: lines.append(line)
        doubles=[]
    print("number of SNPs after filtering:", len(lines))
    matrix.close()
    file_out.close()
		
def filter_singletons(last, tmp_matrix):
    """find positions in the matrix where there are at least 2 types of nucleotides
    in the alignment with a minimum frequency of 2 - this is the Mega definition
    of Parsimony informative SNP)"""
    matrix = open(tmp_matrix, "rU")
    file_out=open("tmp2.matrix", "w")
    firstLine = matrix.readline()
    file_out.write(firstLine,)
    lines = [ ]
    for line in matrix:
        fields=line.split("\t")
        counter=collections.Counter(fields[1:last])
        values=counter.values()
        new_values=list(sorted(values, key=int))
        if len(values)==1: sys.exc_clear()
        else:
            for i in range(2,5):
                if len(values)==int(i) and values[i-2]>1: file_out.write(line,)
                if len(values)==int(i) and values[i-2]>1: lines.append(line)
        values=[]
    print("number of parsimony-informative SNPs:", len(lines))
    matrix.close()
    file_out.close()
       
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
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
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.matrix_in, options.prefix)
