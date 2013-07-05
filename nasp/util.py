#!/usr/bin/env python

from collections import deque
import collections

def get_field_index(matrix_in):
    """find where the sequence data stops"""
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    return last
    matrix.close()

def get_genome_index(matrix_in, genome_name):
    """find the index of your genome of interest"""
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    my_genome=first_fields.index(genome_name)
    return my_genome
    matrix.close()

def get_singleton_snps(matrix_in, my_genome, last):
    """find those lines in your matrix, where a unique
    SNP state is observed in your genome of interest"""
    matrix = open(matrix_in, "rU")
    firstLine=matrix.readline()
    output = open("unique_snps.txt", "w")
    #print >> output, firstLine,
    for line in matrix:
        hits = [ ]
        fields = line.split("\t")
        for x in fields[1:last]:
            if x in fields[my_genome]:
                hits.append("1")
        if int(len(hits))==1:
            print >> output, fields[0],"\n",

def inverse_filter_matrix(in_matrix, prefix, coords):
    """filter NASP matrix for a user provided list of genome
    coordinates"""
    infile=open(in_matrix, "U")
    coords_file = open(coords, "rU").read().splitlines()
    out_file=open("%s.coord_filtered.matrix" % prefix, "w")
    firstLine = infile.readline()
    sorted(coords_file)
    print >> out_file, firstLine,
    lines = [ ]
    for line in infile:
	fields=line.split("\t")
	if fields[0] not in coords_file:
            lines.append(line)
	    print >> out_file, line,
    print "Number of matching SNPs:", len(lines)

def raw_matrix(matrix_in):
    """bring in ISG matrix and report the total number of SNPs"""
    matrix=open(matrix_in, "U")
    next(matrix)
    lines=[ ]
    for line in matrix:
	    lines.append(line)
    print "Total SNPs:",len(lines)

def matrix_to_fasta(matrix_in, prefix, type, last):
    """converts a NASP matrix to fasta format"""
    reduced = [ ]
    out_fasta = open("%s.%s.fasta" % (prefix, type), "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])

def filter_matrix(matrix_in, last):
    """filter a NASP matrix position if it contains either missing
    or ambiguous data"""
    matrix = open(matrix_in, "rU")
    file_out=open("tmp.matrix", "w")
    firstLine = matrix.readline()
    print >> file_out, firstLine,
    lines = [ ]
    for line in matrix:
	fields = line.split("\t")
        if "X" not in fields[1:last] and "N" not in fields[1:last]: print >> file_out, line,
	if "X" not in fields[1:last] and "N" not in fields[1:last]: lines.append(line)
    print "number of SNPs after filtering:", len(lines)
    matrix.close()

def filter_singletons(last, tmp_matrix):
    """find positions in the matrix where there are at least 2 types of nucleotides
    in the alignment with a minimum frequency of 2 - this is the Mega definition
    of Parsimony informative SNP)"""
    matrix = open(tmp_matrix, "rU")
    file_out=open("tmp2.matrix", "w")
    firstLine = matrix.readline()
    print >> file_out, firstLine,
    lines = [ ]
    for line in matrix:
	fields=line.split("\t")
	counter=collections.Counter(fields[1:last])
        values=counter.values()
        values.sort(key=int)
        if len(values)==1: sys.exc_clear()
        else:
            for i in range(2,5):
                if len(values)==int(i) and values[i-2]>1: print >> file_out, line,
                if len(values)==int(i) and values[i-2]>1: lines.append(line)
        values=[]
    print "number of parsimony-informative SNPs:", len(lines)
    matrix.close()

def filter_genomes(genomes, in_matrix):
    in_matrix = open(in_matrix, "rU")
    firstLine = in_matrix.readline()
    first_fields = firstLine.split()
    last=first_fields.index("#SNPcall")
    all_genomes=first_fields[:last]
    genomes_file = open(genomes, "rU").read().splitlines()
    to_keep = [ ]
    for x in all_genomes[2:]:
        if x in genomes_file:
            to_keep.append(all_genomes.index(x))
    return to_keep
    in_matrix.close()

def filter_matrix_from_genome(to_keep, in_matrix, prefix):
    matrix = open(in_matrix, "rU")
    outfile = open("%s_genomes.matrix" % prefix, "w")
    for line in matrix:
        fields = line.split()
        deque((list.pop(fields, i) for i in sorted(to_keep, reverse=True)), maxlen=0)
        print >> outfile, "\t".join(fields)
    outfile.close()

def filter_matrix_by_coord(in_matrix, prefix, coords):
    """filter NASP matrix for a user provided list of genome
    coordinates"""
    infile=open(in_matrix, "U")
    coords_file = open(coords, "U").read().splitlines()
    out_file=open("%s.coord_filtered.matrix" % prefix, "w")
    firstLine = infile.readline()
    sorted(coords_file)
    print >> out_file, firstLine,
    lines = [ ]
    for line in infile:
	fields=line.split("\t")
	if fields[0] in coords_file:
	   lines.append(line)
	   print >> out_file, line,
    print "Number of matching SNPs:", len(lines)
    infile.close()
    out_file.close()
