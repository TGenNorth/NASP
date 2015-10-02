#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.9"
__email__ = "dlemmer@tgen.org"

'''
Created on Mar 4, 2014

@author: dlemmer
'''

"""filter matrix based on a minimum distance between SNPs"""

from Bio import SeqIO

def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Return a reduced matrix of SNPs that are distant from each other.")
    parser.add_argument("-m", "--matrix", required=True, help="Path to the starting SNP matrix.")
    parser.add_argument("-r", "--reference", required=True, help="Path to the reference fasta file.")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output files (matrix and fasta).")
    parser.add_argument("-d", "--distance", default=5000, type=int, help="Minimum distance two SNPs can be from each other. Default: 5000")
    return parser.parse_args()

def _get_contigs(reference):
    contigs = {}
    for seq_record in SeqIO.parse(reference, "fasta"):
        contigs[seq_record.id] = len(seq_record)
    return contigs

def _get_index_range(in_matrix):
    matrix=open(in_matrix, "rU")
    firstLine = matrix.readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    matrix.close()
    return last

def _filter_matrix(matrix, prefix, distance, contigs):
    """filter NASP matrix for SNPS that are far apart"""
    import math
    out_file = "%s_distance_filtered.matrix" % prefix
    in_matrix=open(matrix, "U")
    out_matrix=open(out_file, "w")
    first_line = in_matrix.readline()
    out_matrix.write(first_line,)
    first_fields = first_line.split("\t")
    contig_index = first_fields.index("Contig")
    position_index = first_fields.index("Position")
    edge_distance = math.ceil(distance/2)
    current_contig = ""
    current_position = 0
    contig_size = 0
    for line in in_matrix:
        fields=line.split("\t")
        contig = fields[contig_index]
        position = int(fields[position_index])
        if( contig != current_contig ):
            current_contig = contig
            current_position = 0
            contig_size = contigs[contig]
        if position < edge_distance or position > contig_size - edge_distance:
            continue
        if current_position == 0:
            out_matrix.write(line,)
            current_position = position
        else:
            if position - current_position >= distance:
                out_matrix.write(line,)
                current_position = position
    in_matrix.close()
    out_matrix.close()
    return out_file

def _matrix_to_fasta(in_file, prefix):
    """converts a NASP matrix to fasta format"""
    reduced = []
    in_matrix = open(in_file, "rU")
    out_file = "%s_distance_filtered.fasta" % prefix
    out_fasta = open(out_file, "w")
    last = _get_index_range(in_file)
    for line in in_matrix:
        fields = line.split("\t")
        reduced.append(fields[2:last])
    test=map(list, zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:])+"\n")
    in_matrix.close()
    out_fasta.close()
    return out_file

def main():
    commandline_args = _parse_args()
    contigs = _get_contigs(commandline_args.reference)
    out_matrix = _filter_matrix(commandline_args.matrix, commandline_args.prefix, commandline_args.distance, contigs)
    out_fasta = _matrix_to_fasta(out_matrix, commandline_args.prefix)
    print("%s has been filtered for SNPs that are at least %s bases apart. Output is in %s and %s\n" % (commandline_args.matrix, commandline_args.distance, out_matrix, out_fasta))

if __name__ == "__main__":
    main()

