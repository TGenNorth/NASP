#!/usr/bin/env python

"""find autapomorphic SNPs for a single isolate
from a SNP matrix created by NASP"""

from __future__ import print_function
from optparse import OptionParser

def get_field_index(matrix_in):
    """find where the sequence data stops"""
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split()
    last = first_fields.index("#SNPcall")
    return last

def get_genome_index(matrix_in, genome_name):
    """find the index of your genome of interest"""
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split()
    my_genome = first_fields.index(genome_name)
    return my_genome

def get_singleton_snps(matrix_in, my_genome, last):
    """find those lines in your matrix, where a unique
    SNP state is observed in your genome of interest"""
    firstLine = open(matrix_in).readline()
    output = open("unique_snps.txt", "w")
    for line in open(matrix_in):
        if line.startswith("LocusID"):
            pass
        else:
            hits = []
            fields = line.split()
            for x in fields[1:last]:
                if fields[my_genome] != "-" or fields[my_genome] != "X" or fields[my_genome] != "N":
                    if x == fields[my_genome]:
                        hits.append("1")
            if int(len(hits)) == 1:
                output.write(fields[0]+"\n")
            else:
                continue
    output.close()


def test_file(option, opt_str, value, parser):
    try:
        with open(value):
            setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        import sys
        sys.exit()

def main(matrix_in, genome_name):
    last = get_field_index(matrix_in)
    my_genome = get_genome_index(matrix_in, genome_name)
    get_singleton_snps(matrix_in, my_genome, last)


if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--input_matrix", dest="matrix_in",
                      help="/path/to/NASP_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-g", "--genome_name", dest="genome_name",
                      help="genome to look for singleton SNPs [REQUIRED]",
                      action="store", type="string")

    options, args = parser.parse_args()

    mandatories = ["matrix_in", "genome_name"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" % m)
            parser.print_help()
            exit(-1)

    main(options.matrix_in, options.genome_name)
