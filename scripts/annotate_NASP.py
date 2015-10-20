#!/usr/bin/env python

"""controller for SNPeff annotation of a NASP matrix"""

from __future__ import print_function
import sys
import os
from optparse import OptionParser
from subprocess import Popen


def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()

def make_map_dict(map):
    my_dict = {}
    for line in open(map, "rU"):
        new_line = line.strip()
        fields = new_line.split()
        my_dict.update({fields[0]:fields[1]})
    return my_dict

def rename_vcf_file(vcf, my_map):
    outfile = open("tmp.vcf_xyx", "w")
    for line in open(vcf, "rU"):
        if line.startswith("#"):
            outfile.write(line)
        else:
            fields = line.split()
            fields[0] = my_map.get(fields[0])
            outfile.write("\t".join(fields)+"\n")
    outfile.close()

def run_snpeff(snp_eff, reference):
    args = ['java', '-jar', '%s' % snp_eff, 'eff', '-download',
    '-no-downstream','-no-upstream','%s' % reference,'tmp.vcf_xyx']
    try:
        vcf_fh = open('tmp.snpeff.out', 'w')
    except:
        print('could not open snpeff file')
    try:
        log_fh = open('tmp.snpeff.log', 'w')
    except:
        print('could not open log file')
    snpeff_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
    snpeff_run.wait()

def parse_output_file():
    outfile = open("tmp.matrix.xyx", "w")
    line_list = []
    tmp_list = []
    tmp_list.append("type")
    tmp_list.append("locus")
    tmp_list.append("ncbi_id")
    line_list.append(tmp_list)
    for line in open("tmp.snpeff.out", "rU"):
        fields_of_interest = []
        new_line = line.strip()
        if new_line.startswith("#"):
            pass
        else:
            fields = new_line.split()
            good_fields = fields[7].split("|")
            if len(good_fields)>6:
                fields_of_interest.append(good_fields[1])
                fields_of_interest.append(good_fields[3])
                fields_of_interest.append(good_fields[4])
            else:
                fields_of_interest.append("N/A")
                fields_of_interest.append("N/A")
                fields_of_interest.append("N/A")
        line_list.append(fields_of_interest)

    for alist in line_list:
        outfile.write("\t".join(alist)+"\n")
    outfile.close()

def main(in_matrix, vcf, snp_eff, map, reference):
    my_map = make_map_dict(map)
    rename_vcf_file(vcf, my_map)
    run_snpeff(snp_eff, reference)
    parse_output_file()
    """removes the white space lines"""
    os.system("awk '$NF>0' tmp.matrix.xyx > pruned.matrix.xyx")
    os.system("paste %s pruned.matrix.xyx > annotated_%s" % (in_matrix, in_matrix))
    os.system("rm tmp.matrix.xyx pruned.matrix.xyx tmp.snpeff.out tmp.vcf_xyx tmp.snpeff.log")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--in_matrix", dest="in_matrix",
                      help="/path/to/NASP_matrix_to_annotate [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-v", "--vcf", dest="vcf",
                      help="/path/to/VCF file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-s", "--snp_eff", dest="snp_eff",
                      help="/path/to/snpEff jar file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--reference", dest="reference",
                      help="name of snpEff reference [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-m", "--map", dest="map",
                      help="/path/to/mapping file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    mandatories = ["in_matrix", "vcf", "snp_eff", "reference", "map"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.vcf,options.snp_eff,options.map,options.reference)

