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

def find_reference_genome(reference, snp_eff):
    try:
        os.system("java -jar %s databases | grep %s > tmp.snpeff.out" % (snp_eff,reference))
    except:
        print("problem finding your genome, refine search string")
        sys.exit()
    num_lines = sum(1 for line in open('tmp.snpeff.out'))
    if num_lines == 0:
        print("no hits found for your search criteria, try again with the strain name")
        sys.exit()
    elif num_lines > 1:
        print("multiple hits found for your search criteria, refine search criteria")
        os.system("cat tmp.snpeff.out")
        sys.exit()
    else:
        for line in open("tmp.snpeff.out"):
            fields = line.split()
            reference_id = fields[0]
    os.system("rm tmp.snpeff.out")
    return reference_id

def make_map_file(ref_id, snp_eff, vcf):
    os.system("java -jar %s eff -v -download -no-downstream -no-upstream %s %s > tmp.snpeff.out 2>&1" % (snp_eff,ref_id,vcf))
    with open("tmp.snpeff.out") as myFile:
        for num, line in enumerate(myFile, 1):
            if "# Number of chromosomes" in line:
                fields = line.split()
                num_chromes = int(fields[-1])
                ref_line = num
    line_end = num_chromes+ref_line
    print("")
    print("create mapping file with the following information:")
    with open("tmp.snpeff.out") as myFile:
        for line in myFile.readlines()[ref_line:line_end+1]:
            print(line)
    chromosomes = []
    for line in open(vcf):
        fields = line.split()
        if line.startswith("#"):
            pass
        else:
            if fields[0].split("::") in chromosomes:
                pass
            else:
                chromosomes.append(fields[0].split("::"))
    nr=[x for i, x in enumerate(chromosomes) if x not in chromosomes[i+1:]]
    print("")
    print("Your chromosome names are:")
    for chromosome in nr:
        print("".join(chromosome))
    os.system("rm tmp.snpeff.out")
    sys.exit()

def main(in_matrix, vcf, snp_eff, map, reference):
    ref_id = find_reference_genome(reference, snp_eff)
    if "NULL" in map:
        make_map_file(ref_id, snp_eff, vcf)
    else:
        my_map = make_map_dict(map)
    rename_vcf_file(vcf, my_map)
    run_snpeff(snp_eff, ref_id)
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
                      help="/path/to/mapping file",
                      action="callback", callback=test_file, default="NULL", type="string")
    options, args = parser.parse_args()
    mandatories = ["in_matrix", "vcf", "snp_eff", "reference"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.vcf,options.snp_eff,options.map,options.reference)
