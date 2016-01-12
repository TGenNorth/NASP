#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

"""runs complete nasp to plink pipeline"""
import optparse
import os
import sys
import subprocess
import random

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def get_field_index(matrix_in):
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    return last

def transpose_matrix(matrix, last):
    out_matrix = open("tmp.matrix", "w")
    reduced = [ ]
    my_matrix = open(matrix, "rU")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split("\t")
    reduced.append(first_fields)
    for line in my_matrix:
        values = []
        line=line.strip("\n")
        fields = line.split("\t")
        for field in fields[:last]:
            values.append(field.upper())
        reduced.append(values)
    test=map(list, zip(*reduced))
    for x in test:
        out_matrix.write("\t".join(x)+"\n")
    out_matrix.close()

def create_ped_file(matrix, groups_file, prefix, last):
    outfile = open("%s.ped" % prefix, "w")
    groups_dict = {}
    for line in open(groups_file, "U"):
        group_fields = line.split()
        str1=group_fields[0]
        str2=group_fields[1]
        if str1 == "0":
            pass
        else:
            groups_dict.update({str1:str2})
    ref_fields = []
    for line in open(matrix, "rU"):
        fields = line.split()
        if fields[0] == "LocusID":
            pass
        elif fields[0] == "Reference":
            split_fields = line.split()
            for ref_field in split_fields:
                ref_fields.append(ref_field)
        else:
            sample_list = []
            fixed_line = line.strip()
            sample_fields = fixed_line.split()
            if sample_fields[0] in groups_dict:
                # added the arbitrary:   1   0   0   1
                sample_list.append(sample_fields[0]+'\t'+'1'+'\t'+'0'+'\t'+'0'+'\t'+'1')
                num_samples = len(fields)
                try:
                    my_group_info = groups_dict.get(sample_fields[0])
                except:
                    print("names don't match between matrix and groups file, exiting!")
                    sys.exit()
                sample_list.append(my_group_info)
                for i in range(1,num_samples):
                    if sample_fields[i] == ref_fields[i]:
                        sample_list.append("A"+"\t"+"A")
                    elif sample_fields[i] == 'X' or sample_fields[i] == 'N' or sample_fields[i] == '-':
                        sample_list.append("0"+"\t"+"0")
                    else:
                        sample_list.append("C"+"\t"+"C")
                outfile.write("\t".join(sample_list)+"\n")
     outfile.close()

def create_map_file(matrix, prefix):
    outfile = open('%s.map' % prefix,'w')
    chromosome = {}
    count = 0
    for line in open(matrix, 'U'):
        if line.startswith('LocusID') or line.startswith('Reference'):
            pass
        else:
            map_line = []
            fields = line.split()
            coord_split = fields[0].split("::")
            my_coord = coord_split[1]
            if coord_split[0] in chromosome:
                chromo_number = chromosome.get(coord_split[0])
            else:
                count = count + 1
                chromosome[coord_split[0]] = count
                chromo_number = chromosome.get(coord_split[0])
            map_line.append(str(chromo_number)+'\t'+fields[0]+'\t'+'0'+'\t'+my_coord)
            outfile.write("\n".join(map_line)+"\n")
    outfile.close()

def make_groups_dict(groups_file):
    groups_dict = {}
    for line in open(groups_file, "U"):
        fields = line.split()
        try:
            groups_dict[fields[1]].append(fields[0])
        except KeyError:
            groups_dict[fields[1]] = [fields[0]]
    return groups_dict

def convert_nasp_matrix(matrix):
    in_matrix = open(matrix, "U")
    out_matrix = open("redux.matrix", "w")
    firstLine = in_matrix.readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    out_matrix.write("\t".join(first_fields[0:last])+"\n")
    for line in in_matrix:
        new_fields = []
        fields = line.split()
        ref = fields[1]
        new_fields.append(fields[0])
        for field in fields[1:last]:
            if field==ref:
                new_fields.append("0")
            else:
                new_fields.append("1")
        out_matrix.write("\t".join(new_fields)+"\n")
    in_matrix.close()
    out_matrix.close()

def assign_genomes_to_groups(genomes, num_affected):
    sim_affected=random.sample(set(genomes), int(num_affected))
    sim_unaffected = []
    for genome in genomes:
        if genome in sim_affected:
            pass
        else:
            sim_unaffected.append(genome)
    return sim_affected, sim_unaffected

def make_temp_groups_file(affected, unaffected):
    outfile = open("tmp.map", "w")
    for genome in affected:
        outfile.write(str(genome)+"\t"+"1"+"\n")
    for genome in unaffected:
        outfile.write(str(genome)+"\t"+"2"+"\n")
    outfile.close()

def parse_plink_file(plink_out, alpha):
    sig_hits = []
    for line in open(plink_out, "rU"):
        if line.startswith("CHR"):
            pass
        else:
            fields = line.split()
            if fields[0] == "CHR":
                pass
            else:
                if float(fields[8])<float(alpha):
                    sig_hits.append("1")
    return len(sig_hits)

def main(nasp_matrix, groups_file, prefix, alpha):
    last = get_field_index(nasp_matrix)
    transpose_matrix(nasp_matrix, last)
    create_ped_file("tmp.matrix", groups_file, prefix, last)
    create_map_file(nasp_matrix, prefix)
    ab = subprocess.call(['which', 'plink'])
    if ab == 0:
        print("Running Plink")
        subprocess.check_call("plink --noweb --adjust --file %s --assoc --out out > /dev/null 2>&1" % prefix, shell=True)
        print("Finished")
        print("")
        os.system("mv out.assoc.adjusted reference.assoc")
    else:
        print("plink isn't in your PATH and won't be run")
        sys.exit()
    """this part of the script looks for the effect of random associations
    first, I need to parse the reference plink file"""
    sig_hits = parse_plink_file("reference.assoc", alpha)
    if sig_hits == 0:
        print("No sig hits found, exiting..")
        os.system("rm tmp.matrix")
        sys.exit()
    else:
        affected = []
        unaffected = []
        for line in open(groups_file, "U"):
            fields = line.split()
            if fields[1] == "0":
                pass
            elif fields[1] == "1":
                affected.append(fields[0])
            elif fields[1] == "2":
                unaffected.append(fields[0])
            else:
                pass
        all_genomes = affected+unaffected
        total_random_hits = []
        better_hits_than_random = []
        for i in range(1,100):
            """for each iteration, this randomly picks genomes, the sample
            number as those in your affected set"""
            random_affect=random.sample(set(all_genomes), int(len(affected)))
            """how many of these genomes are already in the affected group???"""
            if len(set(random_affect).intersection(affected)) == len(affected):
                pass
            else:
                sim_affected = list(random_affect)
                sim_unaffected = []
                outfile = open("random_genome_ids.txt", "a")
                """make sure we have a unique ped file for each iteration"""
                os.system("rm %s.ped" % prefix)
                for genome in all_genomes:
                    if genome not in sim_affected:
                        sim_unaffected.append(genome)
                outfile.write("\t".join(sim_affected)+"\n")
                make_temp_groups_file(sim_affected, sim_unaffected)
                create_ped_file("tmp.matrix", "tmp.map", prefix, last)
                subprocess.check_call("plink --noweb --adjust --file %s --assoc --out out > /dev/null 2>&1" % prefix, shell=True)
                positive_hits = parse_plink_file("out.assoc.adjusted", alpha)
                total_random_hits.append(positive_hits)
                if positive_hits>=sig_hits:
                    better_hits_than_random.append("1")
        print("True number of associated SNPs at alpha(%s): %s" % (str(alpha),str(sig_hits)))
        print("Average number of random hits: %s" % str(sum(total_random_hits)/int(len(total_random_hits))))
        print("#random iterations better than or equal to true: %s" % len(better_hits_than_random))
        print("p-value: %s" % str(len(better_hits_than_random)/len(total_random_hits)))
    os.system("mv random_genome_ids.txt randomly_selected_genome_ids.txt")
    os.system("rm tmp.matrix")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-m", "--nasp_matrix", dest="nasp_matrix",
                      help="/path/to/nasp matrix [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-g", "--groups_file", dest="groups_file",
                      help="/path/to/groups file [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-p", "--prefix", dest="prefix",
                      help="prefix for output files [REQUIRED]",
                      type="string", action="store")
    parser.add_option("-a", "--alpha", dest="alpha",
                      help="significance alpha, defaults to 0.05",
                      type="float", default = 0.05, action="store")
    options, args = parser.parse_args()

    mandatories = ["nasp_matrix", "groups_file", "prefix"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.nasp_matrix,options.groups_file,options.prefix,options.alpha)
