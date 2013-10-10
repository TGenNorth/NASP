#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "0.6"
__email__ = "dsmith@tgen.org"

import re


def _parse_args():
    'Parses command line arguments.'
    import argparse
    parser = argparse.ArgumentParser( description="Meant to be called from the pipeline automatically." )
    parser.add_argument( "--input-files", nargs="+", required=True, help="Input VCF files for matrix conversion." )
    parser.add_argument( "--output-matrix", required=True, help="Name of output matrix to create." )
    parser.add_argument( "--minimum-coverage", type=int, default=10, help="Minimum coverage depth at a position; call will be changed to 'X' if coverage is below this value." )
    parser.add_argument( "--minimum-proportion", type=float, default=0.9, help="Minimum proportion of reads that much match the call at a position; call will be changed to 'N' if proportion is below this value." )
    parser.add_argument( "--filter-logfile", help="Name of filter logfile to create." )
    parser.add_argument( "--allcall-matrix", help="Name of output matrix with all callable positions included." )
    parser.add_argument( "--nonx-matrix", help="Name of output matrix with positions that contain any N, X, or indel filtered." )
    parser.add_argument( "--output-fasta", help="Name of output snpfasta to create." )
    parser.add_argument( "--nonx-fasta", help="Name of output snpfasta with positions that contain any N, X, or indel filtered." )
    return( parser.parse_args() )

def read_vcf_file( vcf_filename, snp_matrix ):
    'Reads in a VCF file, and stores the data in the passed SNP matrix object.'
    with open( vcf_filename, 'r' ) as vcf_filehandle:
        import vcf
        vcf_nickname = re.search( r'([^/]+?)(?:\.[Vv][Cc][Ff])?$', vcf_filename )
        if vcf_nickname: vcf_nickname = vcf_nickname.group(1)
        else: vcf_nickname = vcf_filename
        vcf_data_handle = vcf.Reader( vcf_filehandle )
        vcf_samples = vcf_data_handle.samples
        print( vcf_samples )
        for vcf_record in vcf_data_handle:
            

def main():
    commandline_args = _parse_args()
    import snp_matrix_object
    snp_matrix = snp_matrix_object.snp_matrix()
    for vcf_filename in commandline_args.input_files:
        read_vcf_file( vcf_filename, snp_matrix )

if __name__ == "__main__": main()


