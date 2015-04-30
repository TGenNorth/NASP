__author__ = 'jtravis'

import os

def get_path(filename):
    return os.path.join('./test_data', filename)

REFERENCE_FASTA = get_path('reference.fasta')
REFERENCE_DUPS = get_path('duplicates.txt')

GENERAL_STATS = get_path('general_stats.tsv')

GATK_VCF = get_path('gatk.vcf')

PARSE_FASTA = get_path('example.fasta')
