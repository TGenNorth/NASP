__author__ = 'jtravis'

from pkg_resources import resource_filename

PARSE_FASTA = resource_filename(__name__, 'example.fasta')

REFERENCE_FASTA = resource_filename(__name__, 'reference.fasta')
REFERENCE_DUPS = resource_filename(__name__, 'duplicates.txt')
REFERENCE_DELTA = resource_filename(__name__, 'reference.delta')

GENERAL_STATS = resource_filename(__name__, 'general_stats.tsv')

GATK_VCF = resource_filename(__name__, 'gatk.vcf')
