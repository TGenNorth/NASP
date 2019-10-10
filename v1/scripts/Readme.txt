###NASP post-matrix scripts###

-All scripts are written for use with Python3

1. report_single_snps_single_isolate.py

-What does it do?

If given a genome name, which must match the name in the NASP matrix, this script will
generate a list of SNPs that are unique to that genome

Example:

python ../scripts/report_single_snps_single_isolate.py -m bestsnp.tsv -g Reference

-Output file: "unique_snps.txt"

-These SNPs can then be further investigated or removed using scripts described below

2. matrix_to_fasta.py

-What does it do?

The script converts a NASP matrix into a series of FASTA files. The input is a NASP matrix,
and an optional parameter controls how much missing data makes it into your FASTA file
as missing data "-". The script also creates a FASTA file that only contains parsimony-
informative characters, which may be useful for certain applications.

Example:

python ../scripts/matrix_to_fasta.py -m missingdata.tsv -p 90 -f 0.90

Total SNPs: 4467
number of SNPs after filtering: 1749
number of parsimony-informative SNPs: 1209

-Output files:

-prefix.raw.fasta (Direct conversion of NASP to FASTA)
-prefix.filtered.fasta (FASTA where SNPs greater than a threshold are kept. Any ambiguous, or
 missing data is converted to a dash "-"
-prefix.filtered_PI_snps_only.fasta (FASTA where SNPs greater than a threshold are kept). Any
 ambiguous or missing data is converted to a dash "-". Only PI SNPs, or those that have two
 or more states in more than two genomes, are retained
-clean_matrix.txt (Converted NASP matrix that corresponds to prefix.filtered.fasta)
-clean_PI_matrix.txt (Converted NASP matrix that corresponds to prefix.filtered_PI_snps_only.fasta)

3. filter_matrix_by_genome.py

-What does it do?

Given a NASP matrix, this script either removes genomes or keeps genomes from a NASP matrix.
Fields following the DNA calls aren't altered with this script.

Example:

python ../scripts/filter_matrix_by_genome.py -m bestsnp.tsv -p filtered -g to_prune.txt -a keep

-Output file:

-Filtered_genomes.matrix: New NASP matrix with genomes removed or kept

4. filter_matrix_by_coord.py

-What does it do?

Given a NASP matrix and a set of coordinates, this script will produce a new matrix with only
those coordinates either kept or removed

Example:

python ../scripts/filter_matrix_by_coord.py -i bestsnp.tsv -p filtered -c snp_coords.txt -a remove
Number of SNPs matching in your list: 1
number of SNPs after missing regions removed: 1748
number of parsimony-informative SNPs: 1488

5. filter_matrix_by_distance.py

-What does it do?

Given a NASP SNP matrix, corresponding reference fasta, and a minimum distance, this script will produce
a new matrix with only SNPs that are guaranteed to be at least that distance apart from each other.
This can be used to create a random sampling of SNPs that are spread throughout the genome.

Example:

python ../scripts/filter_matrix_by_distance.py -m bestsnp.tsv -r reference.fasta -p prefix -d 5000
bestsnp.tsv has been filtered for SNPs that are at least 5000 bases apart. Output is in
prefix_distance_filtered.matrix and prefix_distance_filtered.fasta

6. annotate_NASP.py

-What does it do?

Provides functional information on a set of SNPs using snpEff (tested version is 4.2). For the script to work, you will
need the name of the snpEff database and a mapping file that correlates the names of your input
FASTA with the names in the snpEff database

Example

-First, I run the script without a mapping file, using Y. pestis CO92 as the reference:

[js2829@wind /scratch/js2829/nasp_bowtie/snpeff_test ]$ python annotate_NASP.py -i bestsnp.tsv -v bestsnp.vcf -s /common/contrib/tools/snpEff/snpEff.jar -r CO92

create mapping file with the following information:
# Chromosomes                : Format 'chromo_name size codon_table'

#		'NC_003143'	4653728	Standard
#		'NC_003134'	96210	Standard
#		'NC_003131'	70305	Standard
#		'NC_003132'	9612	Standard

Your chromosome names are:
gi|16120353|ref|NC_003143.1|
gi|5763810|emb|AL109969.1|

-Now you must create a mapping file that contains your names and snpEff's names:

gi|16120353|ref|NC_003143.1|	NC_003143
gi|5834685|emb|AL117211.1|	NC_003134
gi|5832423|emb|AL117189.1|	NC_003131
gi|5763810|emb|AL109969.1|	NC_003132

-Now you can annotate your matrix:

python annotate_NASP.py -i bestsnp.tsv -v bestsnp.vcf -s /common/contrib/tools/snpEff/snpEff.jar -r CO92 -m map.txt

-The output is a new NASP matrix ("annotated_bestsnp.tsv") with three extra columns: type    locus   ncbi_id

7. nasp_to_plink.py

-What does it do?

If provided a groups file, this script creates the input files for Plink (tested version is 1.07), runs Plink, then randomly shuffles the genomes into groups and
calculates how many signficant SNPs are identified between phenotypes

Example:

The groups files looks like:

Reference 0
Genome1 1
Genome2 1
Genome3 2
Genome4 2

-"1" is unaffected, "2" is affected, and "0" is unknown

python ../scripts/nasp_to_plink.py -m bestsnp.tsv -g groups.txt -p test
/common/contrib/bin/plink
Running Plink
Finished
True number of associated SNPs at alpha(0.05): 270
Average number of random hits: 43.2222222222
#random iterations better than or equal to true: 4
p-value: 0.040404040404

-Output files:

1. reference.assoc : These are SNPs that are positively associated at the given alpha
2. random_genome_ids.txt : The IDs of genomes selected at random in p-value iterations
