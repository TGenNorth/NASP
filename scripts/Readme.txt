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

5. annotate_NASP.py

-What does it do?

Provides functional information on a set of SNPs using SNPEff
