[![Actions Status](https://github.com/TGenNorth/NASP/workflows/.github/workflows/main.yml/badge.svg)](https://github.com/TGenNorth/NASP/actions)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/nasp/badges/platforms.svg)](https://anaconda.org/bioconda/nasp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/nasp/badges/license.svg)](https://anaconda.org/bioconda/nasp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/nasp/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

# DEVELOPMENT BRANCH

# Overview

NASP is a pipeline for analysis of genomic data. It is a suite of tools meant to collect and report on statistically-relevant high-confidence positions in a collection of genomes, with emphasis on variant positions, especially single nucleotide polymorphisms (SNPs). NASP expects some combination of files in FASTA, FASTQ, SAM, BAM, and VCF format as input, and will produce output files also in those formats. As NASP is a pipeline, it expects to link a set of external tools (usually installed separately) to complete specific analysis tasks.

# Installation

`conda install -c bioconda 'nasp>1.1.2'`

# Snakemake Workflow

The NASP dispatcher is in the process of being replaced with a Snakemake workflow. The following is a proposed usage:

```bash
git clone --branch workflows https://github.com/TGenNorth/NASP.git workflows
snakemake -j -s workflows/Snakemake --config reference=reference.fasta assembly_dir=assemblies pe_reads=pe_reads matrix
```

https://katacoda.com/corburn/scenarios/nasp

# Publication
Please read our paper for more information:

Jason W. Sahl, Darrin Lemmer, Jason Travis, James M. Schupp, John D. Gillece, Maliha Aziz, Elizabeth M. Driebe, Kevin P. Drees, Nathan D. Hicks, Charles Hall Davis Williamson, Crystal M. Hepp, David Earl Smith, Chandler Roe, David M. Engelthaler, David M. Wagner, Paul Keim

"NASP: an accurate, rapid method for the identification of SNPs in WGS datasets that supports flexible input and output formats". Published Ahead of Print: 21 June, 2016 Microbial Genomics doi: 10.1099/mgen.0.000074 (http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000074)

CONTACT:
TGen North
3051 W Shamrell Blvd Ste 106
Flagstaff, AZ 86001-9435
Darrin Lemmer
dlemmer@tgen.org
+1-928-226-6374
