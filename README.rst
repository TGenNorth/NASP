.. image:: https://landscape.io/github/TGenNorth/NASP/tests/landscape.svg?style=flat
   :target: https://landscape.io/github/TGenNorth/NASP/tests
   :alt: Code Health

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

The Northern Arizona SNP Pipeline (NASP)
========================================

OVERVIEW:
---------

NASP is a pipeline for analysis of genomic data. It is a suite of tools
meant to collect and report on statistically-relevant high-confidence
positions in a collection of genomes, with emphasis on variant
positions, especially single nucleotide polymorphisms (SNPs). NASP
expects some combination of files in FASTA, FASTQ, SAM, BAM, and VCF
format as input, and will produce output files also in
those formats. As NASP is a pipeline, it expects to link a set of
external tools (usually installed separately) to complete specific
analysis tasks.

USAGE:
------

Usage depends upon the installation method used on your system, and the
user interface you select. For standard installations with the
command-line interface, you would collect (or symbolically link) all of
your input files into a folder, and then invoke the command-line
interface from that folder. Expected usage for the command-line
interface is:

`nasp.py [output\_folder]`

You will then be prompted to answer a few questions about your analysis.

Optionally, if you are re-running a previous analysis with the same (or
similar) options, you can pass in an xml-based configuration file (this
is written out to your output\_folder after running the command-line
interface) using the format:

`nasp.py --config`

INSTALLATION:
-------------

See the main page for documentation (http://tgennorth.github.io/NASP/).

DEPENDENCIES:
-------------

For information about external tools that are required, or can be
utilized, and those versions that have been tested to work with NASP,
refer to the main NASP page (http://tgennorth.github.io/NASP/)

LICENSE:
--------

Copyright |copy| The Translational Genomics Research Institute See the
included "LICENSE" document.

CONTACT:
--------

| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

| Darrin Lemmer
| dlemmer@tgen.org
| +1-928-226-6374
