#!/usr/bin/env perl

use strict;
use warnings;
use lib '/media/lumberyard/bin';
use Cwd;
use nasp;

if( @ARGV < 2 || @ARGV > 4 )
{
  print <<EOF;
This is the distance pipeline, for determining phylogenetic distance of unknown samples
against known samples with a well-determined phylogenetic relationship, using a matrix
of known SNP positions.

The distance pipeline requires a SNP matrix for the known samples, and a folder of 
[currently, already-aligned bam/sam files] for the unknowns.  It will produce a distance
matrix that relates SNP counts between every unknown sample and all others.

Usage:
distance_pipeline.pl <mastermatrix.tsv> <reference.fasta> [bamfolder [outputfolder]]
EOF
  exit();
}

my $currentdirectory = getcwd();
$currentdirectory =~ s/\/+$//;
my $mastermatrixfile = shift();
if( $mastermatrixfile =~ /^[^\/~]/ ){ $mastermatrixfile = "$currentdirectory/$mastermatrixfile"; }
my $referencefastafile = shift();
if( $referencefastafile =~ /^[^\/~]/ ){ $referencefastafile = "$currentdirectory/$referencefastafile"; }
my $bamfilefolder = $currentdirectory;
if( @ARGV > 0 ){ $bamfilefolder = shift(); }
$bamfilefolder =~ s/\/+$//;
if( $bamfilefolder =~ /^[^\/~]/ ){ $bamfilefolder = "$currentdirectory/$bamfilefolder"; }
my $outputfilefolder = "$bamfilefolder/distance_pipeline_results";
if( @ARGV > 0 ){ $outputfilefolder = shift(); }
$outputfilefolder =~ s/\/+$//;
if( $outputfilefolder =~ /^[^\/~]/ ){ $outputfilefolder = "$currentdirectory/$outputfilefolder"; }

nasp( $referencefastafile, $bamfilefolder, $outputfilefolder, $mastermatrixfile );

