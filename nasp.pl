#!/usr/bin/perl -w

use strict;
use warnings;
use lib '/media/lumberyard/bin';
use Cwd;
use nasp;

if( @ARGV < 1 || @ARGV > 3 )
{
  print <<EOF;
This is the experimental TGen North GATK/SolSNP pipeline.

Usage:
nasp <reference.fasta> [read_folder [output_folder]]
EOF
  exit();
}

my $currentdirectory = getcwd();
$currentdirectory =~ s/\/+$//;
my $referencefastafile = shift();
if( $referencefastafile =~ /^[^\/~]/ ){ $referencefastafile = "$currentdirectory/$referencefastafile"; }
my $bamfilefolder = $currentdirectory;
if( @ARGV > 0 ){ $bamfilefolder = shift(); }
$bamfilefolder =~ s/\/+$//;
if( $bamfilefolder =~ /^[^\/~]/ ){ $bamfilefolder = "$currentdirectory/$bamfilefolder"; }
my $outputfilefolder = "$bamfilefolder/nasp_results";
if( @ARGV > 0 ){ $outputfilefolder = shift(); }
$outputfilefolder =~ s/\/+$//;
if( $outputfilefolder =~ /^[^\/~]/ ){ $outputfilefolder = "$currentdirectory/$outputfilefolder"; }

nasp( $referencefastafile, $bamfilefolder, $outputfilefolder, '' );

