#!/usr/bin/perl -w

use strict;
use warnings;
#use Data::Dumper;

if( @ARGV != 2 )
{
  print <<EOF;
This script compares a SNP matrix to its reference fasta file, and tries to recreate the strand information.
Because several SNP matrix generators do a very bad job of naming columns, you need to modify your matrix
so that there is a column called "Reference" that contains the reference calls.  If you have more than one
chromosome in your file, you'll also need a column called "Chromosome".

Usage:
make_strand_column.pl reference.fasta snpmatrix.tsv > fixedmatrix.tsv
EOF
  exit();
}

my $referencefastafile = shift();
my $snpmatrixfile = shift();

# Read in fasta file
if( open( my $fastafilehandle, '<', $referencefastafile ) )
{
  my $referencechromosomes = {};
  my $referencefastadata = {};
  my $currentchromosome = 'X';
  my $linefromfile = '';
  while( $linefromfile = <$fastafilehandle> )
  {
    if( $linefromfile =~ /^>/ )
    {
      if( $linefromfile =~ /^>([A-Za-z0-9._-]+)/ )
      {
        $currentchromosome = $1;
        $referencechromosomes->{$currentchromosome} = 1;
        $referencefastadata->{$currentchromosome} = ();
      }
    } else
    {
      if( $linefromfile =~ /^([A-Za-z]+)$/ ){ push( @{$referencefastadata->{$currentchromosome}}, split( '', $1 ) ); }
    }
  }
  close( $fastafilehandle );

  # Read in SNP matrix
  if( open( my $matrixfilehandle, '<', $snpmatrixfile ) )
  {
    my $matrixcolumnposition = -1;
    my $matrixcolumnchrom = -1;
    my $matrixcolumnref = -1;
    $linefromfile = <$matrixfilehandle>;
    chomp( $linefromfile );
    my @matrixcolumns = split( /\t/, $linefromfile );
    my $i = 0;
    while( $i < scalar( @matrixcolumns ) )
    {
      if( $matrixcolumns[$i] =~ /^position$/i ){ $matrixcolumnposition = $i; }
      if( $matrixcolumns[$i] =~ /^chrom(?:osome)?$/i ){ $matrixcolumnchrom = $i; }
      if( $matrixcolumns[$i] =~ /^ref(?:erence)?$/i ){ $matrixcolumnref = $i; }
      $i++;
    }
    if( $matrixcolumnposition >= 0 )
    {
      if( $matrixcolumnref >= 0 )
      {
        print $linefromfile . "\tStrand\n";
        while( $linefromfile = <$matrixfilehandle> )
        {
          chomp( $linefromfile );
          my @columnsfromfile = split( /\t/, $linefromfile );
          $currentchromosome = ( ( $matrixcolumnchrom == -1 ) ? 'X' : $columnsfromfile[$matrixcolumnchrom] );
          my $strandidentifier = "?";
          if( $columnsfromfile[$matrixcolumnref] =~ /^[A-Za-z]$/ )
          {
            my $matrixrefcall = $columnsfromfile[$matrixcolumnref];
            my $matrixrefcomplement = $matrixrefcall;
            $matrixrefcomplement =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            if( !defined( $referencefastadata->{$currentchromosome} ) ){ $currentchromosome = ${keys($referencefastadata)}[0]; }
            if( defined( $referencefastadata->{$currentchromosome}[$columnsfromfile[$matrixcolumnposition]-1] ) )
            {
              if( uc( $referencefastadata->{$currentchromosome}[$columnsfromfile[$matrixcolumnposition]-1] ) eq uc( $matrixrefcall ) ){ $strandidentifier = "+"; }
              elsif( uc( $referencefastadata->{$currentchromosome}[$columnsfromfile[$matrixcolumnposition]-1] ) eq uc( $matrixrefcomplement ) ){ $strandidentifier = "-"; }
            }
          }
          print $linefromfile . "\t$strandidentifier\n";
        }
  
      } else { print STDERR "Could not find 'reference' column in '$snpmatrixfile'!\n"; }
    } else { print STDERR "Could not find 'position' column in '$snpmatrixfile'!\n"; }
    close( $matrixfilehandle );
  } else { print STDERR "Could not open '$snpmatrixfile'!\n"; }
} else { print STDERR "Could not open '$referencefastafile'!\n"; }



