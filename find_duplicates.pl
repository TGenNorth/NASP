#!/usr/bin/env perl

use strict;
use warnings;
#use Data::Dumper;

my $nucmerpath = "nucmer";

if( @ARGV != 2 )
{
  print <<EOF;
Meant to be called from the pipeline automatically.

Usage:
find_duplicates.pl <reference.fasta> <output.txt>
EOF
  exit();
}

my $referencefilename = shift();
my $outputfilename = shift();

$referencefilename =~ s/'/\\'/g;
# Make delta file
print `$nucmerpath --prefix=reference --maxmatch --nosimplify '$referencefilename' '$referencefilename' 2>&1`;

# Parse delta file into dups string
my $chromosomea = "";
my $chromosomeb = "";
my $isduplicate = {};
if( open( my $deltahandle, "<", "reference.delta" ) )
{
  my $linefromfile = '';
  while( $linefromfile = <$deltahandle> )
  {
    if( $linefromfile =~ /^>([^ ]+) ([^ ]+) (\d+) (\d+)\s*$/ )
    {
      $chromosomea = $1;
      $chromosomeb = $2;
      if( !defined( $isduplicate->{$chromosomea} ) ){ $isduplicate->{$chromosomea} = ""; }
      if( !defined( $isduplicate->{$chromosomeb} ) ){ $isduplicate->{$chromosomeb} = ""; }
      if( length( $isduplicate->{$chromosomea} ) < $3 ){ $isduplicate->{$chromosomea} .= ( '0' x ( $3 - length( $isduplicate->{$chromosomea} ) ) ); }
      if( length( $isduplicate->{$chromosomeb} ) < $4 ){ $isduplicate->{$chromosomeb} .= ( '0' x ( $4 - length( $isduplicate->{$chromosomeb} ) ) ); }
    }
    elsif( $linefromfile =~ /^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$/ )
    {
      my $chromosomeastart = $1;
      my $chromosomeaend = $2;
      my $chromosomebstart = $3;
      my $chromosomebend = $4;
      if( ( $chromosomea ne $chromosomeb ) || ( $chromosomeastart != $chromosomebstart ) )
      {
        if( $chromosomeaend < $chromosomeastart )
        {
          my $temp = $chromosomeastart;
          $chromosomeastart = $chromosomeaend;
          $chromosomeaend = $temp;
        }
        if( $chromosomebend < $chromosomebstart )
        {
          my $temp = $chromosomebstart;
          $chromosomebstart = $chromosomebend;
          $chromosomebend = $temp;
        }
        substr( $isduplicate->{$chromosomea}, ( $chromosomeastart - 1 ), ( $chromosomeaend - $chromosomeastart + 1 ), ( '1' x ( $chromosomeaend - $chromosomeastart + 1 ) ) );
        substr( $isduplicate->{$chromosomeb}, ( $chromosomebstart - 1 ), ( $chromosomebend - $chromosomebstart + 1 ), ( '1' x ( $chromosomebend - $chromosomebstart + 1 ) ) );
      }
    }
  }
  close( $deltahandle );
  
  # Write data to final output
  if( open( my $outputfilehandle, ">", $outputfilename ) )
  {
    foreach my $currentchrom ( sort( keys( %{$isduplicate} ) ) )
    {
      print $outputfilehandle ">$currentchrom\n";
      print $outputfilehandle "$isduplicate->{$currentchrom}\n";
    }
  } else { print STDERR "Could not open '$outputfilename'!\n"; }

} else { print STDERR "Could not open delta file!\n"; }


