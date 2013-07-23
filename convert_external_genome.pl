#!/usr/bin/perl -w

use strict;
use warnings;
#use Data::Dumper;

my $nucmerpath = "nucmer";
my $deltafilterpath = "delta-filter";

if( @ARGV != 3 )
{
  print <<EOF;
Meant to be called from the pipeline automatically.

Usage:
convert_external_genome.pl <reference.fasta> <external_genome.fasta> <output.frankenfasta>
EOF
  exit();
}

my $referencefilename = shift();
my $externalfilename = shift();
my $outputfilename = shift();

$referencefilename =~ s/'/\\'/g;
if( open( my $externalfilehandle, "<", $externalfilename ) )
{
  my $externalnickname = "external_genome_" . sprintf( '%06u', int( rand() * 1000000 ) );
  if( $externalfilename =~ /^(?:.*\/)?([^\/]+?)(?:\.[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?)?$/ ){ $externalnickname = $1; }
  $externalnickname =~ s/'//g;
  if( open( my $concatenatedhandle, ">", "$externalnickname.concatenated.fasta" ) )
  {
    # Concatenate all chromosomes in external fasta
    print $concatenatedhandle ">$externalnickname\n";
    my $externalgenome = "";
    my $linefromfile = "";
    while( $linefromfile = <$externalfilehandle> )
    {
      if( $linefromfile =~ /^([A-Za-z.-]+)\s*$/ )
      {
        my $fastachunk = $1;
        print $concatenatedhandle "$fastachunk\n";
        $externalgenome .= $fastachunk;
      }
    }
    close( $concatenatedhandle );

    # Make delta file
    print `$nucmerpath '--prefix=$externalnickname' '$referencefilename' '$externalnickname.concatenated.fasta' 2>&1`;
    print `$deltafilterpath -r -o 100 '$externalnickname.delta' 2>&1 > '$externalnickname.filtered.delta'`;
    #print `$deltafilterpath -q -r -o 100 '$externalnickname.delta' 2>&1 > '$externalnickname.filtered.delta'`;

    # Parse delta file into call string
    my $currentchrom = "";
    my $externalcalls = {};
    my $referencecounter;
    my $externalcounter;
    my $finalposition;
    my $externalreverse;
    my $chromsizes = {};
    open( my $deltahandle, "<", "$externalnickname.filtered.delta" );
    while( $linefromfile = <$deltahandle> )
    {
      if( $linefromfile =~ /^>([^ ]+) [^ ]+ (\d+) \d+\s*$/ )
      {
        $currentchrom = $1;
        if( !defined( $chromsizes->{$currentchrom} ) || ( $chromsizes->{$currentchrom} < $2 ) ){ $chromsizes->{$currentchrom} = $2; }
        if( !defined( $externalcalls->{$currentchrom} ) ){ $externalcalls->{$currentchrom} = ""; }
      }
      elsif( $linefromfile =~ /^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$/ )
      {
        $referencecounter = $1;
        $finalposition = $2;
        $externalcounter = $3;
        $externalreverse = ( ( $3 <= $4 ) ? 0 : 1 );
        # Pad any skipped sections
        if( length( $externalcalls->{$currentchrom} ) < ( $referencecounter - 1 ) )
        {
          $externalcalls->{$currentchrom} .= ( 'N' x ( $referencecounter - length( $externalcalls->{$currentchrom} ) - 1 ) );
        }
      }
      elsif( $linefromfile =~ /^(\-?)(\d+)\s*$/ )
      {
        my $distancechange = $2;
        my $externalinsert = ( ( $1 eq '-' ) ? 0 : 1 );
        if( $distancechange == 0 )
        {
          # Last segment
          $distancechange = ( $finalposition - $referencecounter + 1 );
          if( $externalreverse )
          {
            my $alignedsegment = substr( $externalgenome, ( $externalcounter - $distancechange ), $distancechange );
            $alignedsegment = scalar( reverse( $alignedsegment ) );
            $alignedsegment =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            substr( $externalcalls->{$currentchrom}, ( $referencecounter - 1 ), $distancechange, $alignedsegment );
          }
          else { substr( $externalcalls->{$currentchrom}, ( $referencecounter - 1 ), $distancechange, substr( $externalgenome, ( $externalcounter - 1 ), $distancechange ) ); }
        } else
        {
          $distancechange -= 1;
          if( $externalreverse )
          {
            my $alignedsegment = substr( $externalgenome, ( $externalcounter - $distancechange ), $distancechange );
            $alignedsegment = scalar( reverse( $alignedsegment ) );
            $alignedsegment =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
            substr( $externalcalls->{$currentchrom}, ( $referencecounter - 1 ), $distancechange, $alignedsegment );
          }
          else { substr( $externalcalls->{$currentchrom}, ( $referencecounter - 1 ), $distancechange, substr( $externalgenome, ( $externalcounter - 1 ), $distancechange ) ); }
          $referencecounter += $distancechange;
          $externalcounter += ( $externalreverse ? ( 0 - $distancechange ) : $distancechange );
          if( $externalinsert )
          {
            substr( $externalcalls->{$currentchrom}, ( $referencecounter - 1 ), 1, 'N' );
            $referencecounter += 1;
          } else { $externalcounter += ( $externalreverse ? -1 : 1 ); }
        }
      }
    }
    close( $deltahandle );
    # Pad any unaligned ends
    foreach my $currentchrom ( sort( keys( %{$externalcalls} ) ) )
    {
      if( length( $externalcalls->{$currentchrom} ) < $chromsizes->{$currentchrom} )
      {
        $externalcalls->{$currentchrom} .= ( 'N' x ( $chromsizes->{$currentchrom} - length( $externalcalls->{$currentchrom} ) ) );
      }
    }

    # Write data to final output
    if( open( my $outputfilehandle, ">", $outputfilename ) )
    {
      foreach my $currentchrom ( sort( keys( %{$externalcalls} ) ) )
      {
        print $outputfilehandle ">franken::$currentchrom\n";
        print $outputfilehandle "$externalcalls->{$currentchrom}\n";
      }
    } else { print STDERR "Could not open '$outputfilename'!\n"; }
 
  } else { print STDERR "Could not open temp file!\n"; }
  close( $externalfilehandle );
} else { print STDERR "Could not open '$externalfilename'!\n"; }


