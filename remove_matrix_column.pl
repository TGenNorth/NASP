#!/usr/bin/perl -w

use strict;
use warnings;
#use Data::Dumper;

if( @ARGV < 6 )
{
  print <<EOF;
Best results if you provide an allcallable matrix.
This script is not indel-safe, and will completely throw out any sample that contains an indel.

Usage:
remove_matrix_column.pl <original_matrix.tsv> [sample1toremove [sample2toremove ...]] <output_matrix.tsv> <noNXindel_matrix.tsv> <allcallable_matrix.tsv> <output_matrix.snpfasta> <noNXindel_matrix.snpfasta>
EOF
  exit();
}

my $oldmatrixfile = shift();
my @samplestoremove = ();
while( my $sampletoremove = shift() ){ push( @samplestoremove, $sampletoremove ); }
my $nonxsnpfastafile = pop( @samplestoremove );
my $snpfastafile = pop( @samplestoremove );
my $allcallmatrixfile = pop( @samplestoremove );
my $nonxmatrixfile = pop( @samplestoremove );
my $outputmatrixfile = pop( @samplestoremove );

if( open( my $oldmatrixhandle, '<', $oldmatrixfile ) )
{

  # Read in original SNP matrix
  my $positioncolumn = -1;
  my $referencecolumn = -1;
  my $chromosomecolumn = -1;
  my $strandcolumn = -1;
  my $linefromfile = <$oldmatrixhandle>;
  my @columnheadings = split( /\t/, $linefromfile );
  my $snppositions = {};
  my $samplecallsdata = {};
  my $columnswedontwant = {};
  my $referencecalls = {};
  my $i = 0;
  while( $i < scalar( @columnheadings ) )
  {
    if( $columnheadings[$i] =~ /^position$/i ){ $positioncolumn = $i; }
    if( $columnheadings[$i] =~ /^ref(?:erence)?$/i ){ $referencecolumn = $i; }
    if( $columnheadings[$i] =~ /^chrom(?:osome)?$|^contig$/i ){ $chromosomecolumn = $i; }
    if( $columnheadings[$i] =~ /^strand$/i ){ $strandcolumn = $i; }
    $i++;
  }
  if( ( $positioncolumn >= 0 ) && ( $referencecolumn >= 0 ) )
  {
    foreach my $sampletoremove (@samplestoremove){ $columnswedontwant->{$sampletoremove} = 1; }
    $columnswedontwant->{$columnheadings[$positioncolumn]} = 1;
    $columnswedontwant->{$columnheadings[$referencecolumn]} = 1;
    if( $chromosomecolumn >= 0 ){ $columnswedontwant->{$columnheadings[$chromosomecolumn]} = 1; }
    if( $strandcolumn >= 0 ){ $columnswedontwant->{$columnheadings[$strandcolumn]} = 1; }
    while( $linefromfile = <$oldmatrixhandle> )
    {
      chomp( $linefromfile );
      my @columnsfromfile = split( /\t/, $linefromfile );
      my $currentchromosome = ( ( $chromosomecolumn == -1 ) ? 'X' : $columnsfromfile[$chromosomecolumn] );
      if( !defined( $snppositions->{$currentchromosome} ) ){ $snppositions->{$currentchromosome} = {}; }
      $snppositions->{$currentchromosome}{$columnsfromfile[$positioncolumn]} = 1;
      if( !defined( $referencecalls->{$currentchromosome} ) ){ $referencecalls->{$currentchromosome} = {}; }
      $referencecalls->{$currentchromosome}{$columnsfromfile[$positioncolumn]} = $columnsfromfile[$referencecolumn];
      $i = 0;
      while( $i < scalar( @columnheadings ) )
      {
        if( !defined( $columnsfromfile[$i] ) || $columnsfromfile[$i] !~ /^[A-Za-z]$/ )
        {
          $columnswedontwant->{$columnheadings[$i]} = 1;
          if( defined( $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]} ) ){ delete( $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]} ); }
        } elsif( !defined( $columnswedontwant->{$columnheadings[$i]} ) )
        {
          if( !defined( $samplecallsdata->{$oldmatrixfile} ) ){ $samplecallsdata->{$oldmatrixfile} = {}; }
          if( !defined( $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]} ) ){ $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]} = {}; }
          if( !defined( $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]}{$currentchromosome} ) ){ $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]}{$currentchromosome} = {}; }
          $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]}{$currentchromosome}{$columnsfromfile[$positioncolumn]} = $columnsfromfile[$i];
          if( ( $strandcolumn >= 0 ) && ( $columnsfromfile[$strandcolumn] eq "-" ) )
          {
            $samplecallsdata->{$oldmatrixfile}{$columnheadings[$i]}{$currentchromosome}{$columnsfromfile[$positioncolumn]} =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
          }
        }
        $i++;
      }
    }
    #print Dumper( $samplecallsdata );

    # Make new SNP matrix
    my $snpfastadata = {};
    my $nonxsnpfastadata = {};
    my $refsnpfastadata = '';
    my $nonxrefsnpfastadata = '';
    if( open( my $matrixfilehandle, '>', $outputmatrixfile ) && open( my $nonxfilehandle, '>', $nonxmatrixfile ) && open( my $allcallfilehandle, '>', $allcallmatrixfile ) )
    {
      my $allpatternarray = {};
      my $allnextpatternnum = 1;
      my $snppatternarray = {};
      my $snpnextpatternnum = 1;
      my $nonxpatternarray = {};
      my $nonxnextpatternnum = 1;
      print $matrixfilehandle "LocusID\tReference\t";
      print $nonxfilehandle "LocusID\tReference\t";
      print $allcallfilehandle "LocusID\tReference\t";
      foreach my $samplefile ( sort keys( $samplecallsdata ) )
      {
        foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
        {
          print $matrixfilehandle "$samplecolumn\t";
          print $nonxfilehandle "$samplecolumn\t";
          print $allcallfilehandle "$samplecolumn\t";
        }
      }
      print $matrixfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tNotes\t\n";
      print $nonxfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tNotes\t\n";
      print $allcallfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tNotes\t\n";
      foreach my $currentchromosome ( sort keys( %$snppositions ) )
      {
        foreach my $currentposition ( sort {$a <=> $b} keys( $snppositions->{$currentchromosome} ) )
        {
          my $linetoprint = '';
          my $allelecounts = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0, 'ref' => 0, 'snp' => 0, 'sum' => 0 };
          if( !defined( $referencecalls->{$currentchromosome}{$currentposition} ) ){ $referencecalls->{$currentchromosome}{$currentposition} = "N"; }
          my $patternassignments = {};
          my $callpattern = '';
          if( $referencecalls->{$currentchromosome}{$currentposition} =~ /^[Nn]$/ ){ $callpattern .= 'N'; } else
          {
            $patternassignments->{$referencecalls->{$currentchromosome}{$currentposition}} = 1;
            $callpattern .= '1';
          }
          my $nextassignmentnum = 2;
          $linetoprint .= "${currentchromosome}::${currentposition}\t$referencecalls->{$currentchromosome}{$currentposition}\t";
          foreach my $samplefile ( sort keys( $samplecallsdata ) )
          {
            foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
            {
              if( !defined( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} ) ){ $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} = "N"; }
              $linetoprint .= "$samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}\t";
              if( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[ACGTacgt]$/ ){ $allelecounts->{uc($samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition})}++; }
              elsif( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[Uu]$/ ){ $allelecounts->{"T"}++; }
              else { $allelecounts->{"N"}++; }
              if( ( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[ACGTUacgtu]$/ ) && ( $referencecalls->{$currentchromosome}{$currentposition} =~ /^[ACGTUacgtu]$/ ) )
              {
                if( ( uc( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} ) eq uc( $referencecalls->{$currentchromosome}{$currentposition} ) ) || ( ( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[TUtu]$/ ) && ( $referencecalls->{$currentchromosome}{$currentposition} =~ /^[TUtu]$/ ) ) )
                {
                  $allelecounts->{"ref"}++;
                } else { $allelecounts->{"snp"}++; }
              }
              $allelecounts->{"sum"}++;
              if( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[NnXx]$/ )
              {
                $callpattern .= "N";
              } else
              {
                if( !defined( $patternassignments->{$samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}} ) )
                {
                  $patternassignments->{$samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}} = $nextassignmentnum;
                  if( $nextassignmentnum ne '+' ){ $nextassignmentnum = ( ( $nextassignmentnum < 9 ) ? ( $nextassignmentnum + 1 ) : '+' ); }
                }
                $callpattern .= $patternassignments->{$samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}};
              }
            }
          }
          $linetoprint .= "$allelecounts->{'snp'}\t$allelecounts->{'ref'}\t$allelecounts->{'A'}\t$allelecounts->{'C'}\t$allelecounts->{'G'}\t$allelecounts->{'T'}\t$allelecounts->{'N'}\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"snp"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"ref"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"A"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"C"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"G"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"T"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= sprintf( "%.2f", ( $allelecounts->{"N"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          $linetoprint .= "$currentchromosome\t$currentposition\t";
          if( $allelecounts->{'snp'} > 0 )
          {
            if( !defined( $snppatternarray->{$callpattern} ) ){ $snppatternarray->{$callpattern} = $snpnextpatternnum++; }
            print $matrixfilehandle $linetoprint . "'$callpattern'\t$snppatternarray->{$callpattern}\t\t\n";
            $refsnpfastadata .= ( ( length( $referencecalls->{$currentchromosome}{$currentposition} ) == 1 ) ? $referencecalls->{$currentchromosome}{$currentposition} : 'N' );
            foreach my $samplefile ( sort keys( $samplecallsdata ) )
            {
              foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
              {
                if( !defined( $snpfastadata->{$samplefile} ) ){ $snpfastadata->{$samplefile} = {}; }
                if( !defined( $snpfastadata->{$samplefile}{$samplecolumn} ) ){ $snpfastadata->{$samplefile}{$samplecolumn} = ''; }
                if( length( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} ) == 1 )
                {
                  $snpfastadata->{$samplefile}{$samplecolumn} .= $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition};
                } else { $snpfastadata->{$samplefile}{$samplecolumn} .= 'N'; }
              }
            }
            if( $allelecounts->{"N"} == 0 )
            {
              if( !defined( $nonxpatternarray->{$callpattern} ) ){ $nonxpatternarray->{$callpattern} = $nonxnextpatternnum++; }
              print $nonxfilehandle $linetoprint . "'$callpattern'\t$nonxpatternarray->{$callpattern}\t\t\n";
              $nonxrefsnpfastadata .= ( ( length( $referencecalls->{$currentchromosome}{$currentposition} ) == 1 ) ? $referencecalls->{$currentchromosome}{$currentposition} : 'N' );
              foreach my $samplefile ( sort keys( $samplecallsdata ) )
              {
                foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
                {
                  if( !defined( $nonxsnpfastadata->{$samplefile} ) ){ $nonxsnpfastadata->{$samplefile} = {}; }
                  if( !defined( $nonxsnpfastadata->{$samplefile}{$samplecolumn} ) ){ $nonxsnpfastadata->{$samplefile}{$samplecolumn} = ''; }
                  $nonxsnpfastadata->{$samplefile}{$samplecolumn} .= $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition};
                }
              }
            }
          }
          if( !defined( $allpatternarray->{$callpattern} ) ){ $allpatternarray->{$callpattern} = $allnextpatternnum++; }
          print $allcallfilehandle $linetoprint . "'$callpattern'\t$allpatternarray->{$callpattern}\t\t\n";
        }
      }
      close( $matrixfilehandle );
      close( $nonxfilehandle );
      close( $allcallfilehandle );
    } else { print STDERR "Could not open '$outputmatrixfile'!\n"; }
    
    if( open( my $snpfastahandle, '>', $snpfastafile ) && open( my $nonxsnpfastahandle, '>', $nonxsnpfastafile ) )
    {
      if( length( $refsnpfastadata ) )
      {
        print $snpfastahandle ">SNP::Reference\n";
        print $snpfastahandle "$refsnpfastadata\n";
        foreach my $samplefile ( sort keys( $samplecallsdata ) )
        {
          foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
          {
            print $snpfastahandle ">SNP::$samplecolumn\n";
            print $snpfastahandle "$snpfastadata->{$samplefile}{$samplecolumn}\n";
          }
        }
      }
      if( length( $nonxrefsnpfastadata ) )
      {
        print $nonxsnpfastahandle ">SNP::Reference\n";
        print $nonxsnpfastahandle "$nonxrefsnpfastadata\n";
        foreach my $samplefile ( sort keys( $samplecallsdata ) )
        {
          foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
          {
            print $nonxsnpfastahandle ">SNP::$samplecolumn\n";
            print $nonxsnpfastahandle "$nonxsnpfastadata->{$samplefile}{$samplecolumn}\n";
          }
        }
      }
      close( $snpfastahandle );
      close( $nonxsnpfastahandle );
    } else { print STDERR "Could not open '$snpfastafile'!\n"; }

  } else { print STDERR "Could not find 'position' or 'reference' column in '$oldmatrixfile'!\n"; }
  close( $oldmatrixhandle );
} else { print STDERR "Could not open '$oldmatrixfile'!\n"; }

