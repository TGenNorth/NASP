#!/usr/bin/perl -w
#/usr/bin/perl -w -d:SmallProf

use strict;
use warnings;
use Vcf;
#use Data::Dumper;

if( @ARGV < 10 )
{
  print <<EOF;
Meant to be called from the pipeline automatically.

Usage:
vcf_to_matrix.pl <min_coverage> <min_proportion> <reference.fasta> <samples.vcf> [samples2.vcf ...] <variant_matrix.tsv> <bestsnp_matrix.tsv> <allcallable_matrix.tsv> <variant_matrix.snpfasta> <bestsnp_matrix.snpfasta> <statistics_file.tsv>
EOF
  exit();
}

my $mincoverage = shift();
my $minproportion = shift();
my $referencefastafile = shift();
my @samplecallsfiles = ();
while( my $samplecallsfile = shift() ){ push( @samplecallsfiles, $samplecallsfile ); } 
my $statisticsfile = pop( @samplecallsfiles );
my $bestsnpfastafile = pop( @samplecallsfiles );
my $snpfastafile = pop( @samplecallsfiles );
my $allcallmatrixfile = pop( @samplecallsfiles );
my $bestsnpmatrixfile = pop( @samplecallsfiles );
my $outputmatrixfile = pop( @samplecallsfiles );

# Parse extra data from filenames, if present
my @vcffiles = ();
my @externalfiles = ();
my $samplefileinfo = {};
my $duplicatesfile = '';
foreach my $samplecallsfile (@samplecallsfiles)
{
  if( $samplecallsfile =~ /^((?:[A-Za-z0-9._-]+,)+)::(.*)$/ )
  {
    my $filepath = $2;
    my @fileinfo = split( /,/, $1 );
    if( $fileinfo[0] eq 'external' )
    {
      push( @externalfiles, $filepath );
      $samplefileinfo->{$filepath} = "external," . $fileinfo[1];
    }
    elsif( $fileinfo[0] eq 'vcf' )
    {
      push( @vcffiles, $filepath );
      $samplefileinfo->{$filepath} = $fileinfo[1] . "," . $fileinfo[2];
    }
    elsif( $fileinfo[0] eq 'dups' ){ $duplicatesfile = $filepath; }
  } else
  {
    push( @vcffiles, $samplecallsfile );
    $samplefileinfo->{$samplecallsfile} = "pre-aligned,pre-called";
  }
}

my $samplecallsdata = {};
my $seenchromosomes = {};
my $referencecalls = {};
my $indellist = {};
my $dupscalls = {};
my $samplefilterdata = {};
$samplefilterdata->{'called'} = {};
$samplefilterdata->{'coverage'} = {};
$samplefilterdata->{'proportion'} = {};

# Read in reference file
if( open( my $referencefilehandle, '<', $referencefastafile ) )
{
  my $linefromfile = "";
  my $currentchromosome = "";
  while( $linefromfile = <$referencefilehandle> )
  {
    if( $linefromfile =~ /^>([^\s]+)(?:\s|$)/ )
    {
      $currentchromosome = $1;
      $seenchromosomes->{$currentchromosome} = 1;
      $referencecalls->{$currentchromosome} = '';
    }
    elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([A-Za-z.-]+)\s*$/ ) ){ $referencecalls->{$currentchromosome} .= $1; }
  }
  close( $referencefilehandle );
} else { print STDERR "Could not open '$referencefastafile'!\n"; }

# Read in duplicate region calls file
if( length( $duplicatesfile ) && open( my $duplicateshandle, '<', $duplicatesfile ) )
{
  my $linefromfile = "";
  my $currentchromosome = "";
  while( $linefromfile = <$duplicateshandle> )
  {
    if( $linefromfile =~ /^>([^\s]+)(?:\s|$)/ )
    {
      $currentchromosome = $1;
      $dupscalls->{$currentchromosome} = '';
    }
    elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([01]+)\s*$/ ) ){ $dupscalls->{$currentchromosome} .= $1; }
  }
  close( $duplicateshandle );
} else { print STDERR "Could not open '$duplicatesfile'!\n"; }

# Read in external genomes
foreach my $externalfile (@externalfiles)
{
  if( $externalfile =~ /^(?:.*\/)?([^\/]+?)\.frankenfasta$/ )
  {
    my $externalnickname = $1;
    $samplecallsdata->{$externalfile} = {};
    $samplecallsdata->{$externalfile}{$externalnickname} = {};
    if( open( my $externalfilehandle, '<', $externalfile ) )
    {
      my $linefromfile = "";
      my $currentchromosome = "";
      while( $linefromfile = <$externalfilehandle> )
      {
        if( $linefromfile =~ /^>franken::(.*[^\s])\s*$/ )
        {
          $currentchromosome = $1;
          $samplecallsdata->{$externalfile}{$externalnickname}{$currentchromosome} = '';
        }
        elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([A-Za-z.-]+)\s*$/ ) ){ $samplecallsdata->{$externalfile}{$externalnickname}{$currentchromosome} .= $1; }
      }
      close( $externalfilehandle );
    } else { print STDERR "Could not open '$externalfile'!\n"; }
  }
}

# Read in called SNPs
my $vcfcallcounts = {};
$vcfcallcounts->{'depthsum'} = {};
$vcfcallcounts->{'breadthpositions'} = {};
$vcfcallcounts->{'depthfiltered'} = {};
$vcfcallcounts->{'propfiltered'} = {};
foreach my $samplecallsfile ( @vcffiles )
{
  if( -e( $samplecallsfile ) && ( my $samplefilehandle = eval { Vcf->new( 'file'=>$samplecallsfile ); } ) )
  {
    eval { $samplefilehandle->parse_header(); };
    my (@samplelist) = eval { $samplefilehandle->get_samples(); };
    $samplecallsdata->{$samplecallsfile} = {};
    $samplefilterdata->{'called'}{$samplecallsfile} = {};
    $samplefilterdata->{'coverage'}{$samplecallsfile} = {};
    $samplefilterdata->{'proportion'}{$samplecallsfile} = {};
    foreach my $currentsample (@samplelist)
    {
      $samplecallsdata->{$samplecallsfile}{$currentsample} = {};
      $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample} = {};
      $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample} = {};
      $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample} = {};
    }
    my $positiondata = {};
    while( $positiondata = eval { $samplefilehandle->next_data_hash(); } )
    {
      #print Dumper( $positiondata );
      my $currentchromosome = $positiondata->{'CHROM'};
      if( defined( $referencecalls->{$currentchromosome} ) && ( length( $referencecalls->{$currentchromosome} ) >= ( $positiondata->{'POS'} + length( $positiondata->{'REF'} ) - 1 ) ) ) # Ignore positions outside the reference
      {
        my $i = 0;
        while( $i < length( $positiondata->{'REF'} ) ) # Delections may be called as insertions into reference
        {
          if( ( uc( substr( $referencecalls->{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1 ) ) ne 'N' ) && ( uc( substr( $referencecalls->{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1 ) ) ne uc( substr( $positiondata->{'REF'}, $i, 1 ) ) ) )
          {
            if( ( substr( $referencecalls->{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1 ) !~ /^[TtUu]$/ ) || ( substr( $positiondata->{'REF'}, $i, 1 ) !~ /^[TtUu]$/ ) )
            {
              print STDERR "Are you crazy?! This data was generated from disagreeing references! '" . substr( $referencecalls->{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), length( $positiondata->{'REF'} ) ) . "' is not '$positiondata->{'REF'}' at position $positiondata->{'POS'} on chromosome '$currentchromosome' from file '$samplecallsfile'\n";
              die( "Are you crazy?! This data was generated from disagreeing references! '" . substr( $referencecalls->{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), length( $positiondata->{'REF'} ) ) . "' is not '$positiondata->{'REF'}' at position $positiondata->{'POS'} on chromosome '$currentchromosome' from file '$samplecallsfile'\n" );
            }
          }
          $i++;
        }
        foreach my $currentsample (@samplelist)
        {
          if( !defined( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome} = ''; }
          if( !defined( $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome} = ''; }
          if( !defined( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome} = ''; }
          if( !defined( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome} = ''; }
          if( length( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome} ) < ( $positiondata->{'POS'} - 1 ) ){ $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome} .= ( 'X' x ( $positiondata->{'POS'} - length( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome} ) - 1 ) ); }
          if( length( $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) < ( $positiondata->{'POS'} - 1 ) ){ $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome} .= ( 'N' x ( $positiondata->{'POS'} - length( $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) - 1 ) ); }
          if( length( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) < ( $positiondata->{'POS'} - 1 ) ){ $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome} .= ( '?' x ( $positiondata->{'POS'} - length( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) - 1 ) ); }
          if( length( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) < ( $positiondata->{'POS'} - 1 ) ){ $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome} .= ( '?' x ( $positiondata->{'POS'} - length( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome} ) - 1 ) ); }
          my $calledallele;
          my $genotypenum = $positiondata->{'gtypes'}{$currentsample}{'GT'};
          if( $genotypenum =~ /^(\d+)[\/|]/ ){ $genotypenum = $1; } # For SNP callers that brokenly call everything diploid
          if( $genotypenum =~ /^\d{1,3}$/ )
          {
            if( ( $positiondata->{'ALT'}[0] eq '.' ) || ( $genotypenum == 0 ) ){ $calledallele = $positiondata->{'REF'}; }
            else { $calledallele = $positiondata->{'ALT'}[ $genotypenum - 1 ]; }
          }
          elsif( $genotypenum =~ /^\./ ){ $calledallele = ''; }
          else { print STDERR "Could not understand genotype number '$positiondata->{'gtypes'}{$currentsample}{'GT'}' at position $positiondata->{'POS'} on chromosome '$currentchromosome' from file '$samplecallsfile'\n"; }
          my $positiondepth = ( defined( $positiondata->{'INFO'}{'DP'} ) ? ( $positiondata->{'INFO'}{'DP'} / scalar( @samplelist ) ) : ( defined( $positiondata->{'INFO'}{'ADP'} ) ? ( $positiondata->{'INFO'}{'ADP'} / scalar( @samplelist ) ) : 0 ) );
          if( defined( $positiondata->{'gtypes'}{$currentsample}{'DP'} ) ){ $positiondepth = $positiondata->{'gtypes'}{$currentsample}{'DP'}; }
          my $callfrequency = -1;
          if( defined( $positiondata->{'gtypes'}{$currentsample}{'AD'} ) )
          {
            if( $positiondata->{'gtypes'}{$currentsample}{'AD'} =~ /,/ ){ $callfrequency = ( split( /,/, $positiondata->{'gtypes'}{$currentsample}{'AD'} ) )[$genotypenum]; }
            elsif( $genotypenum ){ $callfrequency = $positiondata->{'gtypes'}{$currentsample}{'AD'}; }
          }
          if( defined( $positiondata->{'gtypes'}{$currentsample}{'RD'} ) && ( $genotypenum == 0 ) ){ $callfrequency = $positiondata->{'gtypes'}{$currentsample}{'RD'}; }
          if( defined( $positiondata->{'INFO'}{'PL'} ) )
          {
            my $callfrequencypl = scalar( grep( ( uc( $_ ) eq uc( $calledallele ) ), split( '', $positiondata->{'INFO'}{'PL'} ) ) );
            if( ( $callfrequency == -1 ) || ( $callfrequencypl < $callfrequency ) ){ $callfrequency = $callfrequencypl; }
          }
          if( !defined( $vcfcallcounts->{'depthsum'}{$samplecallsfile} ) ){ $vcfcallcounts->{'depthsum'}{$samplecallsfile} = {}; }
          if( !defined( $vcfcallcounts->{'depthsum'}{$samplecallsfile}{$currentsample} ) ){ $vcfcallcounts->{'depthsum'}{$samplecallsfile}{$currentsample} = 0; }
          if( !defined( $vcfcallcounts->{'breadthpositions'}{$samplecallsfile} ) ){ $vcfcallcounts->{'breadthpositions'}{$samplecallsfile} = {}; }
          if( !defined( $vcfcallcounts->{'breadthpositions'}{$samplecallsfile}{$currentsample} ) ){ $vcfcallcounts->{'breadthpositions'}{$samplecallsfile}{$currentsample} = 0; }
          if( !defined( $vcfcallcounts->{'depthfiltered'}{$samplecallsfile} ) ){ $vcfcallcounts->{'depthfiltered'}{$samplecallsfile} = {}; }
          if( !defined( $vcfcallcounts->{'depthfiltered'}{$samplecallsfile}{$currentsample} ) ){ $vcfcallcounts->{'depthfiltered'}{$samplecallsfile}{$currentsample} = 0; }
          if( !defined( $vcfcallcounts->{'propfiltered'}{$samplecallsfile} ) ){ $vcfcallcounts->{'propfiltered'}{$samplecallsfile} = {}; }
          if( !defined( $vcfcallcounts->{'propfiltered'}{$samplecallsfile}{$currentsample} ) ){ $vcfcallcounts->{'propfiltered'}{$samplecallsfile}{$currentsample} = 0; }
          if( $positiondepth > 0 ){ $vcfcallcounts->{'depthsum'}{$samplecallsfile}{$currentsample} += $positiondepth; }
          if( $positiondepth < $mincoverage )
          {
            $vcfcallcounts->{'depthfiltered'}{$samplecallsfile}{$currentsample}++;
            substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' );
          } elsif( ( $callfrequency != -1 ) && ( $callfrequency < $mincoverage ) )
          {
            $vcfcallcounts->{'depthfiltered'}{$samplecallsfile}{$currentsample}++;
            substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' );
          } else
          {
            $vcfcallcounts->{'breadthpositions'}{$samplecallsfile}{$currentsample}++;
            substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'Y' );
          }
          if( ( $callfrequency != -1 ) && ( ( $callfrequency / $positiondepth ) < $minproportion ) )
          {
            $vcfcallcounts->{'propfiltered'}{$samplecallsfile}{$currentsample}++;
            substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' );
          } else { substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'Y' ); }
          if( length( $positiondata->{'REF'} ) > 1 ) # There's a delete
          {
            $calledallele = substr( $positiondata->{'REF'}, 0, 1 );
            substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, $calledallele );
            if( !defined( $indellist->{$samplecallsfile} ) ){ $indellist->{$samplecallsfile} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample} ) ){ $indellist->{$samplecallsfile}{$currentsample} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} = {}; }
            $i = 1;
            if( ( uc( substr( $positiondata->{'REF'}, 0, 1 ) ) eq uc( substr( $calledallele, 0, 1 ) ) ) && ( uc( substr( $positiondata->{'REF'}, 1 ) ) ne uc( $calledallele ) ) ) # The format of the delete makes sense
            {
              my $j = 1;
              while( $i < length( $positiondata->{'REF'} ) )
              {
                if( uc( substr( $positiondata->{'REF'}, $i, 1 ) ) eq uc( substr( $calledallele, $j, 1 ) ) )
                {
                  substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, substr( $calledallele, $j, 1 ) );
                  $j++;
                } else
                {
                  $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome}{( $positiondata->{'POS'} + $i )} = '.';
                  substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, '.' );
                }
                $i++;
              }
            } else # The format of the delete does not make sense, blindly deleting the entire range
            {
              while( $i < length( $positiondata->{'REF'} ) )
              {
                $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome}{( $positiondata->{'POS'} + $i )} = '.';
                substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, '.' );
                $i++;
              }
            }
          } elsif( length( $calledallele ) > 1 ) # There's an insert
          {
            if( !defined( $indellist->{$samplecallsfile} ) ){ $indellist->{$samplecallsfile} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample} ) ){ $indellist->{$samplecallsfile}{$currentsample} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} = {}; }
            $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome}{$positiondata->{'POS'}} = $calledallele;
            $calledallele = substr( $calledallele, 0, 1 );
          } elsif( length( $calledallele ) == 0 ) # There's a delete
          {
            if( !defined( $indellist->{$samplecallsfile} ) ){ $indellist->{$samplecallsfile} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample} ) ){ $indellist->{$samplecallsfile}{$currentsample} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} = {}; }
            $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome}{$positiondata->{'POS'}} = '.';
            $calledallele = '.';
          }
          substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, $calledallele );
          substr( $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'Y' );
        }
      } else { print STDERR "Call omitted for being outside of reference at position $positiondata->{'POS'} on chromosome '$currentchromosome from file '$samplecallsfile'.\n"; }
    }
    eval { $samplefilehandle->close(); };
  } else { print STDERR "Could not open '$samplecallsfile'!\n"; }
}

#print Dumper( $samplecallsdata );

# Make SNP matrix
my $snpfastadata = {};
my $bestsnpfastadata = {};
my $refsnpfastadata = '';
my $bestrefsnpfastadata = '';
my $coregenomesize = 0;
my $numtotalvariant = 0;
my $numbestsnps = 0;
my $badcallcounts = {};
$badcallcounts->{'N'} = {};
$badcallcounts->{'X'} = {};
$badcallcounts->{'?'} = {};
my $badcallsum = {};
$badcallsum->{'N'} = 0;
$badcallsum->{'X'} = 0;
$badcallsum->{'?'} = 0;
if( open( my $matrixfilehandle, '>', $outputmatrixfile ) && open( my $bestsnpfilehandle, '>', $bestsnpmatrixfile ) && open( my $allcallfilehandle, '>', $allcallmatrixfile ) )
{
  my $allpatternarray = {};
  my $allnextpatternnum = 1;
  my $snppatternarray = {};
  my $snpnextpatternnum = 1;
  my $bestsnppatternarray = {};
  my $bestsnpnextpatternnum = 1;
  print $matrixfilehandle "LocusID\tReference\t";
  print $bestsnpfilehandle "LocusID\tReference\t";
  print $allcallfilehandle "LocusID\tReference\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) )
  {
    $badcallcounts->{'N'}{$samplefile} = {};
    $badcallcounts->{'X'}{$samplefile} = {};
    $badcallcounts->{'?'}{$samplefile} = {};
    foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
    {
      print $matrixfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      print $bestsnpfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      print $allcallfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      $badcallcounts->{'N'}{$samplefile}{$samplecolumn} = 0;
      $badcallcounts->{'X'}{$samplefile}{$samplecolumn} = 0;
      $badcallcounts->{'?'}{$samplefile}{$samplecolumn} = 0;
    }
  }
  print $matrixfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tInDupRegion\tNotes\t\n";
  print $bestsnpfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tInDupRegion\tNotes\t\n";
  print $allcallfilehandle "#SNPcall\t#Refcall\t#A\t#C\t#G\t#T\t#NXindel\t%SNPcall\t%Refcall\t%A\t%C\t%G\t%T\t%NXindel\tChromosome\tPosition\tPattern\tPattern#\tInDupRegion\tNotes\t\n";
  foreach my $currentchromosome ( sort keys( %$seenchromosomes ) )
  {
    my $currentposition = 1;
    while( $currentposition <= length( $referencecalls->{$currentchromosome} ) )
    {
      my $linetoprint = '';
      my $allelecounts = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0, 'ref' => 0, 'snp' => 0, 'sum' => 0 };
      my $referencecall = substr( $referencecalls->{$currentchromosome}, ( $currentposition - 1 ), 1 );
      my $patternassignments = {};
      my $callpattern = '';
      my $simplifiedrefcall = "N";
      if( $referencecall =~ /^[ACGTUacgtu]$/ )
      {
        $simplifiedrefcall = uc( $referencecall );
        if( $referencecall =~ /^[Uu]$/ ){ $simplifiedrefcall = "T"; }
        $patternassignments->{$referencecall} = 1;
        $callpattern .= '1';
      } else { $callpattern .= 'N'; }
      my $nextassignmentnum = 2;
      my $dupscall = ( defined( $dupscalls->{$currentchromosome} ) ? ( ( substr( $dupscalls->{$currentchromosome}, ( $currentposition - 1 ), 1 ) ) ? 'Yes' : 'No' ) : 'Unchecked' );
      my $intersectioncalls = {};
      $linetoprint .= "${currentchromosome}::${currentposition}\t$referencecall\t";
      my $somethingcalledn = 0;
      my $somethingcalledx = 0;
      my $somethingcalledoth = 0;
      foreach my $samplefile ( sort keys( $samplecallsdata ) )
      {
        foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
        {
          if( length( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} ) < ( $currentposition ) ){ $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} .= ( 'X' x ( $currentposition - length( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} ) ) ); }
          my $samplecall = substr( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 );
          if( defined( $indellist->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} ) ){ $samplecall = $indellist->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}; }
          $linetoprint .= "$samplecall\t";
          my $simplifiedbasecall = "N";
          if( substr( $samplecall, 0, 1 ) =~ /^[ACGTUacgtu]$/ )
          {
            $simplifiedbasecall = uc( $samplecall );
            if( $samplecall =~ /^[Uu]$/ ){ $simplifiedbasecall = "T"; }
          }
          $allelecounts->{$simplifiedbasecall}++;
          if( ( $simplifiedrefcall ne 'N' ) && ( $simplifiedbasecall ne 'N' ) ){ if( $simplifiedrefcall eq $simplifiedbasecall ){ $allelecounts->{"ref"}++; } else { $allelecounts->{"snp"}++; } }
          $allelecounts->{"sum"}++;
          if( $simplifiedbasecall eq "N" ){ $callpattern .= "N"; } else
          {
            if( !defined( $patternassignments->{$simplifiedbasecall} ) )
            {
              $patternassignments->{$simplifiedbasecall} = $nextassignmentnum;
              if( $nextassignmentnum ne '+' ){ $nextassignmentnum = ( ( $nextassignmentnum < 9 ) ? ( $nextassignmentnum + 1 ) : '+' ); }
            }
            $callpattern .= $patternassignments->{$simplifiedbasecall};
            if( defined( $intersectioncalls->{$samplecolumn} ) )
            {
              if( ( $intersectioncalls->{$samplecolumn} eq 'N' ) || ( $simplifiedbasecall eq 'N' ) || ( $simplifiedbasecall ne $intersectioncalls->{$samplecolumn} ) ){ $intersectioncalls->{$samplecolumn} = 'N'; }
            } else { $intersectioncalls->{$samplecolumn} = ( ( ( $simplifiedrefcall ne 'N' ) && ( $simplifiedbasecall ne 'N' ) ) ? $simplifiedbasecall : 'N' ); }
          }
          if( uc( $samplecall ) eq 'N' )
          {
            $badcallcounts->{'N'}{$samplefile}{$samplecolumn}++;
            $somethingcalledn = 1;
          }
          elsif( uc( $samplecall ) eq 'X' )
          {
            $badcallcounts->{'X'}{$samplefile}{$samplecolumn}++;
            $somethingcalledx = 1;
          }
          elsif( $simplifiedbasecall eq 'N' )
          {
            $badcallcounts->{'?'}{$samplefile}{$samplecolumn}++;
            $somethingcalledoth = 1;
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
      if( ( $allelecounts->{'snp'} > 0 ) && ( $dupscall ne 'Yes' ) )
      {
        if( !defined( $snppatternarray->{$callpattern} ) ){ $snppatternarray->{$callpattern} = $snpnextpatternnum++; }
        print $matrixfilehandle $linetoprint . "'$callpattern'\t$snppatternarray->{$callpattern}\t$dupscall\t\t\n";
        $refsnpfastadata .= $simplifiedrefcall;
        foreach my $samplefile ( sort keys( $samplecallsdata ) )
        {
          foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
          {
            if( !defined( $snpfastadata->{$samplefile} ) ){ $snpfastadata->{$samplefile} = {}; }
            if( !defined( $snpfastadata->{$samplefile}{$samplecolumn} ) ){ $snpfastadata->{$samplefile}{$samplecolumn} = ''; }
            $snpfastadata->{$samplefile}{$samplecolumn} .= substr( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 );
          }
        }
        $numtotalvariant++;
        my $intersectionpresent = 1;
        foreach my $samplename ( keys( $intersectioncalls ) ){ if( $intersectioncalls->{$samplename} eq 'N' ){ $intersectionpresent = 0; } }
        if( ( $allelecounts->{"N"} == 0 ) && ( $intersectionpresent == 1 ) )
        {
          if( !defined( $bestsnppatternarray->{$callpattern} ) ){ $bestsnppatternarray->{$callpattern} = $bestsnpnextpatternnum++; }
          print $bestsnpfilehandle $linetoprint . "'$callpattern'\t$bestsnppatternarray->{$callpattern}\t$dupscall\t\t\n";
          $bestrefsnpfastadata .= $simplifiedrefcall;
          foreach my $samplefile ( sort keys( $samplecallsdata ) )
          {
             foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
            {
              if( !defined( $bestsnpfastadata->{$samplefile} ) ){ $bestsnpfastadata->{$samplefile} = {}; }
              if( !defined( $bestsnpfastadata->{$samplefile}{$samplecolumn} ) ){ $bestsnpfastadata->{$samplefile}{$samplecolumn} = ''; }
              $bestsnpfastadata->{$samplefile}{$samplecolumn} .= substr( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 );
            }
          }
          $numbestsnps++;
        }
      }
      if( !defined( $allpatternarray->{$callpattern} ) ){ $allpatternarray->{$callpattern} = $allnextpatternnum++; }
      print $allcallfilehandle $linetoprint . "'$callpattern'\t$allpatternarray->{$callpattern}\t$dupscall\t\t\n";
      if( ( $allelecounts->{"N"} == 0 ) && ( $dupscall ne 'Yes' ) ){ $coregenomesize++; }
      if( $somethingcalledn ){ $badcallsum->{'N'}++; }
      if( $somethingcalledx ){ $badcallsum->{'X'}++; }
      if( $somethingcalledoth ){ $badcallsum->{'?'}++; }
      $currentposition++;
    }
  }
  close( $matrixfilehandle );
  close( $bestsnpfilehandle );
  close( $allcallfilehandle );
} else { print STDERR "Could not open '$outputmatrixfile'!\n"; }

if( open( my $snpfastahandle, '>', $snpfastafile ) && open( my $bestsnpfastahandle, '>', $bestsnpfastafile ) )
{
  if( length( $refsnpfastadata ) )
  {
    print $snpfastahandle ">SNP::Reference\n";
    print $snpfastahandle "$refsnpfastadata\n";
    foreach my $samplefile ( sort keys( $samplecallsdata ) )
    {
      foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
      {
        print $snpfastahandle ">SNP::${samplecolumn}::$samplefileinfo->{$samplefile}\n";
        print $snpfastahandle "$snpfastadata->{$samplefile}{$samplecolumn}\n";
      }
    }
  }
  if( length( $bestrefsnpfastadata ) )
  {
    print $bestsnpfastahandle ">SNP::Reference\n";
    print $bestsnpfastahandle "$bestrefsnpfastadata\n";
    foreach my $samplefile ( sort keys( $samplecallsdata ) )
    {
      foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
      {
        print $bestsnpfastahandle ">SNP::${samplecolumn}::$samplefileinfo->{$samplefile}\n";
        print $bestsnpfastahandle "$bestsnpfastadata->{$samplefile}{$samplecolumn}\n";
      }
    }
  }
  close( $snpfastahandle );
  close( $bestsnpfastahandle );
} else { print STDERR "Could not open '$snpfastafile'!\n"; }

# Do final calculations and write the statistics file
if( open( my $statshandle, '>', $statisticsfile ) )
{
  print $statshandle "# This is a first iteration of a proposal to come up with an eventual plan to outline the future possibility of looking into an investigation whose goal might be to someday create a standard format for this file.\n";
  print $statshandle "# It's especially complicated because it needs to be both human-readable and script-parsable, and preferably have enough information for people to understand what all them crazy numbers mean.\n";
  print $statshandle "Stat_ID\tTotal\t";
  my $nopersampleinfo = '';
  foreach my $samplefile ( sort keys( $samplecallsdata ) )
  {
    foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
    {
      $nopersampleinfo .= "\t";
      print $statshandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
    }
  }
  print $statshandle "\n";
  my $referencelength = 0;
  foreach my $currentchromosome ( keys( $seenchromosomes ) ){ $referencelength += length( $referencecalls->{$currentchromosome} ); }
  print $statshandle "reference_length\t$referencelength\t$nopersampleinfo\n";
  my $referencegood = 0;
  #foreach my $currentchromosome ( keys( $seenchromosomes ) ){ $referencegood += scalar( () = ( $referencecalls->{$currentchromosome} =~ /[ACGTUacgtu]/g ) ); }
  foreach my $currentchromosome ( keys( $seenchromosomes ) ){ $referencegood += scalar( grep( /[ACGTUacgtu]/, split( '', $referencecalls->{$currentchromosome} ) ) ); }
  print $statshandle "reference_clean\t$referencegood\t$nopersampleinfo\n";
  my $duppositions = 0;
  foreach my $currentchromosome ( keys( $seenchromosomes ) ){ $duppositions += scalar( grep( /1/, split( '', $dupscalls->{$currentchromosome} ) ) ); }
  print $statshandle "dups_count\t$duppositions\t$nopersampleinfo\n";
  print $statshandle "dups_proportion\t" . sprintf( "%.2f", ( $duppositions * 100 / $referencelength ) ) . "%\t$nopersampleinfo\n";
  print $statshandle "core_genome_size\t$coregenomesize\t$nopersampleinfo\n";
  print $statshandle "core_genome_coverage\t" . sprintf( "%.2f", ( $coregenomesize * 100 / $referencelength ) ) . "%\t$nopersampleinfo\n";
  print $statshandle "num_total_variant\t$numtotalvariant\t$nopersampleinfo\n";
  print $statshandle "num_best_snps\t$numbestsnps\t$nopersampleinfo\n";
  print $statshandle "core_genome_portion_snp\t" . ( ( $coregenomesize > 0 ) ? ( sprintf( "%.2f", ( $numbestsnps * 100 / $coregenomesize ) ) . "%\t" ) : '-\t' ) . "$nopersampleinfo\n";
  print $statshandle "num_n\t$badcallsum->{'N'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "$badcallcounts->{'N'}{$samplefile}{$samplecolumn}\t"; } }
  print $statshandle "\n";
  print $statshandle "num_x\t$badcallsum->{'X'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "$badcallcounts->{'X'}{$samplefile}{$samplecolumn}\t"; } }
  print $statshandle "\n";
  print $statshandle "num_other\t$badcallsum->{'?'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "$badcallcounts->{'?'}{$samplefile}{$samplecolumn}\t"; } }
  print $statshandle "\n";
  print $statshandle "depth_filtered\t-\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle '' . ( defined( $vcfcallcounts->{'depthfiltered'}{$samplefile}{$samplecolumn} ) ? "$vcfcallcounts->{'depthfiltered'}{$samplefile}{$samplecolumn}\t" : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "proportion_filtered\t-\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle '' . ( defined( $vcfcallcounts->{'propfiltered'}{$samplefile}{$samplecolumn} ) ? "$vcfcallcounts->{'propfiltered'}{$samplefile}{$samplecolumn}\t" : "-\t" ); } }
  print $statshandle "\n";
  my $avgavg = 0;
  my $avgcount = 0;
  my $avgstring = '';
  foreach my $samplefile ( sort keys( $samplecallsdata ) )
  {
    foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
    {
      if( defined( $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn} ) )
      {
        $avgcount++;
        $avgavg = ( ( $avgavg * ( $avgcount - 1 ) + $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn} ) / $avgcount );
        $avgstring .= sprintf( "%.2f", ( $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn} / $referencelength ) ) . "\t";
      } else { $avgstring .= "-\t"; }
    }
  }
  print $statshandle "coverage_depth\t" . sprintf( "%.2f", ( $avgavg / $referencelength ) ) . "\t$avgstring\n";
  print $statshandle "coverage_breadth\t-\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle '' . ( defined( $vcfcallcounts->{'breadthpositions'}{$samplefile}{$samplecolumn} ) ? sprintf( "%.2f", ( $vcfcallcounts->{'breadthpositions'}{$samplefile}{$samplecolumn} * 100 / $referencelength ) ) . "%\t" : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "# reference_length: total number of positions found across all chromosomes in the reference\n";
  print $statshandle "# reference_clean: number of positions called A/C/G/T in the reference\n";
  print $statshandle "# dups_count: number of positions that appear to be in duplicated regions of the reference\n";
  print $statshandle "# dups_proportion: dups_count / reference_length * 100%\n";
  print $statshandle "# core_genome_size: number of positions called A/C/G/T across all analyses of all samples\n";
  print $statshandle "# core_genome_coverage: core_genome_size / reference_length * 100%\n";
  print $statshandle "# num_total_variant: number of positions where the reference and at least one analysis of at least one sample differed, were both called A/C/G/T, and were outside known-duplicated regions\n";
  print $statshandle "# num_best_snps: number of positions where the requirements for num_total_variant were met, all analyses of all samples were called A/C/G/T, and all analyses of any one sample were in agreement\n";
  print $statshandle "# core_genome_portion_snp: num_best_snps / core_genome_size * 100%\n";
  print $statshandle "# num_n: number of positions where an analysis tool called an N or the filter cutoffs were not all met\n";
  print $statshandle "# num_x: number of positions where data was missing or an analysis tool made no call or called an X\n";
  print $statshandle "# num_other: number of positions where an insertion or degenerate base was called\n";
  print $statshandle "# depth_filtered: number of positions whose call was changed to N by the depth filter\n";
  print $statshandle "# proportion_filtered: number of positions whose call was changed to N by the proportion filter\n";
  print $statshandle "# coverage_depth: the sum of the depths at each position, whether filtered or not, divided by reference_length\n";
  print $statshandle "# coverage_breadth: number of positions that had any call made and passed depth filtering, divided by reference_length * 100%\n";
} else { print STDERR "Could not open '$statisticsfile'!\n"; }


