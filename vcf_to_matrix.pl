#!/usr/bin/perl -w
#/usr/bin/perl -w -d:SmallProf

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use Vcf;
#use Data::Dumper;

if( @ARGV < 9 )
{
  print <<EOF;
Meant to be called from the pipeline automatically.

Usage:
vcf_to_matrix.pl <min_coverage> <min_proportion> <reference.fasta> <samples.vcf> [samples2.vcf ...] <allvariant_matrix.tsv> <bestsnps_matrix.tsv> <allcallable_matrix.tsv> <bestsnps.snpfasta> <statistics_file.tsv>
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
my $allcallmatrixfile = pop( @samplecallsfiles );
my $bestsnpmatrixfile = pop( @samplecallsfiles );
my $outputmatrixfile = pop( @samplecallsfiles );

# Parse extra data from filenames, if present
my $vcffiles = {};
my $externalfiles = {};
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
      $externalfiles->{$filepath} = 1;
      $samplefileinfo->{$filepath} = "external," . $fileinfo[1];
    }
    elsif( $fileinfo[0] eq 'vcf' )
    {
      $vcffiles->{$filepath} = 1;
      $samplefileinfo->{$filepath} = $fileinfo[1] . "," . $fileinfo[2];
    }
    elsif( $fileinfo[0] eq 'dups' ){ $duplicatesfile = $filepath; }
  } else
  {
    $vcffiles->{$samplecallsfile} = 1;
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
my $statscounters = {};

# Read in reference file
if( open( my $referencefilehandle, '<', $referencefastafile ) )
{
  my $linefromfile = "";
  my $currentchromosome = "";
  $statscounters->{'reflength'} = 0;
  $statscounters->{'refclean'} = 0;
  while( $linefromfile = <$referencefilehandle> )
  {
    if( $linefromfile =~ /^>([^\s]+)(?:\s|$)/ )
    {
      $currentchromosome = $1;
      $seenchromosomes->{$currentchromosome} = 1;
      $referencecalls->{$currentchromosome} = '';
    }
    elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([A-Za-z.-]+)\s*$/ ) )
    {
      my $refchunk = $1;
      $referencecalls->{$currentchromosome} .= $refchunk;
      $statscounters->{'reflength'} += length( $refchunk );
      $statscounters->{'refclean'} += scalar( grep( /[ACGTUacgtu]/, split( '', $refchunk ) ) );
    }
  }
  close( $referencefilehandle );
} else { print STDERR "Could not open '$referencefastafile'!\n"; }

# Read in duplicate region calls file
if( length( $duplicatesfile ) && open( my $duplicateshandle, '<', $duplicatesfile ) )
{
  my $linefromfile = "";
  my $currentchromosome = "";
  $statscounters->{'duppositions'} = 0;
  while( $linefromfile = <$duplicateshandle> )
  {
    if( $linefromfile =~ /^>([^\s]+)(?:\s|$)/ )
    {
      $currentchromosome = $1;
      $dupscalls->{$currentchromosome} = '';
    }
    elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([01]+)\s*$/ ) )
    {
      my $dupchunk = $1;
      $dupscalls->{$currentchromosome} .= $dupchunk;
      $statscounters->{'duppositions'} += scalar( grep( /1/, split( '', $dupchunk ) ) );
    }
  }
  close( $duplicateshandle );
} else { print STDERR "Could not open '$duplicatesfile'!\n"; }

# Read in external genomes
foreach my $externalfile ( sort keys( $externalfiles ) )
{
  if( $externalfile =~ /^(?:.*\/)?([^\/]+?)\.frankenfasta$/ )
  {
    my $externalnickname = $1;
    $samplecallsdata->{$externalfile} = {};
    $samplecallsdata->{$externalfile}{$externalnickname} = {};
    $samplefilterdata->{'called'}{$externalfile} = {};
    $samplefilterdata->{'called'}{$externalfile}{$externalnickname} = {};
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
          $samplefilterdata->{'called'}{$externalfile}{$externalnickname}{$currentchromosome} = '';
        }
        elsif( length( $currentchromosome ) && ( $linefromfile =~ /^([A-Za-z.-]+)\s*$/ ) )
        {
          my $externalchunk = $1;
          $samplecallsdata->{$externalfile}{$externalnickname}{$currentchromosome} .= $externalchunk;
          $externalchunk =~ s/[^N]/Y/g;
          $samplefilterdata->{'called'}{$externalfile}{$externalnickname}{$currentchromosome} .= $externalchunk;
        }
      }
      close( $externalfilehandle );
    } else { print STDERR "Could not open '$externalfile'!\n"; }
  }
}

# Read in called SNPs
my $vcfcallcounts = {};
$vcfcallcounts->{'depthsum'} = {};
foreach my $samplecallsfile ( sort keys( $vcffiles ) )
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
          my $genotypenum = ( defined( $positiondata->{'gtypes'}{$currentsample}{'GT'} ) ? $positiondata->{'gtypes'}{$currentsample}{'GT'} : 0 );
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
          if( $positiondepth > 0 ){ $vcfcallcounts->{'depthsum'}{$samplecallsfile}{$currentsample} += $positiondepth; }
          if( $positiondepth < $mincoverage ){ substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' ); }
          #elsif( ( $callfrequency != -1 ) && ( $callfrequency < $mincoverage ) ){ substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' ); } # This applies the depth-filtering to only reads matching the call
          else { substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'Y' ); }
          if( $callfrequency != -1 )
          {
            if( ( $callfrequency / $positiondepth ) < $minproportion ){ substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'N' ); }
            else { substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, 'Y' ); }
          #} else { substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, ( ( $calledallele eq $positiondata->{'REF'} ) ? '-' : '?' ) ); }
          } else { substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, '-' ); }
          if( length( $positiondata->{'REF'} ) > 1 ) # There's a delete
          {
            $calledallele = substr( $positiondata->{'REF'}, 0, 1 );
            substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, $calledallele );
            if( !defined( $indellist->{$samplecallsfile} ) ){ $indellist->{$samplecallsfile} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample} ) ){ $indellist->{$samplecallsfile}{$currentsample} = {}; }
            if( !defined( $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} ) ){ $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome} = {}; }
            my $coveragefiltercall = substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1 );
            my $proportionfiltercall = substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1 );
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
                substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, $coveragefiltercall );
                substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, $proportionfiltercall );
                $i++;
              }
            } else # The format of the delete does not make sense, blindly deleting the entire range
            {
              while( $i < length( $positiondata->{'REF'} ) )
              {
                $indellist->{$samplecallsfile}{$currentsample}{$currentchromosome}{( $positiondata->{'POS'} + $i )} = '.';
                substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, '.' );
                substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, $coveragefiltercall );
                substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} + $i - 1 ), 1, $proportionfiltercall );
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
            #substr( $samplefilterdata->{'coverage'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, '-' );
            #substr( $samplefilterdata->{'proportion'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, '-' );
            $calledallele = '.';
          }
          substr( $samplecallsdata->{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, $calledallele );
          substr( $samplefilterdata->{'called'}{$samplecallsfile}{$currentsample}{$currentchromosome}, ( $positiondata->{'POS'} - 1 ), 1, ( ( uc( $calledallele ) eq 'N' ) ? 'N' : 'Y' ) );
        }
      } else { print STDERR "Call omitted for being outside of reference at position $positiondata->{'POS'} on chromosome '$currentchromosome from file '$samplecallsfile'.\n"; }
    }
    eval { $samplefilehandle->close(); };
  } else { print STDERR "Could not open '$samplecallsfile'!\n"; }
}

#print Dumper( $samplecallsdata );

# Make SNP matrix
my $bestsnpfastadata = {};
my $bestrefsnpfastadata = '';
if( open( my $matrixfilehandle, '>', $outputmatrixfile ) && open( my $bestsnpfilehandle, '>', $bestsnpmatrixfile ) && open( my $allcallfilehandle, '>', $allcallmatrixfile ) )
{
  my $allpatternarray = {};
  my $allnextpatternnum = 1;
  my $snppatternarray = {};
  my $snpnextpatternnum = 1;
  my $bestsnppatternarray = {};
  my $bestsnpnextpatternnum = 1;
  my $numsamplecolumns = 0;
  $statscounters->{'totalvariant'} = 0;
  $statscounters->{'bestsnps'} = 0;
  $statscounters->{'allcalled'} = 0;
  $statscounters->{'samplecalled'} = {};
  $statscounters->{'allcoveragepass'} = 0;
  $statscounters->{'samplecoveragepass'} = {};
  $statscounters->{'allproportionpass'} = 0;
  $statscounters->{'sampleproportionpass'} = {};
  $statscounters->{'anynxdegen'} = 0;
  $statscounters->{'samplenxdegen'} = {};
  $statscounters->{'allbreadth'} = 0;
  $statscounters->{'samplebreadth'} = {};
  $statscounters->{'consensus'} = 0;
  $statscounters->{'bestpos'} = 0;
  print $matrixfilehandle "LocusID\tReference\t";
  print $bestsnpfilehandle "LocusID\tReference\t";
  print $allcallfilehandle "LocusID\tReference\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) )
  {
    $statscounters->{'samplecalled'}{$samplefile} = {};
    $statscounters->{'samplecoveragepass'}{$samplefile} = {};
    $statscounters->{'sampleproportionpass'}{$samplefile} = {};
    $statscounters->{'samplenxdegen'}{$samplefile} = {};
    $statscounters->{'samplebreadth'}{$samplefile} = {};
    foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
    {
      print $matrixfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      print $bestsnpfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      print $allcallfilehandle "${samplecolumn}::$samplefileinfo->{$samplefile}\t";
      $statscounters->{'samplecalled'}{$samplefile}{$samplecolumn} = 0;
      $statscounters->{'samplecoveragepass'}{$samplefile}{$samplecolumn} = 0;
      $statscounters->{'sampleproportionpass'}{$samplefile}{$samplecolumn} = 0;
      $statscounters->{'samplenxdegen'}{$samplefile}{$samplecolumn} = 0;
      $statscounters->{'samplebreadth'}{$samplefile}{$samplecolumn} = 0;
      $numsamplecolumns++;
    }
  }
  print $matrixfilehandle "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tChromosome\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\tPattern#\t\n";
  print $bestsnpfilehandle "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tChromosome\tPosition\tInDupRegion\tSampleConsensus\tPattern\tPattern#\t\n";
  print $allcallfilehandle "#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tChromosome\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\tPattern#\t\n";
  foreach my $currentchromosome ( sort keys( %$seenchromosomes ) )
  {
    my $currentposition = 1;
    while( $currentposition <= length( $referencecalls->{$currentchromosome} ) )
    {
      my $linetoprint = '';
      my $allelecounts = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'indel' => 0, 'N' => 0, 'refcall' => 0, 'snpcall' => 0, 'indelcall' => 0, 'called' => 0, 'coverage' => 0, 'proportion' => 0 };
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
      my $callmadestring = '';
      my $passedcoveragestring = '';
      my $passedproportionstring = '';
      foreach my $samplefile ( sort keys( $samplecallsdata ) )
      {
        foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
        {
          # Get sample call and calculate pre-filter columns (A/C/G/T/Indel/NXdegen)
          if( length( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} ) < ( $currentposition ) ){ $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} .= ( 'X' x ( $currentposition - length( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome} ) ) ); }
          my $samplecall = substr( $samplecallsdata->{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 );
          my $isindel = 0;
          my $isdelete = 0;
          if( defined( $indellist->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition} ) )
          {
            $samplecall = $indellist->{$samplefile}{$samplecolumn}{$currentchromosome}{$currentposition}; 
            $allelecounts->{'indel'}++;
            $isindel = 1;
            if( substr( $samplecall, 0, 1 ) ne '.' ){ $isdelete = 1; }
          }
          $linetoprint .= "$samplecall\t";
          my $simplifiedbasecall = "N";
          if( $samplecall =~ /^[ACGTUacgtu]/ )
          {
            $simplifiedbasecall = uc( substr( $samplecall, 0, 1 ) );
            if( $simplifiedbasecall eq 'U' ){ $simplifiedbasecall = "T"; }
          }
          if( !( $isdelete ) ){ $allelecounts->{$simplifiedbasecall}++; }
          if( $simplifiedbasecall eq 'N' ){ $statscounters->{'samplenxdegen'}{$samplefile}{$samplecolumn}++; }
          # Check filters
          my $callwasmade = 'N';
          if( defined( $samplefilterdata->{'called'}{$samplefile}{$samplecolumn}{$currentchromosome} ) && ( length( $samplefilterdata->{'called'}{$samplefile}{$samplecolumn}{$currentchromosome} ) >= $currentposition ) ){ $callwasmade = substr( $samplefilterdata->{'called'}{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 ); }
          $callmadestring .= $callwasmade;
          $callwasmade = ( ( $callwasmade eq 'Y' ) ? 1 : 0 );
          $allelecounts->{'called'} += $callwasmade;
          $statscounters->{'samplecalled'}{$samplefile}{$samplecolumn} += $callwasmade;
          my $coveragepassed = '?';
          my $proportionpassed = '?';
          if( defined( $externalfiles->{$samplefile} ) )
          {
            $coveragepassed = '-';
            $proportionpassed = '-';
          } else
          {
            if( defined( $samplefilterdata->{'coverage'}{$samplefile}{$samplecolumn}{$currentchromosome} ) && ( length( $samplefilterdata->{'coverage'}{$samplefile}{$samplecolumn}{$currentchromosome} ) >= $currentposition ) ){ $coveragepassed = substr( $samplefilterdata->{'coverage'}{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 ); }
            if( defined( $samplefilterdata->{'proportion'}{$samplefile}{$samplecolumn}{$currentchromosome} ) && ( length( $samplefilterdata->{'proportion'}{$samplefile}{$samplecolumn}{$currentchromosome} ) >= $currentposition ) ){ $proportionpassed = substr( $samplefilterdata->{'proportion'}{$samplefile}{$samplecolumn}{$currentchromosome}, ( $currentposition - 1 ), 1 ); }
          }
          $passedcoveragestring .= $coveragepassed;
          $passedproportionstring .= $proportionpassed;
          $coveragepassed = ( ( ( $coveragepassed eq 'Y' ) || ( $coveragepassed eq '-' ) ) ? 1 : 0 );
          $proportionpassed = ( ( ( $proportionpassed eq 'Y' ) || ( $proportionpassed eq '-' ) ) ? 1 : 0 );
          $allelecounts->{'coverage'} += $coveragepassed;
          $allelecounts->{'proportion'} += $proportionpassed;
          $statscounters->{'samplecoveragepass'}{$samplefile}{$samplecolumn} += $coveragepassed;
          $statscounters->{'sampleproportionpass'}{$samplefile}{$samplecolumn} += $proportionpassed;
          if( $callwasmade && $coveragepassed){ $statscounters->{'samplebreadth'}{$samplefile}{$samplecolumn}++; }
          # Calculate post-filter columns (SNP/Indel/Ref)
          if( $callwasmade && $coveragepassed && $proportionpassed && ( $simplifiedrefcall ne 'N' ) )
          {
            if( $simplifiedrefcall eq $simplifiedbasecall ){ $allelecounts->{'refcall'}++; }
            elsif( $simplifiedbasecall ne 'N' ){ $allelecounts->{'snpcall'}++; }
            if( $isindel ){ $allelecounts->{'indelcall'}++; }
          } else { $simplifiedbasecall = 'N'; }
          # Make pattern
          if( $simplifiedbasecall eq "N" ){ $callpattern .= "N"; } else
          {
            if( !defined( $patternassignments->{$simplifiedbasecall} ) )
            {
              $patternassignments->{$simplifiedbasecall} = $nextassignmentnum;
              $nextassignmentnum += 1;
            }
            $callpattern .= $patternassignments->{$simplifiedbasecall};
          }
          # Check intersections
          if( $isindel )
          {
            $samplecall = uc( $samplecall );
            $samplecall =~ tr/U/T/;
            if( defined( $intersectioncalls->{$samplecolumn} ) )
            {
              if( $samplecall ne $intersectioncalls->{$samplecolumn} ){ $intersectioncalls->{$samplecolumn} = 'N'; }
            } else { $intersectioncalls->{$samplecolumn} = $samplecall; }
          } else
          {
            if( defined( $intersectioncalls->{$samplecolumn} ) )
            {
              if( ( $intersectioncalls->{$samplecolumn} eq 'N' ) || ( $simplifiedbasecall eq 'N' ) || ( $simplifiedbasecall ne $intersectioncalls->{$samplecolumn} ) ){ $intersectioncalls->{$samplecolumn} = 'N'; }
            } else { $intersectioncalls->{$samplecolumn} = $simplifiedbasecall; }
          }
        }
      }
      $linetoprint .= "$allelecounts->{'snpcall'}\t$allelecounts->{'indelcall'}\t$allelecounts->{'refcall'}\t$allelecounts->{'called'}/$numsamplecolumns\t$allelecounts->{'coverage'}/$numsamplecolumns\t$allelecounts->{'proportion'}/$numsamplecolumns\t$allelecounts->{'A'}\t$allelecounts->{'C'}\t$allelecounts->{'G'}\t$allelecounts->{'T'}\t$allelecounts->{'indel'}\t$allelecounts->{'N'}\t";
      $linetoprint .= "$currentchromosome\t$currentposition\t$dupscall\t";
      my $intersectionpresent = 1;
      foreach my $samplename ( keys( $intersectioncalls ) ){ if( $intersectioncalls->{$samplename} eq 'N' ){ $intersectionpresent = 0; } }
      $linetoprint .= ( $intersectionpresent ? "Yes\t" : "No\t" );
      $statscounters->{'consensus'} += $intersectionpresent;
      my $longlinetoprint = "$callmadestring\t$passedcoveragestring\t$passedproportionstring\t";
      if( ( ( $allelecounts->{'snpcall'} + $allelecounts->{'indelcall'} ) > 0 ) && ( $dupscall ne 'Yes' ) )
      {
        if( !defined( $snppatternarray->{$callpattern} ) ){ $snppatternarray->{$callpattern} = $snpnextpatternnum++; }
        print $matrixfilehandle $linetoprint . $longlinetoprint . "'$callpattern'\t$snppatternarray->{$callpattern}\t\n";
        $statscounters->{'totalvariant'}++;
        if( ( $allelecounts->{'snpcall'} > 0 ) && ( $allelecounts->{'indelcall'} == 0 ) && ( $intersectionpresent == 1 ) && ( $allelecounts->{'called'} == $numsamplecolumns ) && ( $allelecounts->{'coverage'} == $numsamplecolumns ) && ( $allelecounts->{'N'} == 0 ) && ( $allelecounts->{'proportion'} == $numsamplecolumns ) )
        {
          if( !defined( $bestsnppatternarray->{$callpattern} ) ){ $bestsnppatternarray->{$callpattern} = $bestsnpnextpatternnum++; }
          print $bestsnpfilehandle $linetoprint . "'$callpattern'\t$bestsnppatternarray->{$callpattern}\t\n";
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
          $statscounters->{'bestsnps'}++;
        }
      }
      if( !defined( $allpatternarray->{$callpattern} ) ){ $allpatternarray->{$callpattern} = $allnextpatternnum++; }
      print $allcallfilehandle $linetoprint . $longlinetoprint . "'$callpattern'\t$allpatternarray->{$callpattern}\t\n";
      if( $allelecounts->{'called'} == $numsamplecolumns ){ $statscounters->{'allcalled'}++; }
      if( $allelecounts->{'coverage'} == $numsamplecolumns ){ $statscounters->{'allcoveragepass'}++; }
      if( $allelecounts->{'proportion'} == $numsamplecolumns ){ $statscounters->{'allproportionpass'}++; }
      if( $allelecounts->{'N'} > 0 ){ $statscounters->{'anynxdegen'}++; }
      if( ( $allelecounts->{'called'} == $numsamplecolumns ) && ( $allelecounts->{'coverage'} == $numsamplecolumns ) && ( $dupscall ne 'Yes' ) )
      {
        $statscounters->{'allbreadth'}++;
        if( $intersectionpresent && ( $allelecounts->{'N'} == 0 ) && ( $allelecounts->{'proportion'} == $numsamplecolumns ) ){ $statscounters->{'bestpos'}++; }
      }
      $currentposition++;
    }
  }
  close( $matrixfilehandle );
  close( $bestsnpfilehandle );
  close( $allcallfilehandle );
} else { print STDERR "Could not open '$outputmatrixfile'!\n"; }

# Make snpfasta
if( open( my $bestsnpfastahandle, '>', $bestsnpfastafile ) )
{
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
  close( $bestsnpfastahandle );
} else { print STDERR "Could not open '$bestsnpfastafile'!\n"; }

# Do final calculations and write the statistics file
if( open( my $statshandle, '>', $statisticsfile ) )
{
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
  print $statshandle "reference_length\t$statscounters->{'reflength'}\t$nopersampleinfo\n";
  print $statshandle "reference_clean\t$statscounters->{'refclean'}\t$nopersampleinfo\n";
  print $statshandle "dups_count\t$statscounters->{'duppositions'}\t$nopersampleinfo\n";
  print $statshandle "dups_portion\t" . sprintf( "%.2f", ( $statscounters->{'duppositions'} * 100 / $statscounters->{'reflength'} ) ) . "%\t$nopersampleinfo\n";
  print $statshandle "consensus_count\t$statscounters->{'consensus'}\t$nopersampleinfo\n";
  print $statshandle "consensus_portion\t" . sprintf( "%.2f", ( $statscounters->{'consensus'} * 100 / $statscounters->{'reflength'} ) ) . "%\t$nopersampleinfo\n";
  print $statshandle "num_total_variant\t$statscounters->{'totalvariant'}\t$nopersampleinfo\n";
  print $statshandle "num_best_snps\t$statscounters->{'bestsnps'}\t$nopersampleinfo\n";
  print $statshandle "best_position_breadth\t" . sprintf( "%.2f", ( $statscounters->{'bestpos'} * 100 / $statscounters->{'reflength'} ) ) . "%\t$nopersampleinfo\n";
  print $statshandle "core_genome_portion_snp\t" . ( ( $statscounters->{'allbreadth'} > 0 ) ? ( sprintf( "%.2f", ( $statscounters->{'bestsnps'} * 100 / $statscounters->{'allbreadth'} ) ) . "%\t" ) : '-\t' ) . "$nopersampleinfo\n";
  print $statshandle "best_position_portion_snp\t" . ( ( $statscounters->{'bestpos'} > 0 ) ? ( sprintf( "%.2f", ( $statscounters->{'bestsnps'} * 100 / $statscounters->{'bestpos'} ) ) . "%\t" ) : '-\t' ) . "\t$nopersampleinfo\n";
  print $statshandle "num_called\t$statscounters->{'allcalled'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "$statscounters->{'samplecalled'}{$samplefile}{$samplecolumn}\t"; } }
  print $statshandle "\n";
  print $statshandle "portion_called\t" . sprintf( "%.2f", ( $statscounters->{'allcalled'} * 100 / $statscounters->{'reflength'} ) ) . "%\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . sprintf( "%.2f", ( $statscounters->{'samplecalled'}{$samplefile}{$samplecolumn} * 100 / $statscounters->{'reflength'} ) ) . "%\t"; } }
  print $statshandle "\n";
  print $statshandle "num_passed_coverage_filter\t$statscounters->{'allcoveragepass'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . ( defined( $vcffiles->{$samplefile} ) ? "$statscounters->{'samplecoveragepass'}{$samplefile}{$samplecolumn}\t" : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "portion_passed_coverage_filter\t" . sprintf( "%.2f", ( $statscounters->{'allcoveragepass'} * 100 / $statscounters->{'reflength'} ) ) . "%\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . ( defined( $vcffiles->{$samplefile} ) ? ( sprintf( "%.2f", ( $statscounters->{'samplecoveragepass'}{$samplefile}{$samplecolumn} * 100 / $statscounters->{'reflength'} ) ) . "%\t" ) : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "num_passed_proportion_filter\t$statscounters->{'allproportionpass'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . ( defined( $vcffiles->{$samplefile} ) ? "$statscounters->{'sampleproportionpass'}{$samplefile}{$samplecolumn}\t" : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "portion_passed_proportion_filter\t" . sprintf( "%.2f", ( $statscounters->{'allproportionpass'} * 100 / $statscounters->{'reflength'} ) ) . "%\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . ( defined( $vcffiles->{$samplefile} ) ? ( sprintf( "%.2f", ( $statscounters->{'sampleproportionpass'}{$samplefile}{$samplecolumn} * 100 / $statscounters->{'reflength'} ) ) . "%\t" ) : "-\t" ); } }
  print $statshandle "\n";
  print $statshandle "num_NXdegen\t$statscounters->{'anynxdegen'}\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "$statscounters->{'samplenxdegen'}{$samplefile}{$samplecolumn}\t"; } }
  print $statshandle "\n";
  print $statshandle "portion_NXdegen\t" . sprintf( "%.2f", ( $statscounters->{'anynxdegen'} * 100 / $statscounters->{'reflength'} ) ) . "%\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . sprintf( "%.2f", ( $statscounters->{'samplenxdegen'}{$samplefile}{$samplecolumn} * 100 / $statscounters->{'reflength'} ) ) . "%\t"; } }
  print $statshandle "\n";
  my $avgsum = 0;
  my $avgcount = 0;
  my $avgstring = '';
  foreach my $samplefile ( sort keys( $samplecallsdata ) )
  {
    foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) )
    {
      if( defined( $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn} ) )
      {
        $avgcount++;
        $avgsum += $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn};
        $avgstring .= sprintf( "%.2f", ( $vcfcallcounts->{'depthsum'}{$samplefile}{$samplecolumn} / $statscounters->{'reflength'} ) ) . "\t";
      } else { $avgstring .= "-\t"; }
    }
  }
  print $statshandle "coverage_depth\t" . sprintf( "%.2f", ( ( $avgsum / $avgcount ) / $statscounters->{'reflength'} ) ) . "\t$avgstring\n";
  print $statshandle "coverage_breadth\t" . sprintf( "%.2f", ( $statscounters->{'allbreadth'} * 100 / $statscounters->{'reflength'} ) ) . "%\t";
  foreach my $samplefile ( sort keys( $samplecallsdata ) ){ foreach my $samplecolumn ( sort keys( $samplecallsdata->{$samplefile} ) ){ print $statshandle "" . sprintf( "%.2f", ( $statscounters->{'samplebreadth'}{$samplefile}{$samplecolumn} * 100 / $statscounters->{'reflength'} ) ) . "%\t"; } }
  print $statshandle "\n";
  print $statshandle "# reference_length: total number of positions found across all chromosomes in the reference\n";
  print $statshandle "# reference_clean: number of positions called A/C/G/T in the reference\n";
  print $statshandle "# dups_count: number of positions that appear to be in duplicated regions of the reference\n";
  print $statshandle "# dups_portion: dups_count / reference_length * 100%\n";
  print $statshandle "# consensus_count: number of positions in which all analyses of any one sample were in agreement\n";
  print $statshandle "# consensus_portion: consensus_count / reference_length * 100%\n";
  print $statshandle "# num_total_variant: number of positions where the reference and at least one analysis of at least one sample differed, were both called A/C/G/T, passed all filters, and were outside known-duplicated regions\n";
  print $statshandle "# num_best_snps: number of positions where the requirements for num_total_variant were met, there were no indels, all analyses of all samples were called A/C/G/T and passed all filters, and all analyses of any one sample were in agreement\n";
  print $statshandle "# best_position_breadth: number of positions where the requirements for num_best_snps would have been met if a snp were present and indels absent, divided by reference_length * 100%\n";
  print $statshandle "# core_genome_portion_snp: num_best_snps divided by the number of positions meething the requirements for total coverage_breadth\n";
  print $statshandle "# best_position_portion_snp: num_best_snps divided by the number of positions meeting the requirement for best_position_breadth\n";
  print $statshandle "# num_called: number of positions where the SNP caller made any call except 'N', regardless of filters\n";
  print $statshandle "# portion_called: num_called / reference_length * 100%\n";
  print $statshandle "# num_passed_coverage_filter: number of positions that met the depth filter requirements\n";
  print $statshandle "# portion_passed_coverage_filter: num_passed_coverage_filter / reference_length * 100%\n";
  print $statshandle "# num_passed_proportion_filter: number of positions that met the proportion filter requirements\n";
  print $statshandle "# portion_passed_proportion_filter: num_passed_proportion_filter / reference_length * 100%\n";
  print $statshandle "# num_NXdegen: number of positions called N, X, a degeneracy, or not called, regardless of filters\n";
  print $statshandle "# portion_NXdegen: num_NXdegen / reference_length * 100%\n";
  print $statshandle "# coverage_depth: the sum of the depths at each position, whether filtered or not, divided by reference_length\n";
  print $statshandle "# coverage_breadth: number of positions that had any call made and passed depth filtering, divided by reference_length * 100%\n";
} else { print STDERR "Could not open '$statisticsfile'!\n"; }


