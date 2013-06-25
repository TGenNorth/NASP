#!/usr/bin/perl -w

use strict;
use warnings;
use Vcf;
#use Data::Dumper;

if( @ARGV < 7 )
{
  print <<EOF;
Meant to be called from the distance pipeline automatically.

Usage:
matrix_distance.pl <min_coverage> <min_proportion> <mastermatrix.tsv> <samples.vcf> [samples2.vcf ...] <outputsamplematrix.tsv> <outputdistancefile.tsv> <filterlog.txt>
EOF
  exit();
}

my $mincoverage = shift();
my $minproportion = shift();
my $mastermatrixfile = shift();
my @samplecallsfiles = ();
while( my $samplecallsfile = shift() ){ push( @samplecallsfiles, $samplecallsfile ); } 
my $filterlogfile = pop( @samplecallsfiles );
my $outputdistancefile = pop( @samplecallsfiles );
my $outputmatrixfile = pop( @samplecallsfiles );

# Read in master SNP matrix
if( open( my $masterfilehandle, '<', $mastermatrixfile ) )
{
  my $masternickname = ( ( $mastermatrixfile =~ /([^\/]+?)(?:\.txt|\.xls|\.tsv)?$/ ) ? $1 : $mastermatrixfile );
  my $mastercolumnposition = -1;
  my $mastercolumnchrom = -1;
  my $mastercolumnstrand = -1;
  my $linefromfile = <$masterfilehandle>;
  chomp( $linefromfile );
  my @mastercolumns = split( /\t/, $linefromfile );
  my $mastermatrixdata = {};
  my $mastercolumnshavemeta = {};
  my $masterchromosomes = {};
  my $masterpositions = {};
  my $samplecallsdata = {};
  my @samplecolumns = ();
  my $i = 0;
  while( $i < scalar( @mastercolumns ) )
  {
    if( $mastercolumns[$i] =~ /^position$/i ){ $mastercolumnposition = $i; }
    if( $mastercolumns[$i] =~ /^chrom(?:osome)?$|^contig$/i ){ $mastercolumnchrom = $i; }
    if( $mastercolumns[$i] =~ /^strand$/i ){ $mastercolumnstrand = $i; }
    $i++;
  }
  if( $mastercolumnposition >= 0 )
  {
    while( $linefromfile = <$masterfilehandle> )
    {
      chomp( $linefromfile );
      my @columnsfromfile = split( /\t/, $linefromfile );
      my $currentchromosome = ( ( $mastercolumnchrom == -1 ) ? 'X' : $columnsfromfile[$mastercolumnchrom] );
      $masterchromosomes->{$currentchromosome} = 1;
      if( !defined( $masterpositions->{$currentchromosome} ) ){ $masterpositions->{$currentchromosome} = {}; }
      $masterpositions->{$currentchromosome}{$columnsfromfile[$mastercolumnposition]} = 1;
      $i = 0;
      while( $i < scalar( @mastercolumns ) )
      {
        if( !defined( $columnsfromfile[$i] ) || $columnsfromfile[$i] !~ /^[A-Za-z]$/ )
        {
          $mastercolumnshavemeta->{$mastercolumns[$i]} = 1;
          if( defined( $mastermatrixdata->{$mastercolumns[$i]} ) ){ delete( $mastermatrixdata->{$mastercolumns[$i]} ); }
        } elsif( !defined( $mastercolumnshavemeta->{$mastercolumns[$i]} ) )
        {
          if( !defined( $mastermatrixdata->{$mastercolumns[$i]} ) ){ $mastermatrixdata->{$mastercolumns[$i]} = {}; }
          if( !defined( $mastermatrixdata->{$mastercolumns[$i]}{$currentchromosome} ) ){ $mastermatrixdata->{$mastercolumns[$i]}{$currentchromosome} = {}; }
          $mastermatrixdata->{$mastercolumns[$i]}{$currentchromosome}{$columnsfromfile[$mastercolumnposition]} = $columnsfromfile[$i];
          if( ( $mastercolumnstrand > -1 ) && ( $columnsfromfile[$mastercolumnstrand] eq "-" ) )
          {
            $mastermatrixdata->{$mastercolumns[$i]}{$currentchromosome}{$columnsfromfile[$mastercolumnposition]} =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
          }
        }
        $i++;
      }
    }

    # Read in called SNPs
    my $referencecalls = {};
    open( my $filterlogfilehandle, '>', $filterlogfile );
    foreach my $samplecallsfile ( @samplecallsfiles )
    {
      my $samplenickname = ( ( $samplecallsfile =~ /([^\/]+?)(?:\.vcf)?$/ ) ? $1 : $samplecallsfile );
      if( -e( $samplecallsfile ) && ( my $samplefilehandle = eval { Vcf->new( 'file'=>$samplecallsfile ); } ) )
      {
        $samplefilehandle->parse_header();
        my (@samplelist) = $samplefilehandle->get_samples();
        foreach my $currentsample (@samplelist)
        {
          $samplecallsdata->{"${samplenickname}:${currentsample}"} = {};
          push( @samplecolumns, "${samplenickname}:${currentsample}" );
        }
        my $positiondata = {};
        while( $positiondata = $samplefilehandle->next_data_hash() )
        {
          #print Dumper( $positiondata );
          my $currentchromosome = ( ( $mastercolumnchrom == -1 ) ? 'X' : $positiondata->{'CHROM'} );
          if( defined( $masterpositions->{$currentchromosome}{$positiondata->{'POS'}} ) )
          {
            if( !defined( $referencecalls->{$currentchromosome} ) ){ $referencecalls->{$currentchromosome} = {}; }
            if( !defined( $referencecalls->{$currentchromosome}{$positiondata->{'POS'}} ) ){ $referencecalls->{$currentchromosome}{$positiondata->{'POS'}} = $positiondata->{'REF'}; }
            foreach my $currentsample (@samplelist)
            {
              if( !defined( $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome} ) ){ $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome} = {}; }
              my $calledallele;
              if( ( $positiondata->{'ALT'}[0] eq '.' ) || ( $positiondata->{'gtypes'}{$currentsample}{'GT'} == 0 ) ){ $calledallele = $positiondata->{'REF'}; }
              else { $calledallele = $positiondata->{'ALT'}[ $positiondata->{'gtypes'}{$currentsample}{'GT'} - 1 ]; }
              my $positiondepth = ( defined( $positiondata->{'INFO'}{'DP'} ) ? ( $positiondata->{'INFO'}{'DP'} / scalar( @samplelist ) ) : 0 );
              if( defined( $positiondata->{'gtypes'}{$currentsample}{'DP'} ) ){ $positiondepth = $positiondata->{'gtypes'}{$currentsample}{'DP'}; }
              my $callfrequency = -1;
              if( defined( $positiondata->{'gtypes'}{$currentsample}{'AD'} ) ){ $callfrequency = ( split( /,/, $positiondata->{'gtypes'}{$currentsample}{'AD'} ) )[$positiondata->{'gtypes'}{$currentsample}{'GT'}]; }
              if( defined( $positiondata->{'INFO'}{'PL'} ) )
              {
                my $callfrequencypl = scalar( grep( uc( $_ ) eq uc( $calledallele ), split( '', $positiondata->{'INFO'}{'PL'} ) ) );
                if( ( $callfrequency == -1 ) || ( $callfrequencypl < $callfrequency ) ){ $callfrequency = $callfrequencypl; }
              }
              if( $positiondepth < $mincoverage )
              {
                $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome}{$positiondata->{'POS'}} = 'X';
                print $filterlogfilehandle "Call changed to 'X' for low coverage ($positiondepth < $mincoverage) on sample '$currentsample' from file '$samplenickname' at position $positiondata->{'POS'} on chromosome '$currentchromosome'.\n";
              } elsif( ( $callfrequency != -1 ) && ( $callfrequency < $mincoverage ) )
              {
                $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome}{$positiondata->{'POS'}} = 'X';
                print $filterlogfilehandle "Call changed to 'X' for low coverage ($callfrequency < $mincoverage) on sample '$currentsample' from file '$samplenickname' at position $positiondata->{'POS'} on chromosome '$currentchromosome'.\n";
              } elsif( ( $callfrequency != -1 ) && ( ( $callfrequency / $positiondepth ) < $minproportion ) )
              {
                $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome}{$positiondata->{'POS'}} = 'N';
                print $filterlogfilehandle "Call changed to 'N' for insufficient proportion ($callfrequency/$positiondepth < $minproportion) on sample '$currentsample' from file '$samplenickname' at position $positiondata->{'POS'} on chromosome '$currentchromosome'.\n";
              } else { $samplecallsdata->{"${samplenickname}:${currentsample}"}{$currentchromosome}{$positiondata->{'POS'}} = $calledallele; }
            }
          }
        }
        $samplefilehandle->close();
      } else { print STDERR "Could not open '$samplecallsfile'!\n"; }
    }
    close( $filterlogfilehandle );

    # Make SNP matrix
    if( open( my $matrixfilehandle, '>', $outputmatrixfile ) )
    {
      my $patternarray = {};
      my $nextpatternnum = 1;
      print $matrixfilehandle "Chromosome\tPosition\tReference\t";
      foreach my $samplecolumn (@samplecolumns){ print $matrixfilehandle "$samplecolumn\t"; }
      print $matrixfilehandle "#Ref\t#SNP\t#NXindel\t%A\t%C\t%G\t%T\t%NXindel\tpattern\tpat_num\n";
      foreach my $currentchromosome ( sort keys( %$masterpositions ) )
      {
        foreach my $currentposition ( sort {$a <=> $b} keys( $masterpositions->{$currentchromosome} ) )
        {
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
          print $matrixfilehandle "$currentchromosome\t$currentposition\t$referencecalls->{$currentchromosome}{$currentposition}\t";
          foreach my $samplecolumn (@samplecolumns)
          {
            if( !defined( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} ) ){ $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} = "N"; }
            print $matrixfilehandle "$samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition}\t";
            if( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[ACGTacgt]$/ ){ $allelecounts->{uc($samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition})}++; }
            elsif( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[Uu]$/ ){ $allelecounts->{"T"}++; }
            else { $allelecounts->{"N"}++; }
            if( ( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[ACGTUacgtu]$/ ) && ( $referencecalls->{$currentchromosome}{$currentposition} =~ /^[ACGTUacgtu]$/ ) )
            {
              if( ( uc( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} ) eq uc( $referencecalls->{$currentchromosome}{$currentposition} ) ) || ( ( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[TUtu]$/ ) && ( $referencecalls->{$currentchromosome}{$currentposition} =~ /^[TUtu]$/ ) ) )
              {
                $allelecounts->{"ref"}++;
              } else { $allelecounts->{"snp"}++; }
            }
            $allelecounts->{"sum"}++;
            if( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition} =~ /^[NnXx]$/ )
            {
              $callpattern .= "N";
            } else
            {
              if( !defined( $patternassignments->{$samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition}} ) )
              {
                $patternassignments->{$samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition}} = $nextassignmentnum;
                if( $nextassignmentnum ne '+' ){ $nextassignmentnum = ( ( $nextassignmentnum < 9 ) ? ( $nextassignmentnum + 1 ) : '+' ); }
              }
              $callpattern .= $patternassignments->{$samplecallsdata->{$samplecolumn}{$currentchromosome}{$currentposition}};
            }
          }
          print $matrixfilehandle "$allelecounts->{'ref'}\t";
          print $matrixfilehandle "$allelecounts->{'snp'}\t";
          print $matrixfilehandle "$allelecounts->{'N'}\t";
          print $matrixfilehandle sprintf( "%.2f", ( $allelecounts->{"A"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          print $matrixfilehandle sprintf( "%.2f", ( $allelecounts->{"C"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          print $matrixfilehandle sprintf( "%.2f", ( $allelecounts->{"G"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          print $matrixfilehandle sprintf( "%.2f", ( $allelecounts->{"T"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          print $matrixfilehandle sprintf( "%.2f", ( $allelecounts->{"N"} / $allelecounts->{"sum"} * 100 ) ) . "\t";
          print $matrixfilehandle "'$callpattern'\t";
          if( !defined( $patternarray->{$callpattern} ) ){ $patternarray->{$callpattern} = $nextpatternnum++; }
          print $matrixfilehandle "$patternarray->{$callpattern}\n";
        }
      }
    } else { print STDERR "Could not open '$outputmatrixfile'!\n"; }

    # Make distance matrix
    if( open( my $distancefilehandle, '>', $outputdistancefile ) )
    {
      print $distancefilehandle "\t";
      foreach my $samplecolumn (@samplecolumns){ print $distancefilehandle "$samplecolumn - distance\t$samplecolumn - comparisons\t$samplecolumn - mismatches\t"; }
      print $distancefilehandle "\n";
      foreach my $mastercolumn (@mastercolumns)
      {
        if( !defined( $mastercolumnshavemeta->{$mastercolumn} ) )
        {
          print $distancefilehandle "${masternickname}:${mastercolumn}\t";
          foreach my $samplecolumn (@samplecolumns)
          {
            my $numdifferent = 0;
            my $numcompared = 0;
            foreach my $currentchromosome ( keys( %$masterpositions ) )
            {
              foreach my $masterposition ( keys( $masterpositions->{$currentchromosome} ) )
              {
                if( defined( $mastermatrixdata->{$mastercolumn}{$currentchromosome}{$masterposition} ) && defined( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} ) && ( $mastermatrixdata->{$mastercolumn}{$currentchromosome}{$masterposition} =~ /^[ACGTUacgtu]$/ ) && ( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} =~ /^[ACGTUacgtu]$/ ) )
                {
                  $numcompared++;
                  if( uc( $mastermatrixdata->{$mastercolumn}{$currentchromosome}{$masterposition} ) ne uc( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} ) ){ $numdifferent++; }
                }
              }
            }
            print $distancefilehandle '' . ( ( $numcompared > 0 ) ? $numdifferent/$numcompared : "X" ) . "\t$numcompared\t$numdifferent\t";
          }
          print $distancefilehandle "\n";
        }
      }
      foreach my $comparesamplecolumn (@samplecolumns)
      {
        print $distancefilehandle "${comparesamplecolumn}\t";
        foreach my $samplecolumn (@samplecolumns)
        {
          my $numdifferent = 0;
          my $numcompared = 0;
          foreach my $currentchromosome ( keys( %$masterpositions ) )
          {
            foreach my $masterposition ( keys( $masterpositions->{$currentchromosome} ) )
            {
              if( defined( $samplecallsdata->{$comparesamplecolumn}{$currentchromosome}{$masterposition} ) && defined( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} ) && ( $samplecallsdata->{$comparesamplecolumn}{$currentchromosome}{$masterposition} =~ /^[ACGTUacgtu]$/ ) && ( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} =~ /^[ACGTUacgtu]$/ ) )
              {
                $numcompared++;
                if( uc( $samplecallsdata->{$comparesamplecolumn}{$currentchromosome}{$masterposition} ) ne uc( $samplecallsdata->{$samplecolumn}{$currentchromosome}{$masterposition} ) ){ $numdifferent++; }
              }
            }
          }
          print $distancefilehandle '' . ( ( $numcompared > 0 ) ? $numdifferent/$numcompared : "X" ) . "\t$numcompared\t$numdifferent\t";
        }
        print $distancefilehandle "\n";
      }
      close( $distancefilehandle );
    } else { print STDERR "Could not open '$outputdistancefile'!\n"; }

  } else { print STDERR "Could not find 'position' column in '$mastermatrixfile'!\n"; }
  close( $masterfilehandle );
} else { print STDERR "Could not open '$mastermatrixfile'!\n"; }

# SolSNP reference
#          'FORMAT' => [
#                        'GT',
#                        'GQ'
#                      ],
#          'QUAL' => '30.0',
#          'ID' => '.',
#          'CHROM' => 'Ba-Ames_Ancestor',
#          'INFO' => {
#                      'DP' => '50',
#                      'CL' => 'CC',
#                      'PL' => 'cccccccccccccccccccccccccccccccccccccccccccccccccc',
#                      'AR' => '0.00'
#                    },
#          'FILTER' => [
#                        'PASS'
#                      ],
#          'gtypes' => {
#                        'Ba-4599_4.bam' => {
#                                             'GQ' => '30.0',
#                                             'GT' => '0'
#                                           }
#                      },
#          'REF' => 'C',
#          'ALT' => [
#                     '.'
#                   ],
#          'POS' => '14287'

# SolSNP SNP
#          'FORMAT' => [
#                        'GT',
#                        'GQ'
#                      ],
#          'QUAL' => '30.0',
#          'ID' => '.',
#          'CHROM' => 'Ba-Ames_Ancestor',
#          'INFO' => {
#                      'DP' => '60',
#                      'CL' => 'TT',
#                      'PL' => 'tTtatTTttttTtTtTttttTtTTTttTtTtTTttttTtTTTtttTTtTtTTTTttttTT',
#                      'AR' => '1.00'
#                    },
#          'FILTER' => [
#                        'PASS'
#                      ],
#          'gtypes' => {
#                        'Banthracis-A0362_1_noINDEL_baq.bam' => {
#                                                                  'GQ' => '30.0',
#                                                                  'GT' => '1'
#                                                                }
#                      },
#          'REF' => 'C',
#          'ALT' => [
#                     'T'
#                   ],
#          'POS' => '3306423'

# GATK reference
#          'FORMAT' => [
#                        'GT'
#                      ],
#          'QUAL' => '3.01',
#          'ID' => '.',
#          'CHROM' => 'Ba-Ames_Ancestor',
#          'INFO' => {
#                      'MQ' => '57.20',
#                      'DP' => '143',
#                      'MQ0' => '0'
#                    },
#          'FILTER' => [
#                        'LowQual'
#                      ],
#          'gtypes' => {
#                        'Ba-4599_4' => {
#                                         'GT' => '.'
#                                       }
#                      },
#          'REF' => 'G',
#          'ALT' => [
#                     '.'
#                   ],
#          'POS' => '4489'

# GATK SNP
#          'FORMAT' => [
#                        'GT',
#                        'AD',
#                        'DP',
#                        'GQ',
#                        'MLPSAC',
#                        'MLPSAF',
#                        'PL'
#                      ],
#          'QUAL' => '1704',
#          'ID' => '.',
#          'CHROM' => 'Ba-Ames_Ancestor',
#          'INFO' => {
#                      'AC' => '1',
#                      'HaplotypeScore' => '0.8667',
#                      'FS' => '0.000',
#                      'DP' => '66',
#                      'MQ0' => '0',
#                      'QD' => '25.82',
#                      'Dels' => '0.00',
#                      'AF' => '1.00',
#                      'MQ' => '55.73',
#                      'AN' => '1',
#                      'MLEAF' => '1.00',
#                      'SB' => '-7.860e+02',
#                      'MLEAC' => '1'
#                    },
#          'FILTER' => [
#                        '.'
#                      ],
#          'gtypes' => {
#                        'Ba-4599_4' => {
#                                         'GQ' => '99',
#                                         'DP' => '66',
#                                         'MLPSAF' => '1.00',
#                                         'PL' => '1734,0',
#                                         'MLPSAC' => '1',
#                                         'AD' => '0,66',
#                                         'GT' => '1'
#                                       }
#                      },
#          'REF' => 'T',
#          'ALT' => [
#                     'C'
#                   ],
#          'POS' => '2376'


