#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;

# Some constants for tweaking.
# Some of these should be set dynamically... someday.
my $naspversion = "0.8.4";
my $finddupspath = "find_duplicates.pl";
my $gigsofmemforindex = "2";
my $wallhoursforindex = "1";
my $convertexternalpath = "convert_external_genome.pl";
my $gigsofmemforexternal = "2";
my $wallhoursforexternal = "1";
my $bwapath = "bwa";
my $gigsofmemforbwa = "12";
my $numcpusforbwa = "4";
my $wallhoursforbwa = "36";
my $defaultbwaalnargs = "";
my $defaultbwasampeargs = "";
my $defaultbwamemargs = "";
my $novopath = "novoalign";
my $gigsofmemfornovo = "12";
my $numcpusfornovo = "4";
my $wallhoursfornovo = "36";
my $defaultnovoargs = "-r all";
my $novopairedargs = "-i PE 500,100";
my $novoindexpath = "novoindex";
my $gatkpath = "GenomeAnalysisTK.jar";
my $gigsofmemforgatk = "12";
my $numcpusforgatk = "4";
my $wallhoursforgatk = "36";
my $defaultgatkargs = "-stand_call_conf 100 -stand_emit_conf 100";
my $dictgeneratorpath = "CreateSequenceDictionary.jar";
my $samtoolspath = "samtools";
my $solsnppath = "SolSNP.jar";
my $gigsofmemforsolsnp = "12";
my $numcpusforsolsnp = "4";
my $wallhoursforsolsnp = "36";
my $defaultsolsnpargs = "";
my $varscanpath = "VarScan.jar";
my $gigsofmemforvarscan = "2";
my $numcpusforvarscan = "1";
my $wallhoursforvarscan = "24";
my $defaultvarscanargs = "";
my $bcftoolspath = "bcftools";
my $gigsofmemforbcftools = "2";
my $numcpusforbcftools = "1";
my $wallhoursforbcftools = "24";
my $gigsofmemtomakematrix = "4";
my $wallhourstomakematrix = "48";
my $matrixmakingscript = "vcf_to_matrix.pl";
my $gigsofmemfordistancecalc = "60";
my $wallhoursfordistancecalc = "36";
my $distancecalcscript = "matrix_distance.pl";

# Standalone executables will be checked for in the
# $PATH by default, using the normal methods.
# Jar files will be searched for in the $PATH plus:
my @jarpaths =
(
  "/usr/share/java",
  "/media/lumberyard/bin",
  "/tnorth/bin",
  "~/bin",
  "~/jars",
  "~/tools",
  "."
);

if( @ARGV < 1 || @ARGV > 3 )
{
  print <<EOF;
This is the experimental "Northern Arizona SNP Pipeline", version $naspversion.

Usage:
nasp <reference.fasta> [read_folder [output_folder]]
EOF
  exit();
}

# Find jar files
push( @jarpaths, split( ':', $ENV{"PATH"} ) );
foreach my $jarpath (@jarpaths)
{
  $jarpath =~ s/^\~/$ENV{"HOME"}/;
  if( ( $gatkpath !~ /\// ) && -e( "$jarpath/$gatkpath" ) ){ $gatkpath = "$jarpath/$gatkpath"; }
  if( ( $dictgeneratorpath !~ /\// ) && -e( "$jarpath/$dictgeneratorpath" ) ){ $dictgeneratorpath = "$jarpath/$dictgeneratorpath"; }
  if( ( $solsnppath !~ /\// ) && -e( "$jarpath/$solsnppath" ) ){ $solsnppath = "$jarpath/$solsnppath"; }
  if( ( $varscanpath !~ /\// ) && -e( "$jarpath/$varscanpath" ) ){ $varscanpath = "$jarpath/$varscanpath"; }
}

# Process command line arguments and call subs.
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

# The calling program parses command-line arguments and then calls this function.
# Flow should be changed in the future to allow for more types of user interface.
sub nasp
{
  my $referencefastafile = shift();
  my $readfilefolder = shift();
  my $outputfilefolder = shift();
  my $mastermatrixfile = shift();
  print "Welcome to nasp version $naspversion.\n";
  print "* Starred features might be broken.\n";
  
  # This section is the interactive command-line user input section.
  # A web form or Java interface would replace this section.
  if( -e( $outputfilefolder ) )
  {
    print "\nOutput folder '$outputfilefolder' already exists!\nFiles in it may be overwritten!\nShould we continue anyway [N]? ";
    if( <> !~ /^[yY]/ ){ die( "Operation cancelled!\n" ); }
  }

  # Reference/external section
  my $finddupregions = 0;
  my $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "\nDo you want to check the reference for duplicated regions\nand skip SNPs that fall in those regions [Y]? ";
    $userinput = <>;
    if( $userinput !~ /^[Nn]/ ){ $finddupregions = 1; }
  }
  my $findexternalfastas = 0;
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "\nDo you have fasta files for external genomes you wish to include [Y]? ";
    $userinput = <>;
    if( $userinput !~ /^[Nn]/ ){ $findexternalfastas = 1; }
  }

  # Filters section
  $userinput = '';
  while( $userinput !~ /^\d{1,31}$/ )
  {
    print "\nThis pipeline can do filtering based on coverage.\nIf you do not want filtering based on coverage, enter 0.\nWhat is your minimum coverage threshold [10]? ";
    $userinput = <>;
    chomp( $userinput );
    if( $userinput =~ /^$/ ){ $userinput = 10; }
  }
  my $mincoverage = $userinput;
  $userinput = '';
  while( $userinput !~ /^(?:[01](?:\.0{0,31})?|0?\.\d{1,31})$/ )
  {
    print "\nThis pipeline can do filtering based on the proportion of reads that match the call made by the SNP caller.\nIf you do not want filtering based on proportion, enter 0.\nWhat is the minimum acceptable proportion [0.9]? ";
    $userinput = <>;
    chomp( $userinput );
    if( $userinput =~ /^$/ ){ $userinput = .9; }
  }
  my $minproportion = $userinput;
  
  # Alignment section
  my $runbwa = 0;
  my $runbwamem = 0;
  my $runnovo = 0;
  my $findbams = 0;
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "\nThis pipeline currently supports two aligners, BWA and Novoalign, and you can provide pre-aligned BAM files.\nYou can choose as many options as you want.\nWould you like to run BWA samp/se [N]?* ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ ){ $runbwa = 1; }
  }
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "Would you like to run BWA mem [Y]? ";
    $userinput = <>;
    if( $userinput !~ /^[Nn]/ ){ $runbwamem = 1; }
  }
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "Would you like to run Novoalign [Y]? ";
    $userinput = <>;
    if( $userinput !~ /^[Nn]/ ){ $runnovo = 1; }
  }
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "Would you like to provide pre-aligned SAM/BAM files [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ ){ $findbams = 1; }
  }
  my $findpaired = 0;
  if( $runbwa || $runbwamem || $runnovo )
  {
    $userinput = 'X';
    while( $userinput !~ /^$|^[YNyn]/ )
    {
      print "Are your read files paired [Y]? ";
      $userinput = <>;
      if( $userinput !~ /^[Nn]/ ){ $findpaired = 1; }
    }
    if( $runbwa || $runbwamem )
    {
      print "Would you like to set advanced BWA settings [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  Would you like to use an alternate BWA version [N]? ";
        $userinput = <>;
        if( $userinput =~ /^[Yy]/ )
        {
          print "  What is the path to the BWA runtime you wish to use [system default]? ";
          $userinput = <>;
          chomp( $userinput );
          if( length( $userinput ) ){ $bwapath = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,4}$/ )
        {
          print "  How much memory will BWA require [$gigsofmemforbwa]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $gigsofmemforbwa = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,2}$/ )
        {
          print "  How many CPUs do you want BWA to use [$numcpusforbwa]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $numcpusforbwa = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,3}$/ )
        {
          print "  How many hours will BWA take to run [$wallhoursforbwa]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $wallhoursforbwa = $userinput; }
        }
        if( $runbwa )
        {
          $userinput = '';
          print "  What additional arguments would you like to pass to 'bwa aln' [ $defaultbwaalnargs ]? ";
          $userinput = <>;
          chomp( $userinput );
          # This could use some sanity-checking, and perhaps multiple smarter questions.
          if( length( $userinput ) ){ $defaultbwaalnargs = $userinput; }
          $userinput = '';
          print "  What additional arguments would you like to pass to 'bwa samp/se' [ $defaultbwasampeargs ]? ";
          $userinput = <>;
          chomp( $userinput );
          # This could use some sanity-checking, and perhaps multiple smarter questions.
          if( length( $userinput ) ){ $defaultbwasampeargs = $userinput; }
        }
        if( $runbwamem )
        {
          $userinput = '';
          print "  What additional arguments would you like to pass to 'bwa mem' [ $defaultbwamemargs ]? ";
          $userinput = <>;
          chomp( $userinput );
          # This could use some sanity-checking, and perhaps multiple smarter questions.
          if( length( $userinput ) ){ $defaultbwamemargs = $userinput; }
        }
      }
    }
    if( $runnovo )
    {
      print "Would you like to set advanced Novoalign settings [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  Would you like to use an alternate Novoalign version [N]? ";
        $userinput = <>;
        if( $userinput =~ /^[Yy]/ )
        {
          print "  What is the path to the Novoalign runtime you wish to use [system default]? ";
          $userinput = <>;
          chomp( $userinput );
          if( length( $userinput ) ){ $novopath = $userinput; }
          print "  What is the path to the Novoindex runtime you wish to use [system default]? ";
          $userinput = <>;
          chomp( $userinput );
          if( length( $userinput ) ){ $novoindexpath = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,4}$/ )
        {
          print "  How much memory will Novoalign require [$gigsofmemfornovo]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $gigsofmemfornovo = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,2}$/ )
        {
          print "  How many CPUs do you want Novoalign to use [$numcpusfornovo]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $numcpusfornovo = $userinput; }
        }
        $userinput = 'X';
        while( $userinput !~ /^\d{0,3}$/ )
        {
          print "  How many hours will Novoalign take to run [$wallhoursfornovo]? ";
          $userinput = <>;
          chomp( $userinput );
          if( $userinput !~ /^$/ ){ $wallhoursfornovo = $userinput; }
        }
        if( $findpaired )
        {
          $userinput = '';
          print "  What paired-alignment arguments would you like to pass to Novoalign [ $novopairedargs ]? ";
          $userinput = <>;
          chomp( $userinput );
          # This could use some sanity-checking, and perhaps multiple smarter questions.
          if( length( $userinput ) ){ $novopairedargs = $userinput; }
        }
        $userinput = '';
        print "  What additional arguments would you like to pass to Novoalign [ $defaultnovoargs ]? ";
        $userinput = <>;
        chomp( $userinput );
        # This could use some sanity-checking, and perhaps multiple smarter questions.
        if( length( $userinput ) ){ $defaultnovoargs = $userinput; }
      }
    }
  }
  
  # SNP calling section
  my $rungatk = 0;
  my $runsolsnp = 0;
  my $runvarscan = 0;
  my $runsamtools = 0;
  my $findvcfs = 0;
  if( $runbwa || $runbwamem || $runnovo || $findbams )
  {
    $userinput = 'X';
    while( $userinput !~ /^$|^[YNyn]/ )
    {
      print "\nThis pipeline currently supports four SNP callers: GATK, SolSNP, VarScan, and SAMtools, and you can provide VCF files.\nYou can choose as many options as you want.\nWould you like to run GATK [Y]? ";
      $userinput = <>;
      if( $userinput !~ /^[Nn]/ ){ $rungatk = 1; }
    }
    $userinput = 'X';
    while( $userinput !~ /^$|^[YNyn]/ )
    {
      print "Would you like to run SolSNP [Y]? ";
      $userinput = <>;
      if( $userinput !~ /^[Nn]/ ){ $runsolsnp = 1; }
    }
    $userinput = 'X';
    while( $userinput !~ /^$|^[YNyn]/ )
    {
      print "Would you like to run VarScan [Y]? ";
      $userinput = <>;
      if( $userinput !~ /^[Nn]/ ){ $runvarscan = 1; }
    }
    $userinput = 'X';
    while( $userinput !~ /^$|^[YNyn]/ )
    {
      print "Would you like to run SAMtools [Y]? ";
      $userinput = <>;
      if( $userinput !~ /^[Nn]/ ){ $runsamtools = 1; }
    }
  }
  $userinput = 'X';
  while( $userinput !~ /^$|^[YNyn]/ )
  {
    print "Would you like to provide VCF files [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ ){ $findvcfs = 1; }
  }
  if( $rungatk )
  {
    print "Would you like to set advanced GATK settings [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ )
    {
      print "  Would you like to use an alternate GATK version [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  What is the path to the GATK runtime you wish to use [system default]? ";
        $userinput = <>;
        chomp( $userinput );
        if( length( $userinput ) ){ $gatkpath = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,4}$/ )
      {
        print "  How much memory will GATK require [$gigsofmemforgatk]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $gigsofmemforgatk = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,2}$/ )
      {
        print "  How many CPUs do you want GATK to use [$numcpusforgatk]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $numcpusforgatk = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,3}$/ )
      {
        print "  How many hours will GATK take to run [$wallhoursforgatk]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $wallhoursforgatk = $userinput; }
      }
      $userinput = 'X';
      print "  What additional arguments would you like to pass to GATK [ $defaultgatkargs ]? ";
      $userinput = <>;
      chomp( $userinput );
      # This could use some sanity-checking, and perhaps multiple smarter questions.
      if( length( $userinput ) ){ $defaultgatkargs = $userinput; }
    }
  }
  if( $runsolsnp )
  {
    print "Would you like to set advanced SolSNP settings [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ )
    {
      print "  Would you like to use an alternate SolSNP version [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  What is the path to the SolSNP runtime you wish to use [system default]? ";
        $userinput = <>;
        chomp( $userinput );
        if( length( $userinput ) ){ $solsnppath = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,4}$/ )
      {
        print "  How much memory will SolSNP require [$gigsofmemforsolsnp]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $gigsofmemforsolsnp = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,2}$/ )
      {
        print "  How many CPUs do you want SolSNP to use [$numcpusforsolsnp]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $numcpusforsolsnp = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,3}$/ )
      {
        print "  How many hours will SolSNP take to run [$wallhoursforsolsnp]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $wallhoursforsolsnp = $userinput; }
      }
      $userinput = '';
      print "  What additional arguments would you like to pass to SolSNP [ $defaultsolsnpargs ]? ";
      $userinput = <>;
      chomp( $userinput );
      # This could use some sanity-checking, and perhaps multiple smarter questions.
      if( length( $userinput ) ){ $defaultsolsnpargs = $userinput; }
    }
  }
  if( $runvarscan )
  {
    print "Would you like to set advanced VarScan settings [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ )
    {
      print "  Would you like to use an alternate VarScan version [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  What is the path to the VarScan runtime you wish to use [system default]? ";
        $userinput = <>;
        chomp( $userinput );
        if( length( $userinput ) ){ $varscanpath = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,4}$/ )
      {
        print "  How much memory will VarScan require [$gigsofmemforvarscan]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $gigsofmemforvarscan = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,2}$/ )
      {
        print "  How many CPUs do you want VarScan to use [$numcpusforvarscan]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $numcpusforvarscan = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,3}$/ )
      {
        print "  How many hours will VarScan take to run [$wallhoursforvarscan]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $wallhoursforvarscan = $userinput; }
      }
      $userinput = '';
      print "  What additional arguments would you like to pass to VarScan [ $defaultvarscanargs ]? ";
      $userinput = <>;
      chomp( $userinput );
      # This could use some sanity-checking, and perhaps multiple smarter questions.
      if( length( $userinput ) ){ $defaultvarscanargs = $userinput; }
    }
  }
  if( $runsamtools )
  {
    print "Would you like to set advanced SAMtools settings [N]? ";
    $userinput = <>;
    if( $userinput =~ /^[Yy]/ )
    {
      print "  Would you like to use an alternate SAMtools version [N]? ";
      $userinput = <>;
      if( $userinput =~ /^[Yy]/ )
      {
        print "  What is the path to the 'bcftools' runtime you wish to use [system default]? ";
        $userinput = <>;
        chomp( $userinput );
        if( length( $userinput ) ){ $bcftoolspath = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,4}$/ )
      {
        print "  How much memory will SAMtools require [$gigsofmemforbcftools]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $gigsofmemforbcftools = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,2}$/ )
      {
        print "  How many CPUs do you want SAMtools to use [$numcpusforbcftools]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $numcpusforbcftools = $userinput; }
      }
      $userinput = 'X';
      while( $userinput !~ /^\d{0,3}$/ )
      {
        print "  How many hours will SAMtools take to run [$wallhoursforbcftools]? ";
        $userinput = <>;
        chomp( $userinput );
        if( $userinput !~ /^$/ ){ $wallhoursforbcftools = $userinput; }
      }
    }
  }
  
  # Collect fasta/fastq/bam/vcf files for the pipeline.
  # We do it now so that we can give the user feedback about what we found.
  my @fastafilelist = ();
  my @fastqfilelist = ();
  my @bamfilelist = ();
  my @vcffilelist = ();
  if( opendir( my $readfolderhandle, $readfilefolder ) )
  {
    my $possiblereadfile = "";
    while( $possiblereadfile = readdir( $readfolderhandle ) )
    {
      # fasta
      if( $findexternalfastas && ( $possiblereadfile =~ /^(?:.*\/)?([^\/]+?)\.[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?$/ ) )
      {
        my $possiblefastanickname = $1;
        if( $referencefastafile =~ /^(?:.*\/)?([^\/]+?)(?:\.[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?)?$/ )
        {
          my $referencenickname = $1;
          if( $referencenickname ne $possiblefastanickname ){ push( @fastafilelist, $possiblereadfile ); }
        }
      }
      # fastq
      if( ( $runbwa || $runbwamem || $runnovo ) && ( $possiblereadfile =~ /\.f(?:ast)?q(?:\.gz)?$|_sequence\.txt(?:\.gz)?$/i ) )
      {
        my $matchfound = 0;
        if( $findpaired )
        {
          if( ( $possiblereadfile =~ /\.f(?:ast)?q(?:\.gz)?$/i ) && ( $possiblereadfile =~ /_[Rr][12]_/ ) )
          {
            my $nametomatch = $possiblereadfile;
            $nametomatch =~ s/_[Rr][12]_//g;
            my $i = 0;
            while( !( $matchfound ) && defined( $fastqfilelist[$i] ) )
            {
              if( scalar( @{$fastqfilelist[$i]} ) == 1 )
              {
                my $nametocheck = $fastqfilelist[$i]->[0];
                $nametocheck =~ s/_[Rr][12]_//g;
                if( $nametomatch eq $nametocheck )
                {
                  my @readpair = sort( ( $fastqfilelist[$i]->[0], $possiblereadfile ) );
                  $fastqfilelist[$i] = \@readpair;
                  $matchfound = 1;
                }
              }
              $i++;
            }
          } elsif( $possiblereadfile =~ /_[12]_sequence\.txt(?:\.gz)?$/i )
          {
            my $nametomatch = lc( $possiblereadfile );
            $nametomatch =~ s/_[12](_sequence\.txt(?:\.gz)?)$/$1/;
            my $i = 0;
            while( !( $matchfound ) && defined( $fastqfilelist[$i] ) )
            {
              if( scalar( @{$fastqfilelist[$i]} ) == 1 )
              {
                my $nametocheck = lc( $fastqfilelist[$i]->[0] );
                $nametocheck =~ s/_[12](_sequence\.txt(?:\.gz)?)$/$1/;
                if( $nametomatch eq $nametocheck )
                {
                  my @readpair = sort( ( $fastqfilelist[$i]->[0], $possiblereadfile ) );
                  $fastqfilelist[$i] = \@readpair;
                  $matchfound = 1;
                }
              }
              $i++;
            }
          }
        }
        if( !( $matchfound ) )
        {
          my @readpair = ( $possiblereadfile );
          push( @fastqfilelist, \@readpair );
        }
      }
      # bam
      if( $findbams && ( $possiblereadfile =~ /\.[SsBb][Aa][Mm]$/ ) ){ push( @bamfilelist, $possiblereadfile ); }
      # vcf
      if( $findvcfs && ( $possiblereadfile =~ /\.[Vv][Cc][Ff]$/ ) ){ push( @vcffilelist, $possiblereadfile ); }
    }
    my @finalfilelist = map { "vcf,pre-aligned,pre-called,::$readfilefolder/" . $_ } @vcffilelist;
    closedir( $readfolderhandle );
    if( ( scalar( @fastafilelist ) + scalar( @fastqfilelist ) + scalar( @bamfilelist ) + scalar( @vcffilelist ) ) >= 1 )
    {

      # Log the configuration used in this run, both to a file and to the screen
      if( !( -e( $outputfilefolder ) ) ){ mkdir( $outputfilefolder ); }
      open( my $loghandle, ">", "$outputfilefolder/configuration_log.txt" );
      print $loghandle "nasp version $naspversion begun at " . localtime() . ".\n";
      print "\nConfiguration:\n";
      print $loghandle "Folders:\n  Input folder: $readfilefolder\n  Output folder: $outputfilefolder\n";
      print "Folders:\n  Input folder: $readfilefolder\n  Output folder: $outputfilefolder\n";
      print $loghandle "Fasta data:\n  Reference file: $referencefastafile\n";
      print $loghandle "  Check for duplicated regions: " . ( $finddupregions ? "Yes" : "No" ) . "\n";
      print $loghandle "  External genomes: ";
      if( $findexternalfastas ){ print $loghandle "\n    " . join( "\n    ", @fastafilelist ) . "\n"; } else { print $loghandle "none\n"; }
      print "Fasta data:\n  Reference file: $referencefastafile\n";
      print "  Check for duplicated regions: " . ( $finddupregions ? "Yes" : "No" ) . "\n";
      print "  External genomes: ";
      if( $findexternalfastas ){ print "\n    " . join( "\n    ", @fastafilelist ) . "\n"; } else { print "none\n"; }
      print $loghandle "Filter configuration:\n";
      print $loghandle "  Minimum coverage: " . ( $mincoverage ? "$mincoverage\n" : "Off" );
      print $loghandle "  Minimum proportion: " . ( $minproportion ? "$minproportion\n" : "Off" );
      print "Filter configuration:\n";
      print "  Minimum coverage: " . ( $mincoverage ? "$mincoverage\n" : "Off" );
      print "  Minimum proportion: " . ( $minproportion ? "$minproportion\n" : "Off" );
      print $loghandle "Aligners:\n";
      print $loghandle "  Aligners to run: " . ( ( $runbwa + $runbwamem + $runnovo ) ? ( $runbwa ? 'BWA ' : '' ) . ( $runbwamem ? 'BWA-mem ' : '' ) . ( $runnovo ? 'Novoalign ' : '' ) : "none" ) . "\n";
      if( $runbwa || $runbwamem )
      {
        print $loghandle "  BWA:\n    Executable: '$bwapath'\n";
        if( $runbwa ){ print $loghandle "    BWA aln args: '$defaultbwaalnargs'\n    BWA samp/se args: '$defaultbwasampeargs'\n"; }
        if( $runbwamem ){ print $loghandle "    BWA mem args: '$defaultbwamemargs'\n"; }
      }
      if( $runnovo )
      {
        print $loghandle "  Novoalign:\n    Executable: '$novopath'\n    Args: '$defaultnovoargs'\n";
        if( $findpaired ){ print $loghandle "    Paired-read args: '$novopairedargs'\n"; }
      }
      if( $runbwa + $runbwamem + $runnovo ){ print $loghandle "  Check for paired read files: " . ( $findpaired ? "Yes" : "No" ) . "\n"; }
      print $loghandle "  Fastq files: ";
      if( $runbwa + $runbwamem + $runnovo )
      {
        print $loghandle "\n";
        foreach my $fastqpair (@fastqfilelist)
        {
          if( scalar( @{$fastqpair} ) == 2 )
          {
            print $loghandle "    / $fastqpair->[0]\n";
            print $loghandle "    \\ $fastqpair->[1]\n";
          } else { print $loghandle "    $fastqpair->[0]\n"; }
        }
      } else { print $loghandle "none\n"; }
      print $loghandle "  Pre-aligned files: ";
      if( $findbams ){ print $loghandle "\n    " . join( "\n    ", @bamfilelist ) . "\n"; } else { print $loghandle "none\n"; }
      print "Aligners:\n";
      print "  Aligners to run: " . ( ( $runbwa + $runbwamem + $runnovo ) ? ( $runbwa ? 'BWA ' : '' ) . ( $runbwamem ? 'BWA-mem ' : '' ) . ( $runnovo ? 'Novoalign ' : '' ) : "none" ) . "\n";
      if( $runbwa || $runbwamem )
      {
        print "  BWA:\n    Executable: '$bwapath'\n";
        if( $runbwa ){ print "    BWA aln args: '$defaultbwaalnargs'\n    BWA samp/se args: '$defaultbwasampeargs'\n"; }
        if( $runbwamem ){ print "    BWA mem args: '$defaultbwamemargs'\n"; }
      }
      if( $runnovo )
      {
        print "  Novoalign:\n    Executable: '$novopath'\n    Args: '$defaultnovoargs'\n";
        if( $findpaired ){ print "    Paired-read args: '$novopairedargs'\n"; }
      }
      if( $runbwa + $runbwamem + $runnovo ){ print "  Check for paired read files: " . ( $findpaired ? "Yes" : "No" ) . "\n"; }
      print "  Fastq files: ";
      if( $runbwa + $runbwamem + $runnovo )
      {
        print "\n";
        foreach my $fastqpair (@fastqfilelist)
        {
          if( scalar( @{$fastqpair} ) == 2 )
          {
            print "    / $fastqpair->[0]\n";
            print "    \\ $fastqpair->[1]\n";
          } else { print "    $fastqpair->[0]\n"; }
        }
      } else { print "none\n"; }
      print "  Pre-aligned files: ";
      if( $findbams ){ print "\n    " . join( "\n    ", @bamfilelist ) . "\n"; } else { print "none\n"; }
      print $loghandle "SNP callers:\n";
      print $loghandle "  SNP callers to run: " . ( ( $rungatk + $runsolsnp + $runvarscan + $runsamtools ) ? ( $rungatk ? 'GATK ' : '' ) . ( $runsolsnp ? 'SolSNP ' : '' ) . ( $runvarscan ? 'VarScan ' : '' ) . ( $runsamtools ? 'SAMtools ' : '' ) : "none" ) . "\n";
      if( $rungatk ){ print $loghandle "  GATK:\n    Executable: '$gatkpath'\n    Args: '$defaultgatkargs'\n"; }
      if( $runsolsnp ){ print $loghandle "  SolSNP:\n    Executable: '$solsnppath'\n    Args: '$defaultsolsnpargs'\n"; }
      if( $runvarscan ){ print $loghandle "  VarScan:\n    Executable: '$varscanpath'\n    Args: '$defaultvarscanargs'\n"; }
      if( $runsamtools ){ print $loghandle "  SAMtools:\n    Executable: '$bcftoolspath'\n"; }
      print $loghandle "  Pre-called VCF files: ";
      if( $findvcfs ){ print $loghandle "\n    " . join( "\n    ", @vcffilelist ) . "\n"; } else { print $loghandle "none\n"; }
      print "SNP callers:\n";
      print "  SNP callers to run: " . ( ( $rungatk + $runsolsnp + $runvarscan + $runsamtools ) ? ( $rungatk ? 'GATK ' : '' ) . ( $runsolsnp ? 'SolSNP ' : '' ) . ( $runvarscan ? 'VarScan ' : '' ) . ( $runsamtools ? 'SAMtools ' : '' ) : "none" ) . "\n";
      if( $rungatk ){ print "  GATK:\n    Executable: '$gatkpath'\n    Args: '$defaultgatkargs'\n"; }
      if( $runsolsnp ){ print "  SolSNP:\n    Executable: '$solsnppath'\n    Args: '$defaultsolsnpargs'\n"; }
      if( $runvarscan ){ print "  VarScan:\n    Executable: '$varscanpath'\n    Args: '$defaultvarscanargs'\n"; }
      if( $runsamtools ){ print "  SAMtools:\n    Executable: '$bcftoolspath'\n"; }
      print "  Pre-called VCF files: ";
      if( $findvcfs ){ print "\n    " . join( "\n    ", @vcffilelist ) . "\n"; } else { print "none\n"; }
      print $loghandle "\nCommands submitted:\n\n";

      # This is the main core of the pipeline.
      # This section prepares the jobs, and then submits them to PBS.
      # From there, this script exits and relies on the downstream components to finish the job.
      if( !( -e( $outputfilefolder ) ) ){ mkdir( $outputfilefolder ); }
      if( !( -e( "$outputfilefolder/reference" ) ) ){ mkdir( "$outputfilefolder/reference" ); }
      `ln -s -f $referencefastafile $outputfilefolder/reference/reference.fasta`;
      $referencefastafile = "$outputfilefolder/reference/reference.fasta";
      if( scalar( @fastafilelist ) && $findexternalfastas ){ mkdir( "$outputfilefolder/external" ); }
      if( scalar( @fastqfilelist ) )
      {
        if( $runbwa && !( -e( "$outputfilefolder/bwa" ) ) ){ mkdir( "$outputfilefolder/bwa" ); }
        if( $runbwamem && !( -e( "$outputfilefolder/bwamem" ) ) ){ mkdir( "$outputfilefolder/bwamem" ); }
        if( $runnovo && !( -e( "$outputfilefolder/novoalign" ) ) ){ mkdir( "$outputfilefolder/novoalign" ); }
      }
      if( scalar( @bamfilelist ) && !( -e( "$outputfilefolder/bams" ) ) ){ mkdir( "$outputfilefolder/bams" ); }
      if( ( scalar( @fastqfilelist ) + scalar( @bamfilelist ) ) >= 1 )
      {
        if( $rungatk && !( -e( "$outputfilefolder/gatk" ) ) ){ mkdir( "$outputfilefolder/gatk" ); }
        if( $runsolsnp && !( -e( "$outputfilefolder/solsnp" ) ) ){ mkdir( "$outputfilefolder/solsnp" ); }
        if( $runvarscan && !( -e( "$outputfilefolder/varscan" ) ) ){ mkdir( "$outputfilefolder/varscan" ); }
        if( $runsamtools && !( -e( "$outputfilefolder/samtools" ) ) ){ mkdir( "$outputfilefolder/samtools" ); }
      }
  
      my $indexcommand = "";
      if( ( $runbwa || $runbwamem ) && ( scalar( @fastqfilelist ) >= 1 ) ){ $indexcommand .= "$bwapath index $referencefastafile \n"; }
      if( $runnovo && ( scalar( @fastqfilelist ) >= 1 ) ){ $indexcommand .= "$novoindexpath $referencefastafile.idx $referencefastafile \n"; }
      if( $rungatk && ( ( scalar( @fastqfilelist ) >= 1 ) || ( scalar( @bamfilelist ) >= 1 ) ) )
      {
        $indexcommand .= "java -Xmx${gigsofmemforindex}G -jar $dictgeneratorpath R=$referencefastafile O=$referencefastafile.dict \n";
        if( $referencefastafile =~ /^(?:.*\/)?([^\/]+?)\.([^.]+)$/ ){ `ln -s -f $1.$2.dict $outputfilefolder/reference/$1.dict`; }
        $indexcommand .= "$samtoolspath faidx $referencefastafile \n";
      }
      if( $finddupregions )
      {
        $indexcommand .= "$finddupspath $referencefastafile duplicates.txt \n";
        push( @finalfilelist, "dups,nucmer,::$outputfilefolder/reference/duplicates.txt" );
      }
      if( scalar( @bamfilelist ) >= 1 )
      {
        my @newbamfilelist = ();
        foreach my $bamfilename (@bamfilelist)
        {
          my $bamfilenickname = $bamfilename;
          if( $bamfilename =~ /^(?:.*\/)?([^\/]+?)(?:\.[BbSs][Aa][Mm])?$/ ){ $bamfilenickname = $1; }
          `ln -s -f $readfilefolder/$bamfilename $outputfilefolder/bams/$bamfilenickname.bam`;
          $indexcommand .= "$samtoolspath index $outputfilefolder/bams/$bamfilenickname.bam \n";
          push( @newbamfilelist, "$bamfilenickname.bam" );
        }
        @bamfilelist = @newbamfilelist;
      }
      my @snpcallerqids = ();
      my $possiblereadfile = "";
      my $pipelinestartqid = `echo "sleep 2 \n $indexcommand" | qsub -d '$outputfilefolder/reference' -w '$outputfilefolder/reference' -l ncpus=1,mem=${gigsofmemforindex}gb,walltime=$wallhoursforindex:00:00 -m ab -N 'nasp_index' -x - `;
      chomp( $pipelinestartqid );
      if( $pipelinestartqid =~ /^\d{1,31}\.\w/ )
      {
        print $loghandle "$pipelinestartqid:\n$indexcommand\n";
        if( scalar( @fastafilelist ) >= 1 )
        {
          foreach my $fastafilename (@fastafilelist)
          {
            my $fastafilenickname = $fastafilename;
            if( $fastafilename =~ /^(?:.*\/)?([^\/]+?)(?:\.[Ff][Aa](?:[Ss](?:[Tt][Aa])?)?)?$/ ){ $fastafilenickname = $1; }
            my $externalfastadata = _submit_external_fasta( $pipelinestartqid, $referencefastafile, $readfilefolder, $outputfilefolder, $fastafilename, $fastafilenickname, $loghandle );
            if( scalar( @{$externalfastadata} ) == 2 )
            {
              push( @snpcallerqids, $externalfastadata->[0] );
              push( @finalfilelist, ( "external,nucmer,::" . $externalfastadata->[1] ) );
            }
          }
        }
        if( scalar( @fastqfilelist ) >= 1 )
        {
          foreach my $readfilepair (@fastqfilelist)
          {
            my $readfilenickname = $readfilepair->[0];
            if( $readfilepair->[0] =~ /(?:.*\/)?([^\/]+?)(?:_[12]_)?sequence\.txt(?:\.gz)?$/i ){ $readfilenickname = $1; }
            elsif( $readfilepair->[0] =~ /(?:.*\/)?([^\/]+?)\.f(?:ast)?q(?:\.gz)?$/i )
            {
              $readfilenickname = $1;
              $readfilenickname =~ s/_R[12]_/_/g;
            }
            if( $runbwa )
            {
              my $alignerdata = _submit_bwa( $pipelinestartqid, $referencefastafile, $readfilefolder, $outputfilefolder, $readfilepair, $readfilenickname, $loghandle );
              if( scalar( @{$alignerdata} ) == 2 )
              {
                if( $rungatk )
                {
                  my $snpcallerdata = _submit_gatk( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwa", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwa", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA,GATK,::" . $snpcallerdata->[1] ) );
                }
                if( $runsolsnp )
                {
                  `ln -s -f $alignerdata->[1] $outputfilefolder/bwa/$readfilenickname`;
                  my $snpcallerdata = _submit_solsnp( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwa", $outputfilefolder, $readfilenickname, "$readfilenickname-bwa", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA,SolSnp,::" . $snpcallerdata->[1] ) );
                }
                if( $runvarscan )
                {
                  my $snpcallerdata = _submit_varscan( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwa", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwa", $readfilenickname, $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA,VarScan,::" . $snpcallerdata->[1] ) );
                }
                if( $runsamtools )
                {
                  my $snpcallerdata = _submit_samtools( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwa", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwa", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA,SAMtools,::" . $snpcallerdata->[1] ) );
                }
              }
            }
            if( $runbwamem )
            {
              my $alignerdata = _submit_bwamem( $pipelinestartqid, $referencefastafile, $readfilefolder, $outputfilefolder, $readfilepair, $readfilenickname, $loghandle );
              if( scalar( @{$alignerdata} ) == 2 )
              {
                if( $rungatk )
                {
                  my $snpcallerdata = _submit_gatk( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwamem", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwamem", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA-mem,GATK,::" . $snpcallerdata->[1] ) );
                }
                if( $runsolsnp )
                {
                  `ln -s -f $alignerdata->[1] $outputfilefolder/bwamem/$readfilenickname`;
                  my $snpcallerdata = _submit_solsnp( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwamem", $outputfilefolder, $readfilenickname, "$readfilenickname-bwamem", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA-mem,SolSnp,::" . $snpcallerdata->[1] ) );
                }
                if( $runvarscan )
                {
                  my $snpcallerdata = _submit_varscan( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwamem", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwamem", $readfilenickname, $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA-mem,VarScan,::" . $snpcallerdata->[1] ) );
                }
                if( $runsamtools )
                {
                  my $snpcallerdata = _submit_samtools( $alignerdata->[0], $referencefastafile, "$outputfilefolder/bwamem", $outputfilefolder, $alignerdata->[1], "$readfilenickname-bwamem", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,BWA-mem,SAMtools,::" . $snpcallerdata->[1] ) );
                }
              }
            }
            if( $runnovo )
            {
              my $alignerdata = _submit_novo( $pipelinestartqid, $referencefastafile, $readfilefolder, $outputfilefolder, $readfilepair, $readfilenickname, $loghandle );
              if( scalar( @{$alignerdata} ) == 2 )
              {
                if( $rungatk )
                {
                  my $snpcallerdata = _submit_gatk( $alignerdata->[0], $referencefastafile, "$outputfilefolder/novoalign", $outputfilefolder, $alignerdata->[1], "$readfilenickname-novo", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,Novoalign,GATK,::" . $snpcallerdata->[1] ) );
                }
                if( $runsolsnp )
                {
                  `ln -s -f $alignerdata->[1] $outputfilefolder/novoalign/$readfilenickname`;
                  my $snpcallerdata = _submit_solsnp( $alignerdata->[0], $referencefastafile, "$outputfilefolder/novoalign", $outputfilefolder, $readfilenickname, "$readfilenickname-novo", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,Novoalign,SolSnp,::" . $snpcallerdata->[1] ) );
                }
                if( $runvarscan )
                {
                  my $snpcallerdata = _submit_varscan( $alignerdata->[0], $referencefastafile, "$outputfilefolder/novoalign", $outputfilefolder, $alignerdata->[1], "$readfilenickname-novo", $readfilenickname, $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,Novoalign,VarScan,::" . $snpcallerdata->[1] ) );
                }
                if( $runsamtools )
                {
                  my $snpcallerdata = _submit_samtools( $alignerdata->[0], $referencefastafile, "$outputfilefolder/novoalign", $outputfilefolder, $alignerdata->[1], "$readfilenickname-novo", $loghandle );
                  push( @snpcallerqids, $snpcallerdata->[0] );
                  push( @finalfilelist, ( "vcf,Novoalign,SAMtools,::" . $snpcallerdata->[1] ) );
                }
              }
            }
          }
        }
        if( scalar( @bamfilelist ) >= 1 )
        {
          foreach my $bamfilename (@bamfilelist)
          {
            my $bamfilenickname = $bamfilename;
            if( $bamfilename =~ /^(?:.*\/)?([^\/]+?)(?:\.[BbSs][Aa][Mm])?$/ ){ $bamfilenickname = $1; }
            if( $rungatk )
            {
              my $snpcallerdata = _submit_gatk( $pipelinestartqid, $referencefastafile, "$outputfilefolder/bams", $outputfilefolder, $bamfilename, $bamfilenickname, $loghandle );
              if( scalar( @{$snpcallerdata} ) == 2 )
              {
                push( @snpcallerqids, $snpcallerdata->[0] );
                push( @finalfilelist, ( "vcf,pre-aligned,GATK,::" . $snpcallerdata->[1] ) );
              }
            }
            if( $runsolsnp )
            {
              `ln -s -f $readfilefolder/$bamfilename $outputfilefolder/bams/$bamfilenickname`;
              my $snpcallerdata = _submit_solsnp( $pipelinestartqid, $referencefastafile, "$outputfilefolder/bams", $outputfilefolder, $bamfilenickname, $bamfilenickname, $loghandle );
              if( scalar( @{$snpcallerdata} ) == 2 )
              {
                push( @snpcallerqids, $snpcallerdata->[0] );
                push( @finalfilelist, ( "vcf,pre-aligned,SolSnp,::" . $snpcallerdata->[1] ) );
              }
            }
            if( $runvarscan )
            {
              my $snpcallerdata = _submit_varscan( $pipelinestartqid, $referencefastafile, "$outputfilefolder/bams", $outputfilefolder, $bamfilename, $bamfilenickname, $bamfilenickname, $loghandle );
              if( scalar( @{$snpcallerdata} ) == 2 )
              {
                push( @snpcallerqids, $snpcallerdata->[0] );
                push( @finalfilelist, ( "vcf,pre-aligned,VarScan,::" . $snpcallerdata->[1] ) );
              }
            }
            if( $runsamtools )
            {
              my $snpcallerdata = _submit_samtools( $pipelinestartqid, $referencefastafile, "$outputfilefolder/bams", $outputfilefolder, $bamfilename, $bamfilenickname, $loghandle );
              if( scalar( @{$snpcallerdata} ) == 2 )
              {
                push( @snpcallerqids, $snpcallerdata->[0] );
                push( @finalfilelist, ( "vcf,pre-aligned,SAMtools,::" . $snpcallerdata->[1] ) );
              }
            }
          }
        }
        if( scalar( @finalfilelist ) > 0 )
        {
          my $jobidstowaitfor = join( ':', @snpcallerqids );
          my $finalfilestring = join( ' ', @finalfilelist );
          if( defined( $mastermatrixfile ) && length( $mastermatrixfile ) )
          {
            my $commandtorun = "$distancecalcscript $mincoverage $minproportion $mastermatrixfile $finalfilestring $outputfilefolder/snp_matrix.tsv $outputfilefolder/distance_matrix.tsv $outputfilefolder/filter_log.txt \n";
            my $distancecalcqid = `echo "$commandtorun" | qsub -d '$outputfilefolder' -w '$outputfilefolder' -l ncpus=1,mem=${gigsofmemfordistancecalc}gb,walltime=$wallhoursfordistancecalc:00:00 -m ae -N 'nasp_matrix' -W depend=afterok:$pipelinestartqid -W depend=afterany:$jobidstowaitfor -x - `;
            chomp( $distancecalcqid );
            print $loghandle "$distancecalcqid:\n$commandtorun\n";
            print "\nThe pipeline has been submitted to PBS for batch execution.\nResults can be found in '$outputfilefolder/distance_matrix.tsv' when job '$distancecalcqid' is complete.\n";
          } else
          {
            my $commandtorun = "$matrixmakingscript $mincoverage $minproportion $referencefastafile $finalfilestring $outputfilefolder/allvariant_matrix.tsv $outputfilefolder/bestsnps_matrix.tsv $outputfilefolder/allcallable_matrix.tsv $outputfilefolder/bestsnps.snpfasta $outputfilefolder/statistics.tsv \n";
            my $matrixmakingqid = `echo "$commandtorun" | qsub -d '$outputfilefolder' -w '$outputfilefolder' -l ncpus=1,mem=${gigsofmemtomakematrix}gb,walltime=$wallhourstomakematrix:00:00 -m ae -N 'nasp_matrix' -W depend=afterok:$pipelinestartqid -W depend=afterany:$jobidstowaitfor -x - `;
            chomp( $matrixmakingqid );
            print $loghandle "$matrixmakingqid:\n$commandtorun\n";
            print "\nThe pipeline has been submitted to PBS for batch execution.\nResults can be found in '$outputfilefolder/snp_matrix.tsv' when job '$matrixmakingqid' is complete.\n";
          }
        } else { print STDERR "No call files found for matrix transformation!\n"; }
      } else { print STDERR "Unable to schedule jobs!\n$pipelinestartqid\n"; }

      close( $loghandle );
    } else { print STDERR "Could not find any files to analyze in '$readfilefolder'!\n"; }
  } else { print STDERR "Could not open '$readfilefolder'!\n"; }

}

sub _submit_external_fasta
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $fastafilename = shift();
  my $fastafilenickname = shift();
  my $loghandle = shift();
  my $commandtorun = "$convertexternalpath $referencefastafile $inputfilefolder/$fastafilename $fastafilenickname.frankenfasta \n";
  my $snpcallerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/external' -w '$outputfilefolder/external' -l ncpus=1,mem=${gigsofmemforexternal}gb,walltime=$wallhoursforexternal:00:00 -m a -N 'nasp_external_$fastafilename' -x -W depend=afterok:$jobtodependon - `;
  chomp( $snpcallerqid );
  print $loghandle "$snpcallerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $snpcallerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $snpcallerqid );
    push( @returnarray, "$outputfilefolder/external/$fastafilenickname.frankenfasta" );
  } else { print STDERR "Error scheduling external genome conversion for '$fastafilename'!\n"; }
  return( \@returnarray );
}

sub _submit_bwa
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $readfilepair = shift();
  my $readfilenickname = shift();
  my $loghandle = shift();
  my $ispaired = 0;
  if( scalar( @{$readfilepair} ) == 2 ){ $ispaired = 1; }
  my $oldformatstring = "";
  if( $readfilepair->[0] =~ /(?:.*\/)?[^\/]+?_[12]_sequence\.txt(?:\.gz)?$/i ){ $oldformatstring = "-I"; }
  my $bamstringthingy = '@RG\tID:' . $readfilenickname . '\tSM:' . $readfilenickname;
  my $commandtorun = "";
  if( $ispaired )
  {
    $commandtorun .= "$bwapath aln $oldformatstring $referencefastafile $inputfilefolder/$readfilepair->[0] -t $numcpusforbwa -f $outputfilefolder/bwa/$readfilenickname-R1.sai $defaultbwaalnargs \n";
    $commandtorun .= "$bwapath aln $oldformatstring $referencefastafile $inputfilefolder/$readfilepair->[1] -t $numcpusforbwa -f $outputfilefolder/bwa/$readfilenickname-R2.sai $defaultbwaalnargs \n";
    $commandtorun .= "$bwapath sampe -r '$bamstringthingy' $referencefastafile $outputfilefolder/bwa/$readfilenickname-R1.sai $outputfilefolder/bwa/$readfilenickname-R2.sai $inputfilefolder/$readfilepair->[0] $inputfilefolder/$readfilepair->[1] $defaultbwasampeargs | $samtoolspath view -S -F 4 -b -h -q 5 - | $samtoolspath sort - $readfilenickname-bwa \n";
  } else
  {
    $commandtorun .= "$bwapath aln $oldformatstring $referencefastafile $inputfilefolder/$readfilepair->[0] -t $numcpusforbwa -f $outputfilefolder/bwa/$readfilenickname.sai $defaultbwaalnargs \n";
    $commandtorun .= "$bwapath samse -r '$bamstringthingy' $referencefastafile $outputfilefolder/bwa/$readfilenickname.sai $inputfilefolder/$readfilepair->[0] $defaultbwasampeargs | $samtoolspath view -S -F 4 -b -h -q 5 - | $samtoolspath sort - $readfilenickname-bwa \n";
  }
  $commandtorun .= "$samtoolspath index $readfilenickname-bwa.bam \n";
  my $alignerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/bwa' -w '$outputfilefolder/bwa' -l ncpus=$numcpusforbwa,mem=${gigsofmemforbwa}gb,walltime=$wallhoursforbwa:00:00 -m a -N 'nasp_bwa_$readfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $alignerqid );
  print $loghandle "$alignerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $alignerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $alignerqid );
    push( @returnarray, "$readfilenickname-bwa.bam" );
  } else { print STDERR "Error scheduling BWA for '$readfilepair->[0]'!\n"; }
  return( \@returnarray );
}

sub _submit_bwamem
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $readfilepair = shift();
  my $readfilenickname = shift();
  my $loghandle = shift();
  my $oldformatstring = "";
  my $bamstringthingy = '@RG\tID:' . $readfilenickname . '\tSM:' . $readfilenickname;
  my $readfilestring = "$inputfilefolder/" . join( " $inputfilefolder/", @{$readfilepair} );
  my $commandtorun = "$bwapath mem -R '$bamstringthingy' $defaultbwamemargs -t $numcpusforbwa $referencefastafile $readfilestring | $samtoolspath view -S -F 4 -b -h -q 5 - | $samtoolspath sort - $readfilenickname-bwamem \n $samtoolspath index $readfilenickname-bwamem.bam \n";
  my $alignerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/bwamem' -w '$outputfilefolder/bwamem' -l ncpus=$numcpusforbwa,mem=${gigsofmemforbwa}gb,walltime=$wallhoursforbwa:00:00 -m a -N 'nasp_bwamem_$readfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $alignerqid );
  print $loghandle "$alignerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $alignerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $alignerqid );
    push( @returnarray, "$readfilenickname-bwamem.bam" );
  } else { print STDERR "Error scheduling BWA for '$readfilepair->[0]'!\n"; }
  return( \@returnarray );
}

sub _submit_novo
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $readfilepair = shift();
  my $readfilenickname = shift();
  my $loghandle = shift();
  my $ispaired = 0;
  if( scalar( @{$readfilepair} ) == 2 ){ $ispaired = 1; }
  my $bamstringthingy = '@RG\tID:' . $readfilenickname . '\tSM:' . $readfilenickname;
  my $readfilestring = "$inputfilefolder/" . join( " $inputfilefolder/", @{$readfilepair} );
  my $pairedstring = ( ( $ispaired ) ? "$novopairedargs" : "" );
  my $commandtorun = "$novopath -f $readfilestring $pairedstring -c $numcpusfornovo -o SAM '$bamstringthingy' -d $referencefastafile.idx $defaultnovoargs | $samtoolspath view -S -F 4 -b -h -q 5 - | $samtoolspath sort - $readfilenickname-novo \n $samtoolspath index $readfilenickname-novo.bam \n";
  my $alignerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/novoalign' -w '$outputfilefolder/novoalign' -l ncpus=$numcpusfornovo,mem=${gigsofmemfornovo}gb,walltime=$wallhoursfornovo:00:00 -m a -N 'nasp_novo_$readfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $alignerqid );
  print $loghandle "$alignerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $alignerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $alignerqid );
    push( @returnarray, "$readfilenickname-novo.bam" );
  } else { print STDERR "Error scheduling Novoalign for '$readfilepair->[0]'!\n"; }
  return( \@returnarray );
}

sub _submit_gatk
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $bamfilename = shift();
  my $bamfilenickname = shift();
  my $loghandle = shift();
  my $commandtorun = "java -Xmx${gigsofmemforgatk}G -jar $gatkpath -T UnifiedGenotyper -dt NONE -I $inputfilefolder/$bamfilename -R $referencefastafile -nt $numcpusforgatk -ploidy 1 -o $outputfilefolder/gatk/$bamfilenickname-gatk.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE $defaultgatkargs \n";
  my $snpcallerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/gatk' -w '$outputfilefolder/gatk' -l ncpus=$numcpusforgatk,mem=${gigsofmemforgatk}gb,walltime=$wallhoursforgatk:00:00 -m a -N 'nasp_gatk_$bamfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $snpcallerqid );
  print $loghandle "$snpcallerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $snpcallerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $snpcallerqid );
    push( @returnarray, "$outputfilefolder/gatk/$bamfilenickname-gatk.vcf" );
  } else { print STDERR "Error scheduling GATK for '$bamfilename'!\n"; }
  return( \@returnarray );
}

sub _submit_solsnp
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $bamfilename = shift();
  my $bamfilenickname = shift();
  my $loghandle = shift();
  my $commandtorun = "java -Xmx${gigsofmemforsolsnp}G -jar $solsnppath INPUT=$inputfilefolder/$bamfilename REFERENCE_SEQUENCE=$referencefastafile OUTPUT=$outputfilefolder/solsnp/$bamfilenickname-solsnp.vcf SUMMARY=true CALCULATE_ALLELIC_BALANCE=true MINIMUM_COVERAGE=1 PLOIDY=Haploid STRAND_MODE=None OUTPUT_FORMAT=VCF OUTPUT_MODE=AllCallable $defaultsolsnpargs \n";
  my $snpcallerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/solsnp' -w '$outputfilefolder/solsnp' -l ncpus=$numcpusforsolsnp,mem=${gigsofmemforsolsnp}gb,walltime=$wallhoursforsolsnp:00:00 -m a -N 'nasp_solsnp_$bamfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $snpcallerqid );
  print $loghandle "$snpcallerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $snpcallerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $snpcallerqid );
    push( @returnarray, "$outputfilefolder/solsnp/$bamfilenickname-solsnp.vcf" );
  } else { print STDERR "Error scheduling SolSNP for '$bamfilename'!\n"; }
  return( \@returnarray );
}

sub _submit_varscan
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $bamfilename = shift();
  my $bamfilenickname = shift();
  my $shortenednickname = shift(); # it's in the bam file retard
  my $loghandle = shift();
  my $commandtorun = "echo '$shortenednickname' > $outputfilefolder/varscan/$bamfilenickname.txt \n $samtoolspath mpileup -B -d 10000000 -f $referencefastafile $inputfilefolder/$bamfilename > $inputfilefolder/$bamfilenickname.mpileup \n java -Xmx${gigsofmemforvarscan}G -jar $varscanpath mpileup2cns $inputfilefolder/$bamfilenickname.mpileup --output-vcf 1 --vcf-sample-list $outputfilefolder/varscan/$bamfilenickname.txt > $outputfilefolder/varscan/$bamfilenickname-varscan.vcf $defaultvarscanargs \n";
  my $snpcallerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/varscan' -w '$outputfilefolder/varscan' -l ncpus=$numcpusforvarscan,mem=${gigsofmemforvarscan}gb,walltime=$wallhoursforvarscan:00:00 -m a -N 'nasp_varscan_$bamfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $snpcallerqid );
  print $loghandle "$snpcallerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $snpcallerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $snpcallerqid );
    push( @returnarray, "$outputfilefolder/varscan/$bamfilenickname-varscan.vcf" );
  } else { print STDERR "Error scheduling VarScan for '$bamfilename'!\n"; }
  return( \@returnarray );
}

sub _submit_samtools
{
  my $jobtodependon = shift();
  my $referencefastafile = shift();
  my $inputfilefolder = shift();
  my $outputfilefolder = shift();
  my $bamfilename = shift();
  my $bamfilenickname = shift();
  my $loghandle = shift();
  my $commandtorun = "$samtoolspath mpileup -uD -d 10000000 -f $referencefastafile $inputfilefolder/$bamfilename | $bcftoolspath view -ceg - > $outputfilefolder/samtools/$bamfilenickname-samtools.vcf \n";
  my $snpcallerqid = `echo "$commandtorun" | qsub -d '$outputfilefolder/samtools' -w '$outputfilefolder/samtools' -l ncpus=$numcpusforbcftools,mem=${gigsofmemforbcftools}gb,walltime=$wallhoursforbcftools:00:00 -m a -N 'nasp_samtools_$bamfilenickname' -x -W depend=afterok:$jobtodependon - `;
  chomp( $snpcallerqid );
  print $loghandle "$snpcallerqid:\n$commandtorun\n";
  my @returnarray = ();
  if( $snpcallerqid =~ /^\d{1,31}\.\w/ )
  {
    push( @returnarray, $snpcallerqid );
    push( @returnarray, "$outputfilefolder/samtools/$bamfilenickname-samtools.vcf" );
  } else { print STDERR "Error scheduling SAMtools for '$bamfilename'!\n"; }
  return( \@returnarray );
}


