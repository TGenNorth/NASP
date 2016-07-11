# NASP Qapla Dispatcher Prototype

The `qapla` branch is a prototype for an alternative implementation of the NASP job dispatcher.

TODO: summarize goal of this implementation

The qapla branch contains 5 files. The important ones are `config.yml` and `qapla`.

- `tmpl/` templates qapla will use when submitting jobs to the job scheduler.
- `common.vars.tmpl` template imported by all template jobs that will contain the list of input files and helper functions.
- `config.yml`configuration where the user or sysadmin can customize the program paths and job submission properties
  to match their system configuration.
- `qapla` the main executable
- `.gitignore`

## Installation Guide:

1. Download the qapla branch: 
  `git clone --depth 1 --branch qapla https://github.com/TGenNorth/NASP qapla`
2. Edit the program paths in `config.yml` so they match the paths on your computer.
  (The first two references to GATK and VarScan are probably not used.)
3. Set the nasp program path in `config.yml` as the absolute path to the qapla executable.

## Usage Guide:

The following are examples of different use-cases.

### Basic usage:

```
qapla  run \
  --aligner bowtie2 \
  --snpcaller gatk \
  --reference /path/to/reference.fasta \
  /path/to/fastq/dir /path/to/fasta/dir /path/to/bam/dir ....
```

### Show the usage message:

```
qapla --help
```

### Show usage message for the run subcommand:

```
qapla run --help
```

### Multiple aligners/snpcallers:

```
qapla  run \
  --aligner bowtie2,bwa,samtools,novoalign \
  --snpcaller gatk,varscan,snap,samtools \
  --reference /path/to/reference.fasta \
  /path/to/fastq/dir /path/to/fasta/dir /path/to/bam/dir ....
```

### Create the Job Scripts and Submit them separately:

```
qapla  run \
  --dry-run \
  --aligner bowtie2 \
  --snpcaller gatk \
  --reference /path/to/reference.fasta \
  /path/to/fastq/dir /path/to/fasta/dir /path/to/bam/dir ....
bash RunJobDispatcher
```


### Run NASP with read trimming and duplicate checking:
Before running this command, edit the following lines in the common.vars.tmpl file:

> trimmomatic_pe_adapter_fasta="${trimmomatic%/*}/adapters/TruSeq3-PE-2.fa"
  trimmomatic_se_adapter_fasta="${trimmomatic%/*}/adapters/TruSeq3-SE.fa"
  
```
qapla  run \
  --duplicates \
  --trim-reads \
  --aligner bowtie2 \
  --snpcaller gatk \
  --reference /path/to/reference.fasta \
  /path/to/fastq/dir /path/to/fasta/dir /path/to/bam/dir ....
```
