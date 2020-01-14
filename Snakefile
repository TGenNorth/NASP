#PE_READS, = glob_wildcards(expand('{pe_reads}/{{sample_name}}_R1_001.fastq.gz', pe_reads=config['pe_reads']))
#SE_READS, = glob_wildcards(expand('{se_reads}/{{sample_name}}.fastq.gz', se_reads=config['se_reads']))

from snakemake.utils import min_version, validate

# 5.4.0
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
# 5.6.0
# Add --default-resources flag, that allows to define default resources for jobs (e.g. mem_mb, disk_mb), see docs.
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources
min_version("5.6.0")

# https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration
#
# Schema validation of the config serves to:
# * validate user-input configuration
# * define default values
#
# The `configfile: [FILE]` directive is not used because it requires the file exists in the working directory.
# `snakemake --configfile [FILE]` allows the user to override values in the config, but not the path to the config.
#configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

#cells = pd.read_csv(config["cells"], sep="\t").set_index("id", drop=False)
#validate(cells, schema="schemas/cells.schema.yaml")

# https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#combining-conda-package-management-with-containers
# this container defines the underlying OS for each job when using the workflow
# with snakemake --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3:latest"

# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#constraining-wildcards
# "sorted_reads/{sample,[A-Za-z0-9]+}.bam"
#assemblies, = glob_wildcards(expand("{assemblies_dir}/{{id}}.fasta", assemblies_dir=config['assemblies_dir']))
#assemblies = expand("{assemblies_dir}/{{id}}.fasta", assemblies_dir=config['assemblies_dir'])
#ids, = glob_wildcards("thedir/{id}.fastq.gz")

docstring= r"""
TODO: Help for this script

./pe_reads/
./se_reads/
./ubam/
./assemblies/
./reference.fasta

USAGE

snamekemake -s /path/to/workflow/Snakefile -j $(nproc) --use-conda iqtree
"""

# https://www.biostars.org/p/220268/#220632
# Tutorial: How to make your snakefiles self-documenting
rule help:
  """
  Print list of all targets with help.
  """
  run:
    for rule in workflow.rules:
      print( rule.name )
      print( rule.docstring )

#rule all:
#  input: []
#  output: []
#  params: []
#  log: []
#  benchmark: []
#  message: ''
#  threads: 1
#  resources:
#    mem_mb=lambda wildcards, attempt: attempt * 100
#  run:
#    print(config)

rule config:
  """
  Print workflow configuration
  """
  run:
    import json, sys
    json.dump(config, sys.stdout, indent='  ')

# onstart handler, that allows to add code that shall be only executed before the actual workflow execution (not on dryrun).
# Parameters defined in the cluster config file are now accessible in the job properties under the key “cluster”.


# Error: no Snakefile found, tried Snakefile, snakefile, workflow/Snakefile, workflow/snakefile.

##### setup report #####

report: "report/workflow.rst"

##### load rules #####

include: "rules/testdata.smk"
include: "rules/matrix.smk"
include: "rules/frankenfasta.smk"
include: "rules/bwa.smk"
