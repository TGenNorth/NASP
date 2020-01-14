assemblies, = glob_wildcards(expand("{assemblies_dir}/{{id}}.fasta", assemblies_dir=config['assemblies'])[0])
reads, = glob_wildcards('reads/{id}_1.fq')

# input function for the rule aggregate
def aggregate_input(wildcards):
  # decision based on content of output file
  # Important: use the method open() of the returned file!
  # This way, Snakemake is able to automatically download the file if it is generated in
  # a cloud environment without a shared filesystem.
  with checkpoints.frankenfasta.get(sample=wildcards.sample).output[0].open() as f:
    if f.read().strip() == "a":
      return "post/{sample}.txt"
    else:
      return "alt/{sample}.txt"

rule matrix:
  """
  The NASP matrix aggregates VCF and Frankenfasta files into .tsv matrix and stats file summaries.
  The primary output is the 'bestsnp.tsv' representing high-confidence SNPs across all samples.
  """
  params:
    minimum_coverage=config['minimum_coverage'],
    minimum_proportion=config['minimum_proportion']
  input:
    reference=config['reference'],
    frankenfasta=expand('frankenfasta/{id}.frankenfasta', id=assemblies)
    #vcf=expand('gatk4/{id}.vcf', id=reads)
  params:
    prefix=lambda wildcards, output: output[0][:-4]
  output:
    general_stats="general_stats.tsv",
    sample_stats="sample_stats.tsv",
    bestsnp="bestsnp.tsv",
    master="master.tsv",
    missingdata="missingdata.tsv"
  conda: "envs/mummer.yaml"
  shell: """
    nasp matrix \
      --reference-fasta {input.reference} \
      --minimum-coverage {params.minimum_coverage} \
      --minimum-proportion {params.minimum_proportion} \
      frankenfasta/*.frankenfasta
    """

rule iqtree:
  """
  Generate a Maximum Parsimony tree from the NASP 'bestsnp' matrix.

  The primary output is the 'bestsnp.treefile'
  """
  input: rules.matrix.output.bestsnp
  output:
    temp('bestsnp'),
    expand('bestsnp.{ext}', ext=[
      'bionj',
      'iqtree',
      'log',
      'mldist',
      'model.gz',
      'treefile',
      'ckp.gz'
    ])
  # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads
  # FIXME: iqtree has been observed to abort if given too many threads for too trivial an input.
  # TODO: config parameterize iqtree
  threads: workflow.cores * 0.75
  conda: "envs/iqtree.yaml"
  shell: """
  nasp export --type fasta {input} > bestsnp
  iqtree -nt {threads} -s bestsnp
  """
