assemblies, = glob_wildcards(expand("{assemblies_dir}/{{id}}.fasta", assemblies_dir=config['assemblies'])[0])
reads, = glob_wildcards('reads/{id}_1.fq')

rule matrix:
  params:
    minimum_coverage=config['minimum_coverage'],
    minimum_proportion=config['minimum_proportion']
  input:
    reference=config['reference'],
    frankenfasta=expand('frankenfasta/{id}.frankenfasta', id=assemblies)
    #vcf=expand('gatk4/{id}.vcf', id=reads)
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

# TODO: mark fasta as temporary
rule iqtree:
  input: rules.matrix.output.bestsnp
  conda: "envs/iqtree.yaml"
  shell: """
  nasp export --type fasta {input} > {input}.fasta
  iqtree -nt {threads} -s {input}.fasta
  """
