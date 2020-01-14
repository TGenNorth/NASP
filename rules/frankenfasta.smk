# TODO: what if reference or query is blank/invalid?
rule frankenfasta:
  input:
    reference=config['reference'],
    query=expand('{assemblies_dir}/{{sample_name}}.fasta', assemblies_dir=config['assemblies'])
  output: 'frankenfasta/{sample_name}.frankenfasta'
  params:
    sample_name='{sample_name}'
  log: []
  benchmark: []
  message: ''
  #threads: 1
  #resources: []
  conda: "envs/mummer.yaml"
  # nucmer [options] <Reference> <Query>
  # delta-filter [options] <deltafile>
  #
  # TODO: remove test -s {output}
  # The nasp<=1.1.2 python dispatcher does not forward the nasptool exit code.
  # As a stop-gap, fail if the output is empty.
  shell: '''
  nasp frankenfasta <(
    delta-filter -q -r -o 100 \
      <(nucmer --threads {threads} --delta /dev/stdout {input.reference} {input.query})
  ) > {output}
  '''


