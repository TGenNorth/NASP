rule bwa_index:
  input: config['reference'],
  output: expand('{reference}.{ext}', reference=config['reference'], ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
  conda: "envs/bwa.yaml"
  shell: """
  bwa index {input}
  """

# usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
# http://www.htslib.org/workflow/
#fq=expand('{assemblies_dir}/{{sample_name}}.fasta', assemblies_dir=config['assemblies_dir'])
rule bwa:
  input:
    reference=config['reference'],
    bwa_index=rules.bwa_index.output,
    fq='reads/{sample_name}_1.fq',
    fq2='reads/{sample_name}_2.fq',
  output:
    alignment='bwa/{sample_name}.cram',
    alignment_index='bwa/{sample_name}.cram.crai'
  conda: "envs/bwa.yaml"
  shell: """
  bwa mem \
    -t {threads} \
    -R '@RG\\tID:{wildcards.sample_name}\\tSM:{wildcards.sample_name}\\tLB:{wildcards.sample_name}' \
    {input.reference} \
    {input.fq} \
  | samtools sort -O bam -l 0 -T ./ - \
  | samtools view -T {input.reference} -C -o {output.alignment} -
  samtools index {output.alignment}
  """

rule reference_samtools_faidx:
  input: config['reference']
  output: expand('{reference}.fai', reference=config['reference'])
  conda: 'envs/gatk4.yaml'
  shell: 'samtools faidx {input}'

rule reference_picard_sequence_dictionary:
  input: config['reference']
  output: 'reference.dict'
  conda: 'envs/gatk4.yaml'
  shell: """
  picard CreateSequenceDictionary \
    R={input}
    O={output}
  """

rule gatk4:
  input:
    reference=config['reference'],
    reference_samtools_faidx=rules.reference_samtools_faidx.output,
    reference_picard_sequence_dictionary=rules.reference_picard_sequence_dictionary.output,
    alignment='bwa/{sample_name}.cram',
  output: 'gatk4/{sample_name}.vcf'
  conda: "envs/gatk4.yaml"
  shell: """
  gatk \
    HaplotypeCaller \
    -R {input.reference} \
    -I {input.alignment} \
    -O {output} \
    --output-mode EMIT_ALL_CONFIDENT_SITES
  """

