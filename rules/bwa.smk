pe_samples, = glob_wildcards(os.path.join(config['pe_reads'], '{sample_name}_1.fq'))
reference = config['reference']

rule bwa_index:
  input: config['reference'],
  output: multiext(config['reference'], '.amb', '.ann', '.bwt', '.pac', '.sa')
  conda: 'envs/bwa.yaml'
  shell: 'bwa index {input}'

# usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
# http://www.htslib.org/workflow/
rule bwa_mem:
  input:
    reference=rules.bwa_index.input,
    bwa_index=rules.bwa_index.output,
    fq=os.path.join(config['pe_reads'], '{sample_name}_1.fq'),
    fq2=os.path.join(config['pe_reads'], '{sample_name}_2.fq')
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
  output:
    vcf = 'gatk4/{sample_name}.vcf',
    idx = 'gatk4/{sample_name}.vcf.idx'
  conda: "envs/gatk4.yaml"
  shell: """
  gatk \
    HaplotypeCaller \
    -R {input.reference} \
    -I {input.alignment} \
    -O {output.vcf} \
    --output-mode EMIT_ALL_CONFIDENT_SITES
  """

