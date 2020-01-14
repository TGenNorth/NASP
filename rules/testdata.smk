rule testdata:
  """
  Create a testdata/ directory with randomly generated data to test the workflow.

  Example usage:
    snakemake --use-conda testdata
    snakemake --use-conda -j -d ./testdata/ iqtree
  """
  conda: "envs/bbmap.yaml"
  output:
    'testdata/reference.fasta',
    expand('testdata/ref/genome/1/{ref}', ref=['chr1.chrom.gz', 'info.txt', 'scaffolds.txt.gz', 'summary.txt']),
    expand('testdata/pe_reads/{sample_name}_{read_number}.fq', sample_name=list('abcdefg'), read_number=[1,2]),
    expand('testdata/se_reads/{sample_name}.fq', sample_name=list('abcdefg')),
    expand('testdata/assemblies/{sample_name}.fasta', sample_name=list('abcdefg'))
    # TODO: generate unaligned BAM
    #expand('testdata/ubam/{sample_name}.{ext}', sample_name=list('abcdefgh'), ext=['cram', 'crai'])
  shell: """
    mkdir -pv testdata/{{assemblies,pe_reads,se_reads,ubam}}
    pushd testdata
    randomgenome.sh len=1000 out=reference.fasta
    for sample_name in a b c d e f g; do
      mutate.sh \
        id=0.9 \
        in=reference.fasta \
        out=assemblies/${{sample_name}}.fasta
      randomreads.sh \
        paired=t \
        snprate=0.5 \
        ref=reference.fasta \
        out1=pe_reads/${{sample_name}}_1.fq \
        out2=pe_reads/${{sample_name}}_2.fq \
        len=100 \
        coverage=15
      randomreads.sh \
        paired=f \
        snprate=0.5 \
        ref=reference.fasta \
        out1=se_reads/${{sample_name}}.fq \
        len=100 \
        coverage=15
    done
    popd
  """

