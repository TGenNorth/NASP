matrix
------

::

    usage: nasp matrix --dto-file matrix_dto.xml

Matrix aggregates the VCF and Frankenfasta files from the NASP pipeline into .tsv matrix and stats file summaries
which may be parsed for further analysis.

Given the --dto-file flag, all other flags are optional overrides.
A reference fasta and vcf/frankenfasta files are always required, either from
the dto or listed on the command line.

--dto-file           Path to the matrix_dto.xml file
--num-threads        Max number of CPUs that can be executing simultaneously (default: all)
--reference-fasta    Path to the reference.fasta against which samples are compared
--reference-dups     Path to the duplicates.txt file marking duplicated positions
--stats-folder       Path to the output folder for statistics (default: ./)
--matrix-folder      Path to the output folder for matrices (default: ./)
--minimum-coverage   Filter positions below this coverage/depth threshold (default: 0)
--minimum-proportion Filter positions below this proportion threshold (default: 0.0)
--withallrefpos      Include the withallrefpos.tsv matrix
