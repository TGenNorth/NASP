.. Internal hyperlink target
.. _vcf_to_matrix:

vcf_to_matrix
-------------

::

    vcf_to_matrix --mode=xml --num-threads=1 --dto-file DTO_FILE

.. TODO: Is there a delimiter between the --input-files flag list of INPUT_FILES?

Options:

.. foo -h, --help  show this help message and exit.

--mode mode  Data passing mode
       MODE  can be either 'commandline' or 'xml'.
--reference-fasta file  Path to input reference fasta file.
--reference-dups file  Path to input reference dups file.
--input-files files  Path to input VCF/fasta files for matrix conversion.
              INPUT_FILES is a list of file paths
--matrix-folder path  Name of folder to write output matries to.
--stats-folder path  Name of folder to write statistics files to.
--minimum-coverage number  Minimum coverage depth at a position.
--minimum-proportion number  Minimum proportion of reads that must match the call at a position.
--num-threads number  Number of threads to use when processing input.
--dto-file file  Path to a matrix_dto XML file that defines all the parameters.
