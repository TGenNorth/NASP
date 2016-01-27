Installation
============

Install
-------

.. note:: NASP requires Python 3. If your default version is Python 2, you should be able to use the Python 3 version explicitly by replacing the below 'pip' commands with 'pip3'

If you have virtualenvwrapper installed::

	$ mkvirtualenv nasp
	$ pip install nasp

If you have venv installed::

	$ pyvenv my_env
	$ source my_env/bin/activate
	$ pip install nasp

At the command line::

	$ easy_install nasp

Using the “module load” system:
If you have configured a module system on your high-performance compute cluster, you can create a NASP module for your users to load. Refer to your operating system documentation for specifics. Do not forget to include the software dependencies, as described below.

Update
------

::

	$ pip install --upgrade nasp

Uninstall
---------

::

	$ pip uninstall nasp

Dependencies
------------

- Samtools_ < 1.3
- trimmomatic_
- MUMmer_ >= 3.23

Optional aligners for fastq files
---------------------------------

- BWA_ >= 0.7.5a
- Novocraft_ >= 3.02.04
- Bowtie2_ >= 2.1.0
- Samtools_ < 1.3

Optional snpcallers for bam/fastq files:
----------------------------------------

- VarScan_ >= 2.3.6
- SolSNP_
- GATK_
- picard-tools_ required with GATK_

See the snpcaller documentation for which version of Java_ is required

.. _Samtools: http://samtools.sourceforge.net/
.. _trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic 
.. _MUMmer: http://mummer.sourceforge.net/
.. _BWA: http://bio-bwa.sourceforge.net/
.. _Novocraft: http://www.novocraft.com/main/page.php?s=novoalign
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _VarScan: http://varscan.sourceforge.net/
.. _SolSNP: http://sourceforge.net/projects/solsnp/
.. _GATK: https://www.broadinstitute.org/gatk/
.. _picard-tools: https://broadinstitute.github.io/picard/
.. _Java: http://www.java.com/en/
