Installation
============

Install
-------

.. note:: NASP requires Python 3.


To install system-wide::

    $ pip install nasp

.. _python documentation: https://packaging.python.org/en/latest/installing/

One of the ways to install NASP as a user is in a python virtual environment. See the `python documentation`_ for alternative methods::

    $ pyvenv my_env
    $ source my_env/bin/activate
    $ pip install nasp

If Python 2 is your default version, you should be able to run the Python 3 versions of pyvenv/pip explicitly::

    $ pyvenv-3.4 my_env
    $ source my_env/bin/activate
    $ pip install nasp

The readline module is an optional dependency to enable tab-autocomplete when entering filepaths in NASP::
    
    $ pip install readline

To learn more about installing python packages checkout the `python documentation`_

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
