#!/usr/bin/env python3
#
# Snakemake stores the job properties (e.g. name, rulename, threads, input, output, params etc.) as JSON inside the job
# script (for group jobs, the rulename will be “GROUP”, otherwise it will be the same as the job name).
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# do something useful with the threads
threads = job_properties[threads]

# access property defined in the cluster configuration file (Snakemake >=3.6.0)
job_properties["cluster"]["time"]

os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))
