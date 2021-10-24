"""===========================
Pipeline ATAC
===========================

Overview
========

This pipeline performs various ATAC seq downstream analysis. The pipeline
requires a mapped bam file and macs2 output files.

Usage
=====


Configuration
-------------

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_gat.py config

Input files
-----------

* pipeline.yml config file
* design.csv file to configure the analysis
* bam files
* macs2 output files

Requirements
------------

Pipeline output
===============

Output is a number of RData objects that can be passed into  Rmarkdown
for plotting.


Code
====

"""

import sys
import os
from ruffus import *
import cgatcore.pipeline as P
import cgatcore.experiment as E

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


@transform("design.csv",
           regex("(design).csv"),
           r"DiffBind.Rdata")
def run_diffbind(infile, outfile):
    '''Run diffbind and generate Rdata output'''

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           "R"))

    statement = '''Rscript %(R_SRC_PATH)s/DiffBind.R --design %(infile)s --output %(outfile)s'''

    P.run(statement, job_memory="100G")

@follows(run_diffbind)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))