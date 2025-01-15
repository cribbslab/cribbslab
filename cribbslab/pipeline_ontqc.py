"""Pipeline for quality control of Oxford Nanopore and PacBio long-read sequencing data.

This pipeline runs quality control tools including NanoPlot and LongQC on long-read
sequencing data to assess read quality, length distributions, and other QC metrics.

Configuration parameters are stored in pipeline.yml.
"""

import sys
import os
import sqlite3
import pandas as pd

from ruffus import *
from cgatcore import pipeline as P
from cgatcore import experiment as E

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Helper functions
def get_fastq_files(infile):
    """Return a list of fastq files from a input file or pattern."""
    if isinstance(infile, str):
        if os.path.isfile(infile):
            return [infile]
        else:
            return glob.glob(infile)
    return infile

# Pipeline tasks
@follows(mkdir("nanoplot"))
@transform("*.fastq.gz",
           regex(r"(\S+).fastq.gz"),
           r"nanoplot/\1/NanoPlot-report.html")
def runNanoPlot(infile, outfile):
    """Run NanoPlot on input fastq files."""
    outdir = os.path.dirname(outfile)
    
    statement = """NanoPlot 
                  --fastq %(infile)s 
                  --outdir %(outdir)s 
                  --threads %(nanoplot_threads)s 
                  --maxlength %(nanoplot_maxlength)s 
                  --minlength %(nanoplot_minlength)s 
                  %(nanoplot_options)s"""
    
    P.run(statement)

@follows(mkdir("fastqc"))
@transform("*.fastq.gz",
           regex(r"(\S+).fastq.gz"),
           r"fastqc/\1_fastqc.html")
def runFastQC(infile, outfile):
    """Run FastQC on input fastq files."""
    outdir = os.path.dirname(outfile)
    
    statement = """export _JAVA_OPTIONS="-Xmx%(fastqc_memory)sG" && fastqc 
                  --outdir %(outdir)s 
                  --threads %(fastqc_threads)s 
                  %(fastqc_options)s 
                  %(infile)s"""
    
    P.run(statement, job_memory=str(PARAMS["fastqc"]["memory"]) + "G", job_options='--time=0-02:00:00')

@follows(mkdir("multiqc"))
@follows(runNanoPlot, runFastQC)
@merge([r"nanoplot/*/NanoPlot-report.html",
        r"fastqc/*.html"],
       "multiqc/multiqc_report.html")
def runMultiQC(infiles, outfile):
    """Run MultiQC to collect all QC metrics."""
    outdir = os.path.dirname(outfile)
    
    statement = """multiqc 
                  --force 
                  --outdir %(outdir)s 
                  nanoplot/ 
                  fastqc/"""
    
    P.run(statement)

@follows(runMultiQC)
def full():
    """Run the full pipeline."""
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
