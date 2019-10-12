"""===========================
Pipeline bam2fastq
===========================

Overview
========

This pipeline parallelises the conversion of bam files to
fastq files

Input files
-----------

bam files are required to be located within the directory that
the pipeline is executed.

Pipeline output
===============

The output is a pair of fastq files.


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


SEQUENCESUFFIXES = ("*.bam"
                    )

BAMTARGET = tuple([os.path.join(".", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

@subdivide(BAMTARGET,
         regex("(\S+).bam"),
         [r"\1.fastq.1.gz", r"\1.fastq.2.gz"])
def bam2fastq_paired(infile, outfiles):
    '''Convert a bam file to two fastq files.'''

    outf1, outf2 = outfiles
    out1 = outf1.replace(".gz", "")
    out2 = outf2.replace(".gz", "")

    statement = '''samtools sort %(infile)s -o %(infile)s.tmp &&
                   bedtools bamtofastq -i %(infile)s.tmp -fq %(out1)s -fq2 %(out2)s &&
                   gzip %(out1)s &&
                   gzip %(out2)s &&
                   rm -rf %(infile)s.tmp'''

    P.run(statement)

@transform(BAMTARGET,
           regex("(\S+).bam"),
           r"\1.fastq.gz")
def bam2fastq_single(infile, outfile):
    '''Convert a bam file to one fastq files.'''

    out1 = infile.replace(".gz","")

    statement = '''
                   bedtools bamtofastq -i %(infile)s.tmp -fq %(out1)s &&
                   gzip %(out1)s &&
                   rm -rf %(infile)s.tmp'''

    P.run(statement)

@follows(bam2fastq_paired)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
