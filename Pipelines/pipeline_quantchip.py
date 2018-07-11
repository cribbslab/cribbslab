"""===========================
Pipeline template
===========================

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline was written to convert a bamfile from a quant_chip experiment into a homer tag file that can
be processed and normalised using SP1 spike ins. 

Required a filtered bam file

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import zipfile
import sqlite3
import CGATCore.Experiment as E
import CGAT.BamTools.bamtools as Bamtools
import CGATCore.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelinePeakcalling as PipelinePeakcalling

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


SEQUENCEFILES = tuple([os.path.join(DATADIR, "*.bam")])
BEDFILES = tuple([os.path.join("Bed.dir", "*.bed")])
SEQUENCEFILES_REGEX = regex(
    r"(\S+).(bam)")


###################################################################
###################################################################
###################################################################
# load number of reads
###################################################################


@transform(SEQUENCEFILES,
           suffix(".bam"),
           ".idxstats")
def getIdxstats(infiles, outfile):
    '''gets idxstats for bam file so number of reads per chromosome can
    be plotted later'''
    # I have had to add a sleep to make sure the output is written before
    # the next test.
    statement = '''samtools idxstats %(infiles)s > %(outfile)s && sleep 20'''
    P.run(statement)



@follows(mkdir("bedGraph.dir"))
@transform(SEQUENCEFILES,
           suffix(".bam"),
           add_inputs([getIdxstats]),
           r"bedGraph.dir/\1.bedgraph.gz")
def buildBedGraph(infile, outfile):
    '''build wiggle files from bam files.
    Generate :term:`bigWig` format file from :term:`bam` alignment file
    Parameters
    ----------
    infile : str
       Input filename in :term:`bam` format
    outfile : str
       Output filename in :term:`bigwig` format
    annotations_interface_contigs : str
       :term:`PARAMS`
       Input filename in :term:`bed` format
    '''
    inf = infile[0]
    inf_name = inf.replace(".bam", "")
    idxstats = infile[1]

    # scale by Million reads mapped
    reads_mapped = Bamtools.getNumberOfAlignments(inf)

    for idx in idxstats:
        file_name = idx.replace(".idxstats", "")
        if file_name == inf_name:
            # pass to a function that extracts the number of reads aligned to
            # spike in and human genome
            regex = PARAMS['quant_regex'] + "*"
            scale = PipelinePeakcalling.getSpikeInReads(idx, str(regex))
            contig_sizes = PipelinePeakcalling.getContigSizes(idx)
        else:
            continue

    tmpfile = P.get_temp_filename()
    tmpfile2 = P.get_temp_filename()
    job_memory = "3G"
    statement = '''bedtools genomecov
    -ibam %(inf)s
    -g %(contig_sizes)s
    -bg
    -scale %(scale)f
    > %(tmpfile)s &&
    sort -k1,1 -k2,2n -o %(tmpfile2)s %(tmpfile)s &&
    cat %(tmpfile2)s | grep chr | gzip > %(outfile)s &&
    rm -f %(tmpfile)s %(tmpfile2)s
    '''
    P.run(statement)

@transform(buildBedGraph, regex("bedGraph.dir/(\S+).bedgraph.gz"),
           r"\1/\1.txt")
def makeTagDirectory(infile, outfile):
    '''
    This will create a tag file for each bam file
    for a CHIP-seq experiment
    '''

    bamstrip = infile.replace(".bedgraph.gz", "")
    file_name = bamstrip.replace("bedGraph.dir/","")

    statement = '''
                   makeTagDirectory %(file_name)s -force5th
                   -precision 3 -format bed %(infile)s &> %(file_name)s.makeTagInput.log &&
                   touch %(file_name)s/%(file_name)s.txt'''

    P.run(statement)

@follows(mkdir("Hist.dir"))
@transform(makeTagDirectory,
           regex("\S+/(\S+).txt"),
           r"Hist.dir/\1.txt")
def annotatePeaks(infile, outfile):

    '''
    This will create a histogram txt output of the data centred on tss
    '''
    infile = os.path.basename(infile)
    infile = infile.replace(".txt", "")

    statement = '''
                annotatePeaks.pl tss hg19 -d %(infile)s -hist 50 > %(outfile)s
                '''

    P.run(statement)

@transform(makeTagDirectory,
           regex("\S+/(\S+).txt"),
           add_inputs([BEDFILES]),
           r"Hist.dir/\1_total.txt")
def annotatePeaksBed(infiles, outfile):

    '''
    This will create a histogram txt output of the data centred on tss
    '''
    infile = infiles[0]
    infile = os.path.basename(infile)
    infile = infile.replace(".txt", "")


    bedfile = infiles[1][0][0]
    statement = '''
                annotatePeaks.pl %(bedfile)s hg19 -d %(infile)s -size 4000 -hist 50 > %(outfile)s
                '''

    P.run(statement)
# ---------------------------------------------------
# Generic pipeline tasks
@follows(getIdxstats,buildBedGraph, makeTagDirectory,
         annotatePeaks, annotatePeaksBed)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
