"""===========================
Pipeline quantchip
===========================

Overview
========

This pipeline was written to convert a bamfile from a quant_chip experiment into a homer tag file that can
be processed and normalised using spike ins specified as a regex in the yml file. 

Requires:
 * filtered bam file and an associated index
 * A bed file containing regions of interest (i.e. TSS regions). This file should be placed in
   a folder called Bed.dir/
 * a design file named design.tsv which names the control, treatment and input bam files

Pipeline output
===============

The output is a series of annotaed peak output files from homer that can be plotted using an R script

Code
====

"""
from ruffus import *

import sys
import os
import zipfile
import sqlite3
import cgatcore.experiment as E
import cgat.BamTools.bamtools as Bamtools
import cgatcore.pipeline as P
import cribbslab.ModuleQuantchip as ModuleQuantchip

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


##################################################################
# Load the design file into pipeline 
##################################################################


### IN FUTURE ALL REFERENCES TO PIPELINE PEAKCALLING NEED TO BE REMOVED

if os.path.exists("design.tsv"):
    df = ModuleQuantchip.read_design_table("design.tsv")
    
    CONTROLBAMS = list(set(df["Control"].values))
    TREATMENTBAMS = list(set(df["Treatment"].values))
    INPUTBAMS = list(set(df["Input"].values))
else:
    E.warn("design.tsv is not present within the folder")




SEQUENCEFILES = tuple([os.path.join(".", "*.bam")])
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
    statement = '''samtools index %(infiles)s && samtools idxstats %(infiles)s > %(outfile)s && sleep 20'''
    P.run(statement)



@follows(mkdir("bedGraph.dir"))
@transform(SEQUENCEFILES,
           suffix(".bam"),
           add_inputs([getIdxstats]),
           r"bedGraph.dir/\1.bedgraph")
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
            scale = ModuleQuantchip.getSpikeInReads(idx, str(regex))
            contig_sizes = ModuleQuantchip.getContigSizes(idx)
        else:
            continue

    tmpfile = P.get_temp_filename()
    tmpfile2 = P.get_temp_filename()
    job_memory = "30G"
    statement = '''bedtools genomecov
    -ibam %(inf)s
    -g %(contig_sizes)s
    -bg
    -scale %(scale)f
    > %(tmpfile)s &&
    sort -k1,1 -k2,2n -o %(tmpfile2)s %(tmpfile)s &&
    cat %(tmpfile2)s | grep chr  > %(outfile)s &&
    rm -f %(tmpfile)s %(tmpfile2)s
    '''
    P.run(statement)


@transform(buildBedGraph,
           regex("bedGraph.dir/(\S+).bedgraph"),
           r"\1/\1.txt")
def makeTagDirectory(infile, outfile):
    '''
    This will create a tag file for each bam file
    for a CHIP-seq experiment
    '''

    bamstrip = infile.replace(".bedgraph", "")
    file_name = bamstrip.replace("bedGraph.dir/./","")

    statement = '''
                   makeTagDirectory %(file_name)s -force5th
                   -precision 3 -format bed %(infile)s &> %(file_name)s.makeTagInput.log &&
                   touch %(file_name)s/%(file_name)s.txt'''

    P.run(statement)

@transform(buildBedGraph,
           regex("bedGraph.dir/(\S+).bedgraph"),
           add_inputs([BEDFILES]),
           r"Hist.dir/\1.txt")
def annotatePeaksBed(infiles, outfile):

    '''
    This will create a histogram txt output of the data centred on tss
    '''
    bedGraph = infiles[0]

    try:
        bedfile = infiles[1][0][0]
    except IndexError as error:
        print(error)
        print("No bed file detected")

    statement = '''
                annotatePeaks.pl %(bedfile)s hg19 -bedGraph %(bedGraph)s %(homer_options)s  > %(outfile)s
                '''

    P.run(statement)

@follows(annotatePeaksBed)
@transform(CONTROLBAMS,
           regex("(\S+).bam"),
           r"Hist.dir/\1.eps")
def plot_hist(infile, outfile):
    """
    This function will plot a histogram using the Coverageplot.R
    script
    """

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           "R"))
    control = infile
    dataframe = "design.tsv"

    control, treatment, inputD = ModuleQuantchip.extract_bam_files(dataframe, control)

    statement = """ Rscript %(R_SRC_PATH)s/CoveragePlot.R
                --control=Hist.dir/%(control)s
                --treatment=Hist.dir/%(treatment)s
                --input=Hist.dir/%(inputD)s
                """

    P.run(statement)


@follows(mkdir("annotate.dir"))
@transform(makeTagDirectory,
           regex("(\S+)/(\S+).txt"),
           r"\1/regions.txt")
def find_peaks(infile, outfile):
    """This function will peak call for histone peaks"""

    tag_dir = os.path.dirname(infile)

    statement = """
                findPeaks %(tag_dir)s -region -size 150 -minDist 370 -o auto
                """
    job_memory = "10G"
    P.run(statement)


@follows(mkdir("annotate.dir"))
@transform(buildBedGraph,
           regex("bedGraph.dir/(\S+).bedgraph"),
           r"annotate.dir/\1.xls")
def call_peaks_macs3(infile, outfile):
    """Will call peaks using MACS2 from a bedgraph"""

    out_dir = os.path.dirname(outfile)
    name = os.path.basename(outfile).replace(".xls","")

    outname = outfile.replace("annotate.dir/", "")
    outname = outname.replace(".xls", "")

    job_memory="50G"
    statement = """
                macs3 bdgpeakcall -i %(infile)s  --outdir annotate.dir/ >& annotate.dir/%(name)s.log -o %(outname)s
                """

    P.run(statement, job_condaenv="macs3", job_memory="50G")

#####################################################
# Create bigWig for plotting to IGV
#####################################################

@transform(buildBedGraph,
           regex("(\S+).bedgraph"),
           r"\1.bw")
def bedgraph_to_bw(infile, outfile):
    """This function will generate a bigwig from a bedgraph"""    
    
    contig = os.path.basename(infile)
    contig = contig.replace("bedgraph", "contig")
    
    statement = """bedGraphToBigWig %(infile)s %(contig)s %(outfile)s"""

    P.run(statement)


if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    BW_REGEX = regex(
        r"bedGraph.dir/./%s.bw" % (
            PARAMS["merge_pattern_input"].strip()))

    # the last expression counts number of groups in pattern_input
    BW_OUTPUT = r"mergeBedGraph.dir/%s.bw" % (
            PARAMS["merge_pattern_output"].strip())

    @follows(mkdir("mergeBedGraph.dir"))
    @collate(bedgraph_to_bw,
             BW_REGEX,
             BW_OUTPUT)
    def merge_bw(infiles, outfile):
        """Merge bigWigs using mergeBigWig"""    


        infiles = " ".join(infiles)

        tmpfile = P.get_temp_filename()
        tmpfile2 = P.get_temp_filename()

        statement = '''bigWigMerge %(infiles)s %(tmpfile)s &&
                       LC_COLLATE-C sort -k1,1 -k2,2n -o %(tmpfile2)s %(tmpfile)s &&
                       bedGraphToBigWig %(tmpfile2)s %(contig_file)s %(outfile)s'''

        P.run(statement)
else:
    def merge_bw():
        pass


# ---------------------------------------------------
# Generic pipeline tasks
@follows(getIdxstats,buildBedGraph, makeTagDirectory,
         annotatePeaksBed, plot_hist, bedgraph_to_bw, merge_bw, call_peaks_macs3)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
