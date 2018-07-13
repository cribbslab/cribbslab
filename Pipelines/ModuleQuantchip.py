"""
This module file contains functions that help with running pipeline_quantchip.py

"""

import pandas as pd

def read_design_table(infile):
    """
    This function reads a design table named 'design.tsv' and generates
    objects to be used to match peaks in samples to the appropriate control,
    treatment and input samples.

    This will return a dictionary linking each file to the specific samplt type.
    """
    # read the design file

    df = pd.read_csv(infile, sep="\t")

    return(df)
    
def extract_bam_files(dataframe, control):
    """
    This function takes in a dataframe of the design.tsv and then
    extracts the corresponding Treatment and Input from the Control
    name.
    """
    df = pd.read_csv(dataframe, sep="\t")

    df = df.set_index('Control')
    treatment = df.loc[control]["Treatment"]
    inputD = df.loc[control]["Input"]

    control = control.replace(".bam", ".txt")
    treatment = treatment.replace(".bam", ".txt")
    inputD = inputD.replace(".bam", ".txt")

    return(control, treatment, inputD)

def getSpikeInReads(idx, regex):

    idx_file = pd.read_csv(idx, "\t",
                           header=None,
                           names=["Name", "Length", "Mapped", "Unmapped"])

    spikes = idx_file[idx_file.Name.str.match(regex)]
    spikes = spikes.Mapped.sum()

    mapped_reads = idx_file[-idx_file.Name.str.match(regex)]
    mapped_reads = mapped_reads.Mapped.sum()

    scale = spikes/mapped_reads

    return scale

def getContigSizes(idx):
    idx_file = pd.read_csv(idx, "\t",
                           header=None,
                           names=["Name", "Length", "Mapped", "Unmapped"])

    contigs = idx_file[['Name', 'Length']]

    contig_file = idx.replace(".idxstats", ".contig")

    contigs.to_csv(contig_file, sep="\t", header=None, index=None)

    return contig_file
