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
