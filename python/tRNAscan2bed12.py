""" tRNAscan2bed12.py - convert indexed csv files into bed12 format
======================================================================

Purpose
-------

This script takes csv files of indexed nuclear and mitocondrial tRNA combines
and converts them to bed12 format.
tRNAscan2bed12.pl (Hoffmann et al.) converted to python file, using CGAT tools

Usage
-----

Options
-------

*
*
*

Type:: 

python tRNAscan2bed12.py --help

for command line help

"""

import sys
import re
import cgatcore.iotools as IOTools
import cgatcore.experiment as E
import collections
import pandas as pd
import numpy as np


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    parser.add_option(
        "--info-file-out", dest="info_file", type="str",
        help="name of the info file name[default=%default]")

    (options, args) = E.start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    E.info(options.stdin)


    infile = IOTools.open_file(options.stdin.name)


    # Using pandas module to convert csv into data frame
    df = pd.read_csv(infile, header = None, sep='\t')

    row_number = len(df.index)

    # Loop over lines of csv file
    for i in range(0,row_number):

        # Chromosome number/name    
        chr_num = df.iloc[i,0].rstrip()

        # Score
        score = df.iloc[i,8]*10
        if score > 1000:
            score = 1000 

        # Unique ID
        u_id = str(chr_num +'.tRNA' +str(int(df.iloc[i,1])) + '-' + str(df.iloc[i,4]) +str( df.iloc[i,5]) + '-' + str(df.iloc[i,9]))
        # Last part of unique ID either pseudo or nan, replace nan with space
        u_id = u_id.replace("nan", "")
        block_start = ''
        block_size = ''
    
        # Set strand direction and left/right 
        if df.iloc[i,3] > df.iloc[i,2]:

            left = df.iloc[i,2] -1

            right = df.iloc[i,3]

            strand = "+"

        else:
            left = df.iloc[i,3] -1

            right = df.iloc[i,2]

            strand = "-"
        
        if df.iloc[i,6] and df.iloc[i,7]:

            block_count = 2
        
            if strand == "+":

                block_start = block_start.join('0'+','+str(df.iloc[i,7]-df.iloc[i,2]+1))

                block_size = block_size.join(str(df.iloc[i,6]-df.iloc[i,2])+','+str(df.iloc[i,3]-df.iloc[i,7]))
            
            elif strand == "-":

                block_start = block_start.join('0'+','+str(df.iloc[i,6]-df.iloc[i,3]+1))

                block_size = block_size.join(str(df.iloc[i,7]-df.iloc[i,3])+','+str(df.iloc[i,2]-df.iloc[i,6]))
    
        else:
            block_count = 1

            block_size  = abs(df.iloc[i,3]-df.iloc[i,2])+1

            block_start = 0
    
        # Write line by line to output file
        l = [chr_num, str(left), str(right), u_id, str(score), strand, str(left), str(right), "0", str(block_count), str(block_size), str(block_start)]
        
        options.stdout.write('\t'.join(l[0:]) + '\n')
     
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
