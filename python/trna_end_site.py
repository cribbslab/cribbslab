"""trna_end_site.py - calculates and plots the end of each read
================================================================

Purpose
-------
This script takes in a fasta of the clustered tRNA and calculates the length.
It also takes as an input a bam file and then identified the end of each read,
making a histogram of all of the end sites for each tRNA-species.

Usage
-----


Options
-------
**


Type::


for command line help.

Command line options
--------------------

"""

import sys
import os
import re
import pysam
import pandas as pd
import seaborn as sns
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import cgat.FastaIterator as FastaIterator



def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    parser.add_option("-d", "--outputdir",
                      dest="outdir", type="string",
                      help="output directory to save plots")

    parser.add_option("-f", "--fasta",
                      dest="fasta_file", type="string",
                      help="fasta file containing tRNA cluster fasta seqs")


    parser.set_defaults(fasta_file =None,
                        outdir=None)

    (options, args) = E.start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    E.info(options.stdin)


    dict_trna = {}
    for record in FastaIterator.iterate(IOTools.open_file(options.fasta_file)):
        title = record.title.strip("-")
        length = len(record.sequence)
        dict_trna[title] = length
 
   # For each read in bamfile find end position and then plot this using length of tRNA cluster
    samfile = pysam.AlignmentFile(options.stdin.name, "rb")
    refname = ""
    values = []
    n = 0
    for line in samfile:
        if line.reference_name == refname:
            if line.reference_end is None:
                pass
            else:
                end = int(line.reference_end) - int(line.reference_start)
                values.append(end)
        elif line.reference_name != refname:
            n += 1
            if n > 1:

                values = pd.Series(values)
                percent = values.value_counts() / values.count() * 100
                percent = percent.sort_index()
                percent = pd.DataFrame(percent)
                percent.rename(columns={0: 'Percent'}, inplace=True)
            
                # length of each tRNA from fasta
                length = dict_trna[refname.strip("-")] + 1
            
                temp_df = pd.DataFrame(0, index=range(1,length), columns=['A'])
                temp_df = pd.concat([temp_df, percent], axis=1)
                percent = temp_df.fillna(0)

                refname = options.outdir + refname.strip("-")
                outfile = refname + ".csv"
                outfig = refname + ".png"

                percent.to_csv(outfile)
            
                g = sns.factorplot(x=percent.index, y="Percent", data=percent,
                                   size=8, kind="bar", palette="Blues")
                g.set_xlabels('position from 5\' end')
                g.set_xticklabels(rotation=90)
                g.savefig(outfig)
            
                values = []
                refname = line.reference_name


            else:
            
                refname = line.reference_name

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
