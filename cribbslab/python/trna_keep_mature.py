"""trna_generate_bed.py - generate a bed file from a fasta file of pre-tRNAs
================================================================

Purpose
-------

This script takes as an input a fasta file of pre-tRNAs and will
generate a bed file using the name of each fasta read as the
chromosome name and will output the coordinates of the mature
tRNAs as a bed file.


Usage
-----


Options
-------

-m, --merge-pairs
    Output one region per fragment rather than one region per read,
    thus a single region is create stretching from the start of the
    frist read in pair to the end of the second.

    Read pairs that meet the following criteria are removed:

    * Reads where one of the pair is unmapped
    * Reads that are not paired
    * Reads where the pairs are mapped to different chromosomes
    * Reads where the the insert size is not between the max and
      min (see below)

Type::

   python trna_keep_mature.py --help

for command line help.

Command line options
--------------------

"""

import sys
import re
import cgat.FastaIterator as FastaIterator
import cgatcore.iotools as IOTools
import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    (options, args) = E.start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    E.info(options.stdin)
    fastafile = IOTools.open_file(options.stdin.name)

    fasta = FastaIterator.FastaIterator(fastafile)

    for line in fasta:
        chrom = line.title
        total_len = len(line.sequence)
    
        trna_list = []
        string = None
    
        n = 0
        for letter in line.sequence:
        
            n += 1
            if n == 1:
                string = letter
            else:     
                if string.isupper() and letter.isupper():
                    string = str(string) + str(letter)
                elif string.isupper() and letter.islower():
                    trna_list.append(string)
                    string = letter
                elif string.islower() and letter.islower():
                    string = str(string) + str(letter)
                elif string.islower() and letter.isupper():
                    trna_list.append(string)
                    string = letter
        trna_list.append(string)

        start = 1
        end = 1
        chrom = line.title
        for sequence in trna_list:
            start = end
            end = start + len(sequence)
        
            if sequence.islower():
                strand = chrom.split("(")[1].split(")")[0]
                options.stdout.write(("%s\t%s\t%s\t%s\t%s\t%s\n")%(chrom, start, end, chrom, ".", strand))



    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
