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
import CGAT.FastaIterator as FastaIterator
import CGATCore.IOTools as IOTools
import CGATCore.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    parser.add_option("-f", "--fasta-file", dest="fasta_file",
                      action="store_true",
                      help="merge paired-ended reads and output interval "
                      "for entire fragment [default=%default]. ")

    parser.set_defaults(
    )

    (options, args) = E.start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")


    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
