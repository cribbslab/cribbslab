"""add_cca_tail.py - Adds CCA tails to fasta file sequences
================================================================

Purpose
-------
This script adds CCA tails to the RNA chromosomes and remove pseudogenes. It takes fasta files as input and outputs fasta files.

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
import re
import cgat.FastaIterator as FastaIterator
import cgatcore.iotools as IOTools
import cgatcore.experiment as E
import collections


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


    infile = IOTools.open_file(options.stdin.name)
    iterator = FastaIterator.FastaIterator(infile)

   # outfile_info = IOTools.open_file(options.info_file, "w")

    d = collections.OrderedDict()
    cluster_dict = dict()

    # first iterate over the fasta file and generate a dict
    # with the name (title) as the key and the sequence as the value
    # Remove any pseudo sequences
    for cur_record in iterator:

        key = cur_record.title
        if "pseudo" in key:
            pass

        else:
            d[key] = cur_record.sequence

    # next iterate of over the dict give the cluster a number
    # this will be used to then map back for the info name

    for key, value in d.items():       
        # Add CCA tail
        options.stdout.write((">%s\n%scca\n")%(key, value))

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
