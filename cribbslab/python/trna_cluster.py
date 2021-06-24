"""trna_cluster.py - generate a bed file from a fasta file of pre-tRNAs
================================================================

Purpose
-------
Mature tRNA clustering - only identical tRNAs are clustered.
This script takes a fasta file of mature tRNAs
with CCA tails as input and outputs a fasta file of identical sequences clustered
together with a cluster number and chromsomal coordinates for one of the sequences.

Task: mature_trna_cluster
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

    parser.add_option(
        "--info-file-out", dest="info_file", type="str",
        help="name of the info file name[default=%default]")

    parser.set_defaults(
        info_file="info_file.fa"
    )

    (options, args) = E.start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    E.info(options.stdin)


    infile = IOTools.open_file(options.stdin.name)
    iterator = FastaIterator.FastaIterator(infile)

    outfile_info = IOTools.open_file(options.info_file, "w")

    d = collections.OrderedDict()
    cluster_dict = dict()

    # first iterate over the fasta file and generate a dict
    # with the sequnce as the key and the name as the value
    # only add if the sequence occurs once
    for cur_record in iterator:

        key = cur_record.sequence
        if key in d:
            pass
        else:
            d[key] = cur_record.title
    # next iterate of over the dict give the cluster a number
    # this will be used to then map back for the info name
    n = 0
    for key, value in d.items():
        n +=1
        cluster_dict[key] = n
        # output this to std out

        m = re.match("(\S+)-(\S+)(\(\S+\))", value)

        value = m.group(2)

        options.stdout.write((">cluster%s:%s\n%s\n")%(n, value, key))

    # iterate over the infile again, this time use the 
    # sequence to pull out the cluster it belongs to
    
    infile = IOTools.open_file(options.stdin.name)
    iterator = FastaIterator.FastaIterator(infile)

    for cur_record in iterator:
        cluster = cluster_dict[cur_record.sequence]
        outfile_info.write((">cluster%s:%s\n%s\n")%(cluster, cur_record.title, cur_record.sequence))


    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
