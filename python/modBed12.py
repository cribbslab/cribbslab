"""modBed12.py - modify a bed file
================================================================

Purpose
-------
Create_pre_trna


Usage
-----


Options
-------


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
  
    lines = infile.readlines()

    for line in lines:
       
        column = line.split()

        new_columns = [column[0], str(int(column[1])-50), str(int(column[2]) +50), column[3], column[4], column[5], str(int(column[1])-50), str(int(column[2]) +50), column[8], column[9]]

        if "pseudo" not in column[3]:
           
            if int(column[9]) == 2:
                [c,d] = column[10].split(",")
                block = int(column[2])-int(column[1]) - int(d) + 50
                new_10 = ''.join(str(int(c)+50) + ',' + str(int(d) + 50))
                new_11 = ''.join('0' + ',' + str(block))
                new_columns = new_columns + [new_10, new_11]

            else:
                new_columns = new_columns + [str(int(column[10]) + 100), column[11]]

            options.stdout.write('\t'.join(new_columns[0:]) + '\n')

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
