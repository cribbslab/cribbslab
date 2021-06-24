'''
cribbslab - workflows for cribbslab
========================================================

Our cribbslab workflows can be ran using the following commands:

For this message and a list of available keywords type::

    cribbslab --help


To see the available pipelines within each section type::

    cribbslab <section>

To run a specific pipeline/workflow type the following::

    cribbslab <section> <workflow> [workflow options] [workflow arguments]

To get help for a specify workflow, type::

    cribbslab <section> <workflow> --help
'''

import os
import sys
import re
import glob
import imp
import scpipelines


def printListInColumns(l, ncolumns):
    '''output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % 3 != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):

    argv = sys.argv

    # paths to look for pipelines:
    print(cribbslab.__file__)
    path = os.path.abspath(os.path.dirname(cribbslab.__file__))


    paths = [path]

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":

        print((globals()["__doc__"]))
        print("The list of available sections are:\n")
        return

    elif argv[1] == "main":
        print((globals()["__doc__"]))

        pipelines = []
        pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print("The list of available pipelines is:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))

    try:
        command = argv[2]
        pipeline = "pipeline_{}".format(command)
    except:
        print("No pipeline has been selected under the %s section" % (argv[1]))
        return

    # remove 'scflow' from sys.argv
    del sys.argv[0]

    (file, pathname, description) = imp.find_module(pipeline, paths)

    module = imp.load_module(pipeline, file, pathname, description)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())