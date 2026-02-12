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
import importlib.util
import cribbslab


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
        pipelines = []
        for path in paths:
            pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of available pipelines are:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))
        return

    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)

    # remove 'cribbslab' from sys.argv
    del sys.argv[0]

    # Find the pipeline module in the paths
    module_file = None
    for search_path in paths:
        candidate = os.path.join(search_path, pipeline + ".py")
        if os.path.exists(candidate):
            module_file = candidate
            break

    if module_file is None:
        raise ImportError(f"Cannot find pipeline: {pipeline}")

    # Load the module using importlib
    spec = importlib.util.spec_from_file_location(pipeline, module_file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[pipeline] = module
    spec.loader.exec_module(module)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
