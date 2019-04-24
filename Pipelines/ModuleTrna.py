"""
ModuleTrna.py - Tasks for running trna pipleine

"""

import os
import re
import pysam
import pandas as pd
import seaborn as sns
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import cgatcore.pipeline as P
import cgatcore.database as Database
import cgat.FastaIterator as FastaIterator
import cgat.GTF as GTF


def merge_feature_data(infiles):
    '''will merge all of the input files '''

    final_df = pd.DataFrame()
    for infile in infiles:
        tmp_df = pd.read_table(infile, sep="\t", index_col=0, skiprows=1)
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True)
    final_df = final_df.rename(columns=lambda x: re.sub(".bam","",x))
    final_df = final_df.rename(columns=lambda x: re.sub("mapping.dir/","",x))

    return final_df

def merge_counts_data(infiles):
    '''will merge all of the output of idx stats into a counts
    table so it can be parsed into r for differential analysis'''

    final_df = pd.DataFrame()
    for infile in infiles:
        name = infile.replace(".idxstats", "")
        name = name.replace("post_mapping_bams.dir/","")
        tmp_df = pd.read_table(infile, sep="\t", header=None, names=["length",name,"unmapped"],index_col=0)
        tmp_df = tmp_df[name]
        final_df = final_df.merge(tmp_df.to_frame(), how="outer", left_index=True, right_index=True)
        
    return final_df

def getNumReadsFromReadsFile(infile):
    '''get number of reads from a .nreads file.'''
    with IOTools.open_file(infile) as inf:
        line = inf.readline()
        if not line.startswith("nreads"):
            raise ValueError(
                "parsing error in file '%s': "
                "expected first line to start with 'nreads'")
        nreads = line.split("\t")[1]
        nreads = int(nreads)
    return nreads

def connectToUCSC(host="genome-mysql.cse.ucsc.edu",
                  user="genome",
                  database=None):
    """connect to UCSC database.

    Arguments
    ---------
    host : string
        Host to connect to
    user : string
        Username to connect with
    Database : string
        database to use

    Returns
    -------
    Database handle

    """
    dbhandle = Database.connect(url="mysql://{user}@{host}/{database}".format(**locals()))

    return dbhandle




def getRepeatDataFromUCSC(dbhandle,
                          repclasses,
                          outfile,
                          remove_contigs_regex=None,
                          job_memory="4G"):
    '''download data from UCSC database and write to `outfile` in
    :term:`gff` format.

    This method downloads repeats from the repeatmasker track at
    the UCSC.

    Arguments
    ---------
    dbhandle : object
       Database handle to UCSC mysql database
    repclasses : list
       List of repeat classes to select. If empty, all repeat classes
       will be collected.
    outfile : string
       Filename of output file in :term:`gff` format.
    remove_contigs_regex : string
       If given, remove repeats on contigs matching the regular
       expression given.

    '''
    cc = dbhandle.execute("SHOW TABLES LIKE '%%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    # now collect repeats
    tmpfile = P.get_temp_file(".")

    for table in tables:

        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd,
        '.', strand, '.',
        CONCAT('class \\"', repClass, '\\"; family \\"',
        repFamily, '\\"; repName \\"', repName, '\\";')
        FROM %(table)s"""

        if repclasses:
            repclasses_str = ",".join(
                ["'" + x.strip() + "'" for x in repclasses])
            sql += ''' WHERE repClass in (%(repclasses_str)s) ''' % locals()

        sql = sql % locals()

        E.debug("executing sql statement: %s" % sql)
        cc = dbhandle.execute(sql)
        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")

    tmpfile.close()

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = ['''cat %(tmpfilename)s
    | sort -t$'\\t' -k1,1 -k4,4n
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log ''']

    if remove_contigs_regex:
        statement.append('--contig-pattern="{}"'.format(
            ",".join(remove_contigs_regex)))

    statement.append('''| gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run(statement, job_memory=job_memory)

    os.unlink(tmpfilename)


def process_trimmomatic(infile, outfile, phred, trimmomatic_options):
    """
    Runs the trimmomatic software
    """

    output_prefix = "processed.dir/" + infile.replace(".fastq.gz", "")
    job_threads = 2
    job_memory = "12G"

    statement = """
                trimmomatic SE -threads %(job_threads)s %(phred)s %(infile)s %(outfile)s
                %(trimmomatic_options)s 2>> %(output_prefix)s.log
                """ %locals()

    P.run(statement)

def coverage(idx, coverage, outfile):
    """
    merges idx data with the coverage data and then outputs a data frame
    that can then be plotted using ggplots in r
    """

    coverage = pd.read_csv(coverage, skiprows=3, sep="\t")
    idx = pd.read_csv(idx, names=["#CHROM", "length", "mapped", "unmapped"], sep="\t")

    total = np.sum(idx.loc[:,'mapped':].values)

    idx['percent'] = idx.loc[:,'mapped':].sum(axis=1)/total * 100

    idx = idx.sort_values(by=['percent'], ascending=False).head(50)

    # merge two tables
    df = idx.merge(coverage, left_on='#CHROM', right_on='#CHROM')

    # split post_mapping column into two
    df['GT'], df['read_count'] = df['post_mapping'].str.split(':', 1).str
