"""===========================
Mapping pipeline for bulkRNAseq data
author: Adam Cribbs
===========================
To run locally: python pipeline_featurecounts.py make full --local

Overview:
    - This pipeline uses the hisat2 mapping tool to align reads
      for bulkRNAseq data analysis and then featurecounts to count
      reads to features

Input files:
    - Paired end BulkRNA-seq fastq.gz files labeled sampleA.fastq.1.gz and sampleA.fastq.2.gz
      where possitions 1 & 2 represent read 1 and 2 from paired end sequencing.
    - Zipped fasta (.gtf.gz) reference genome transfer file file for relevant model organism.

Output files:
    - multiqc html reports (before and after pseudoalignment)
    - merged counts matrix

Processing:
    - fastqc for reads
    - mapping transcripts using hisat2
    - counting features using featurecounts 
    - hisat2 QC

Code
====
"""

# import packages necessary to run the script
import sys
import os
from ruffus import *
import cgatcore.pipeline as P

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])



SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    )

SEQUENCEFILES = tuple([os.path.join('.', suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|fastq.gz|fa.gz)")




@follows(mkdir("fastqc"))                       # make a directory called fastqc, and only proceed once this directory is made
@transform(SEQUENCEFILES,                      # (*) - match any character any number of times [global expression, glob]
            SEQUENCEFILES_REGEX,       # (.*) - match any character (.) any number of times (*) [regular expression, regex]
            r"fastqc/fastqc.html")  # r - raw string to python: fastqc/ - create the following file names in fastqc dir. \1 matches the first (.*) in regex
def fastqc(infile, outfile):
    ''' run fastqc on all fastq.gz files
    input:
        string representing gzipped fastq file name e.g. sampleA.fastq.1.gz
    output:
        string representing fastqc html report name for fastq file e.g. fastqc/sampleA_fastqc.1.html
    example_statment:
        fastqc --nogroup -o fastqc sampleA.fastq.1.gz > fastqc/sampleA.fastq.1.gz.log
     '''
    statement = "fastqc --nogroup -o fastqc %(infile)s > %(outfile)s.log"
    P.run(statement)


@follows(fastqc, mkdir("reports"))               # make a directory called reports, and only proceed once this directory is made
@merge(fastqc, "reports/multiqc_report.html")    # take the output of the fastqc task as input and merge it into a single output named multiqc_report.html in reports dir
def multiqc(infiles, outfile):
    ''' run multiqc to collect fastq stats
    input:
        string of all files output from fastqc
    output:
        string of multiqc file html report
        reports/fastqc_report.html
    example statement:
        export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&
        multiqc fastqc/ -f -n reports/fastqc_report.html
    '''
    # Export statements come from error of mutliqc on server and this makes it run
    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -f -d -s -n %(outfile)s""" # -n - specifies the name of the output file, -f - forces to overwrite the report (if it was run before) so the multiqc outfile is always updated
    P.run(statement, job_memory="30G")


@follows (mkdir ("hisat_mapping"))
@transform(SEQUENCEFILES,                      # (*) - match any character any number of times [global expression, glob]
            SEQUENCEFILES_REGEX,
           r'hisat_mapping/\1.bam')
def hisat2_map(infile, outfile):
    ''' run hisat2 index
    input:
        zipped reference trascriptome fasta file (.fa.gz)
    output:
        .idx file
    example statement:
        kallisto index -i %(outfile)s %(infile)s
    '''

    job_threads = PARAMS["hisat2_threads"]
    job_memory = PARAMS["hisat2_memory"]
    
    index_prefix = PARAMS['hisat2_index_dir'] + "/" + PARAMS['genome']


    if PARAMS['hisat2_paired']:

        infile2 = infile.replace('.1.gz','.2.gz')

        statement = '''hisat2 %%(hisat2_options)s -x %(index_prefix)s --threads %%(job_threads)i -1 %(infile)s
                       -2 %(infile2)s  2> %(outfile)s.log | samtools view -bS > %(outfile)s  2>> %(outfile)s.log''' % locals()
    else:

        statement = '''hisat2 %%(hisat2_options)s -x %(index_prefix)s --threads %%(job_threads)i -U %(infile)s 2> %(outfile)s.log | samtools view -bS > %(outfile)s  2>> %(outfile)s.log''' % locals()

    P.run(statement)


@follows(mkdir("featurecounts"))
@transform(hisat2_map,
           regex("hisat_mapping/(\S+).bam"),
           r"featurecounts/\1.features.tsv")
def count_features(infile, outfile):
    """
    runs featurecounts to count reads over features
    """

    gtf = PARAMS['featurecounts_gtf']


    name = os.path.basename(infile)
    intermediate = name.replace(".bam",".tsv.tmp")

    statement = """
                featureCounts -t exon -g gene_id -a %(gtf)s -o featurecounts.dir/%(intermediate)s %(infile)s &&
                cut -f 1,7 featurecounts.dir/%(intermediate)s > %(outfile)s
               """

    P.run(statement)


@collate(count_features,
         regex("featurecounts/(\S+).features.tsv"),
         r"featurecounts.tsv.gz")
def merge_features(infiles, outfile):
    """This function will merge all of the outputs from featurecounts and
    create a single tsv file for all samples"""

    def merge_feature_data(infiles):
        '''will merge all of the input files '''

        final_df = pd.DataFrame()
        for infile in infiles:
            tmp_df = pd.read_table(infile, sep="\t", index_col=0, skiprows=1)
            final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True)
        final_df = final_df.rename(columns=lambda x: re.sub(".bam","",x))
        final_df = final_df.rename(columns=lambda x: re.sub("hisat_mapping/","",x))

        return final_df

    features = merge_feature_data(infiles)

    features.to_csv(outfile, sep="\t", header=True, index=True, compression='gzip')


@follows(hisat2_map)
@merge(hisat2_map, "hisat2_multiqc.html")
def mapping_multiqc(infiles, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc -f -n hisat2_multiqc.html -o .'''
    P.run(statement)


@follows(fastqc, merge_features, mapping_multiqc)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
    
# this main function lets it work from the command line - tells cgatcore pipeline to run the code
if __name__ == '__main__':
    sys.exit( P.main(sys.argv) )
