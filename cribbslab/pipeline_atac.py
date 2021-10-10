"""===========================
ATAC pipeline for bulk ATACseq data
author: Adam Cribbs
===========================
To run locally: python pipeline_atac.py make full --local

Overview:
    - This pipeline uses Bowtie2 to align reads
      for bulkATACseq data analysis

Input files:
    - Paired end BulkATAC-seq fastq.gz files labeled sampleA.fastq.1.gz and sampleA.fastq.2.gz
      where possitions 1 & 2 represent read 1 and 2 from paired end sequencing.
    - A reference genome.

Output files:
    - A mapped bam file
    - 

Processing:
    - fastqc
    - Mapping using Bowtie2
    - Bowtie2 QC
    - 

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


@follows(mkdir("fastqc"))                       # make a directory called fastqc, and only proceed once this directory is made
@transform("*.fastq.*.gz",                      # (*) - match any character any number of times [global expression, glob]
            regex(r"(.*).fastq.(.*).gz"),       # (.*) - match any character (.) any number of times (*) [regular expression, regex]
            r"fastqc/\1.fastq.\2_fastqc.html")  # r - raw string to python: fastqc/ - create the following file names in fastqc dir. \1 matches the first (.*) in regex
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

# Restructure fasta reference transcriptome file into a kallisto-friendly index
@follows (mkdir ("bowtie2_index"))
@transform(PARAMS['genome_fasta'], regex('(.*).fa'), r'bowtie2_index/\1.idx')
def bowtie2_index(infile, outfile):
    ''' Build a bowtie2 index
    input:
        reference genome fasta file (.fa)
    output:
        .bt2 reference
    example statement:
        kallisto index -i %(outfile)s %(infile)s
    '''
    name = infile.replace(".fa", "")
    name = "bowtie2_index" + name
    
    statement = '''bowtie2-build %(outfile)s %(name)s'''
    
    P.run(statement)

@follows (mkdir ("Bowtie2"))
@transform('*.fastq.1.gz', regex(r'(\S+).fastq.1.gz'), add_inputs(kallisto_index), r'Bowtie2/\1.bam')
def bowtie2_map (infiles, outfile):
        ''' Mapping command to align reads to genome and produce bam file.
        input:
            zipped, paired fastq files.
            bowtie2 index file.
        output:
            bowtie2 mapped bam file
        '''

        infile, index_file = infiles
        infile1 = infile
        infile2 = infile.replace(".1.gz",".2.gz")
        output_folder = P.snip(outfile, "/abundance.tsv")
        statement = '''kallisto quant %(kal_quant_options)s
                        -t %(kal_quant_threads)s
                        -b %(kal_quant_bootstraps)s
                        -i %(index_file)s
                        -o %(output_folder)s %(infile1)s %(infile2)s > %(outfile)s.log 2>&1'''
        P.run(statement, job_threads = PARAMS["kal_quant_threads"])


@merge(bowtie2_map, "mapping_qc/Bowtie2_multiqc.html")
def bowtie2_multiqc(infiles, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc -f -n Bowtie2_multiqc.html -o mapping_qc'''
    P.run(statement)

@follows(bowtie2_multiqc)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
    
# this main function lets it work from the command line - tells cgatcore pipeline to run the code
if __name__ == '__main__':
    sys.exit( P.main(sys.argv) )
