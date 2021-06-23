"""===========================
Pseudoalignment pipeline for bulkRNAseq data
author: Alina Kurjan
===========================
To run locally: python pipeline_pseudobulk.py make full -v 5 --local

Overview:
    - This pipeline uses pseudoalignment tool Kallisto to align reads
      for bulkRNAseq data analysis

Input files:
    - Paired end BulkRNA-seq fastq.gz files labeled sampleA.fastq.1.gz and sampleA.fastq.2.gz
      where possitions 1 & 2 represent read 1 and 2 from paired end sequencing.
    - Zipped fasta (.fa.gz) reference transcriptome file for relevant model organism
      (download from https://www.ensembl.org/info/data/ftp/index.html).
    - Chromosome size and GTF files for relevant model organism (if needing gene counts)
      (human g38 download from: https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/)

Output files:
    - multiqc html reports (before and after pseudoalignment)
    - abundance.tsv
    - counts_matrix.csv with mapped genes as rows and samples as columns

Processing:
    - fastqc
    - mapping transcripts using Kallisto
    - translating transcripts into genes

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


@follows(mkdir("fastqc"))                   # make a directory called fastqc, and only proceed once this directory is made
@transform("*.fastq.*.gz",                  # (*) - match any character any number of times [global expression, glob]
            regex(r"(.*).fastq.(.*).gz"),   # (.*) - match any character (.) any number of times (*) [regular expression, regex]
            r"fastqc/\1.fastq.\2_fastqc.html")    # r - raw string to python: fastqc/ - create the following file names in fastqc dir. \1 matches the first (.*) in regex
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


@follows(mkdir("reports"))                      # make a directory called reports, and only proceed once this directory is made
@merge(fastqc, "reports/multiqc_report.html")    # take the output of the fastqc task as input and merge it into a single output named fastqc_report.html in reports dir
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
    P.run(statement)

# Restructure fasta reference transcriptome file into a kallisto-friendly index
@follows (mkdir ("kalindex"))
@transform("*.fa.gz", regex(r'(.*).fa.gz'), r'kalindex/\1.idx')
def kallisto_index(infile, outfile):
    ''' run multiqc to transform fasta format to kallisto format
    input:
        zipped reference trascriptome fasta file (.fa.gz)
    output:
        .idx file
    example statement:
        kallisto index -i %(outfile)s %(infile)s
    '''
    statement = '''kallisto index -i %(outfile)s %(infile)s'''
    P.run(statement)

@follows (mkdir ("quant"))
@transform('*.fastq.1.gz', regex(r'(\S+).fastq.1.gz'), add_inputs(kallisto_index), r'quant/\1/abundance.tsv')
def kal_quant (infiles, outfile):
    ''' pseudoaligment command to align reads to transcripts
        input:
            1) zipped, paired fastq files. Infiles are in the format:
            infiles=(('SampleA.fastq.1.gz', 'kalindex/Reference_transcriptome.idx'),
            ('SampleA.fastq.2.gz', 'kalindex/eference_transcriptome.idx')).
            Need to index to extract them properly!
            2) kallisto index file.
            3) GTF file (describes where the transcripts lie in the genome).
            4) text file containing the length of each chromosome.
        output:
            abundance.tsv file with estimated counts of transcripts per million (TPM)
            abundance.tsv.log file with run log
            run_info.json file with info about the run get_parameters
            abundance.h5 file with main quantification and bootstraps
            .BAM file with gene counts
    '''
    infile, index_file = infiles
    infile1 = infile
    infile2 = infile.replace(".1.gz",".2.gz")
    output_folder = P.snip(outfile, "/abundance.tsv")
    # need to make it so that only one sample (file pair) is supplied at a time!!!
    statement = '''kallisto quant %(kal_quant_options)s
                   -t %(kal_quant_threads)s
                   -b %(kal_quant_bootstraps)s
                   -i %(index_file)s
                   -o %(output_folder)s %(infile1)s %(infile2)s > %(outfile)s.log 2>&1'''
                   #--genomebam --gtf Homo_sapiens.GRCh38.104.gtf.gz --chromosomes hg38.chrom.sizes.txt - not working yet
    P.run(statement, job_threads = PARAMS["kal_quant_threads"])

@merge (kal_quant, "quant/kallisto_multiqc.html")
def kallisto_multiqc(infiles, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc -f -n kallisto_multiqc.html -o quant quant'''
    P.run(statement)

@follows(fastqc, multiqc, kallisto_index, kal_quant, kallisto_multiqc)
def full():
    pass

# this main function lets it work from the command line - tells cgatcore pipeline to run the code
if __name__ == '__main__':
    sys.exit( P.main(sys.argv) )


#https://github.com/cgat-developers/cgat-flow/blob/4e8741781a7a4ca08c8e507bab3f56be078f6e89/cgatpipelines/tools/pipeline_rnaseqdiffexpression.py#L559
#feature counts

#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto
