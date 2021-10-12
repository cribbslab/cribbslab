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


@active_if(PARAMS['bowtie2_index_true'])
@follows (mkdir ("bowtie2_index"))
@transform(PARAMS['genome_fasta'], regex('(.*).fa'), r'bowtie2_index/\1.1.bt2')
def bowtie2_index(infile, outfile):
    ''' Build a bowtie2 index
    input:
        reference genome fasta file (.fa)
    output:
        .bt2 reference
    example statement:
        kallisto index -i %(outfile)s %(infile)s
    '''
    name = os.path.basename(infile.replace(".fa", ""))
    name = "bowtie2_index/" + name
    
    statement = '''bowtie2-build %(infile)s %(name)s'''
    
    P.run(statement)

@follows (mkdir ("Bowtie2"))
@transform('*.fastq.1.gz', regex(r'(\S+).fastq.1.gz'), add_inputs(bowtie2_index), r'Bowtie2/\1.bam')
def bowtie2_map (infiles, outfile):
        ''' Mapping command to align reads to genome and produce bam file.
        input:
            zipped, paired fastq files.
            bowtie2 index file.
        output:
            bowtie2 mapped bam file
        '''

        if not PARAMS['bowtie2_index_true']:
            infile = infiles
        else:
            infile, index_file = infiles
        
        infile1 = "".join(infile)
        infile2 = infile1.replace(".1.gz",".2.gz")
        
        if not PARAMS['bowtie2_index_true']:
            index_file = PARAMS['bowtie2_index_path']

        name = index_file.replace(".1.bt2", "")
        bamname = outfile.replace(".bam", "")

        statement = '''bowtie2 --very-sensitive -k 3 -x %(name)s -1 %(infile1)s -2 %(infile2)s > %(bamname)s.bowtie.sam 2> %(bamname)s.log &&
                       samtools view -S -b %(bamname)s.bowtie.sam > %(bamname)s.bowtie.bam &&
                       samtools sort %(bamname)s.bowtie.bam -o %(outfile)s &&
                       samtools index %(outfile)s &&
                       rm -rf %(bamname)s.bowtie.sam &&
                       rm -rf %(bamname)s.bowtie.bam'''

        P.run(statement,
              job_threads = PARAMS["bowtie2_threads"],
              job_memory = PARAMS['bowtie2_memory'])



@merge(bowtie2_map, "reports/Bowtie2_multiqc.html")
def bowtie2_multiqc(infiles, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc Bowtie2/ -f -n %(outfile)s'''
    P.run(statement)


@follows(mkdir("HMMRATAC"))
@transform(bowtie2_map, regex(r'Bowtie2/(\S+).bam'), r'HMMRATAC/\1.genome.info')
def make_genomeinfo(infile, outfile):
    '''Make a genome.info file for running HMMRATAC'''

    statement = '''samtools view -H %(infile)s| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\\t",$2,"\\n"}' > %(outfile)s'''

    P.run(statement)


###################################################
### Peak callers
###################################################


@transform(make_genomeinfo, regex(r'HMMRATAC/(\S+).genome.info'), r'HMMRATAC/\1_peaks.gappedPeak')
def run_hmmratac(infile, outfile):
    '''Run HMMRATAC '''

    bamfile = infile.replace(".genome.info", ".bam")
    bamfile = bamfile.replace("HMMRATAC/", "Bowtie2/")

    outname = outfile.replace("_peaks.gappedPeak", "")

    statement = '''HMMRATAC -b %(bamfile)s -i %(bamfile)s.bai -g %(infile)s -o %(outname)s --window 1000000'''

    P.run(statement, job_memory="50G")


@transform(run_hmmratac, regex(r'HMMRATAC/(\S+)_peaks.gappedPeak'), r'HMMRATAC/\1_filterpeaks.gappedPeak')
def filter_hmmratac(infile, outfile):


    statement = '''awk -v OFS="\t" '$13>=10 {print}' %(infile)s > %(outfile)s '''

    P.run(statement)


@transform(run_hmmratac, regex(r'HMMRATAC/(\S+)_peaks.gappedPeak'), r'HMMRATAC/\1_filteredSummits.bed')
def filtersummit_hmmratac(infile, outfile):

    infile = infile.replace("_peaks.gappedPeak", "_summits.bed")

    statement = '''awk -v OFS="\t" '$5>=10 {print}' %(infile)s > %(outfile)s'''

    P.run(statement)


@follows(multiqc, bowtie2_multiqc, filter_hmmratac, filtersummit_hmmratac)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
    
# this main function lets it work from the command line - tells cgatcore pipeline to run the code
if __name__ == '__main__':
    sys.exit( P.main(sys.argv) )
