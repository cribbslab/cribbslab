"""===========================
Pipeline trna
===========================

Overview
========

This pipeline was developed to accurately map small RNA sequencing data and then perform
accurate mapping of tRNA reads and qualitatively analyse the resulting data. trna-frag
has an emphasis on profiling nuclear and mitochondrial tRNA fragments. 


Requires:
 * a single end fastq file - if you have paired end data we recoment flashing the reads together
 to make a single file or only using the first read of your paired end data.
 * a bowtie2 indexed genome
 * ensembl gtf: can be downloaded from 

Pipeline output
===============

The output of running this software is the generation of a html report.

Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGATCore.Pipeline as P
import CGATCore.Experiment as E
import ModuleTrna


# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

###########################################################
# Download the gff of rna types from ucsc
###########################################################

# Connect to the ucsc handle
def connectToUCSC():
    return ModuleTrna.connectToUCSC(
        host=PARAMS["ucsc_host"],
        user=PARAMS["ucsc_user"],
        database=PARAMS["ucsc_database"])


@follows(mkdir("gtf.dir"))
@originate("gtf.dir/rna.gff.gz")
def get_repeat_gff(outfile):
    """This task downloads UCSC repetetive RNA types.
    """
    ModuleTrna.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.as_list(PARAMS["ucsc_rnatypes"]),
        outfile=outfile,
        remove_contigs_regex=PARAMS["ucsc_remove_contigs"],
        job_memory="3G")


##############################################################
# Perform quality control of the fastq files
##############################################################

INPUT_FORMATS = ["*.fastq.gz"]

SEQUENCEFILES_REGEX = r"(\S+).(?P<suffix>fastq.gz)"

@follows(mkdir("fastqc_pre.dir"))
@transform(INPUT_FORMATS,
           suffix(".fastq.gz"),
           r"fastqc_pre.dir/\1.fastq")
def fastqc_pre(infile, outfile):
    """
    Runs fastQC on each input file
    """

    statement = """fastqc -q -o fastqc_pre.dir/ %(infile)s
                """

    P.run(statement)


@follows(mkdir("processed.dir"))
@transform(INPUT_FORMATS,
           suffix(".fastq.gz"),
           r"processed.dir/\1_processed.fastq.gz")
def process_reads(infile, outfile):
    """
    Runs trimmomatic quality related trimming
    """

    if PARAMS["trimmomatic_run"]:

        trimmomatic_options = PARAMS["trimmomatic_options"]

        trimmomatic_options = "ILLUMINACLIP:%s:%s:%s:%s" % (
            PARAMS["trimmomatic_adapter"],
            PARAMS["trimmomatic_mismatches"],
            PARAMS["trimmomatic_p_thresh"],
            PARAMS["trimmomatic_c_thresh"]) + "\t" + trimmomatic_options

        phred = PARAMS["trimmomatic_phred"]

        ModuleTrna.process_trimmomatic(infile, outfile, phred,
                                   trimmomatic_options)
    else:
        statement = "ln %(infile)s %(outfile)s"

        P.run(statement)

@follows(mkdir("fastqc_post.dir"))
@transform(process_reads,
           regex("processed.dir/(\S+)_processed.fastq.gz"),
           r"fastqc_post.dir/\1.fastq")
def fastqc_post(infile, outfile):
    """
    Runs fastQC on each of the processed files
    """

    statement = """fastqc -q -o fastqc_post.dir/ %(infile)s
                """

    P.run(statement)

#####################################################
# Count features over a subset of the data
#####################################################

@follows(mkdir("downsample.dir"))
@transform(process_reads,
           regex("processed.dir/(\S+)_processed.fastq.gz"),
           r"downsample.dir/\1.fastq.gz")
def downsample_fastq(infile, outfile):
    """
    downsamples a fastq file to 500,000 reads each
    """

    statement = """
                seqtk sample -s100 %(infile)s 500000 > %(outfile)s
                """

    P.run(statement)

@follows(mkdir("mapping.dir"))
@transform(downsample_fastq,
           regex("downsample.dir/(\S+).fastq.gz"),
           add_inputs(os.path.join(PARAMS["bowtie2_genome_dir"],
                            PARAMS["bowtie2_genome"] + ".fa")),
           r"mapping.dir/\1.bam")
def map_with_bowtie2(infiles, outfile):
    """
    map reads with bowtie2
    """
    fastq, genome = infiles

    genome = genome.replace(".fa", "")

    samfile = outfile.replace(".bam", ".sam")

    statement = """bowtie2 %(bowtie2_options)s -x %(genome)s -U %(fastq)s -S %(samfile)s 2>%(outfile)s_bowtie.log &&
                   samtools view -bS %(samfile)s |
                   samtools sort -o %(outfile)s &&
                   samtools index %(outfile)s"""

    P.run(statement)


@transform(get_repeat_gff,
           regex("gtf.dir/(\S+).gff.gz"),
           add_inputs(PARAMS['gtf_location']),
           r"gtf.dir/full.gtf")
def process_gtf(infiles, outfile):
    """
    process the gff files so that gene_id is set to source
    so that featurecounts can be ran correctly
    """

    repeats, ensembl = infiles

    statement = """
                zcat %(repeats)s | cgat gff2bed --set-name=class |
                cgat bed2gff --as-gtf | gzip > gtf.dir/rna.gtf.gz &&
                zcat %(gtf_location)s | cgat gff2bed --set-name=source | cgat bed2gff --as-gtf | gzip > gtf.dir/ensembl.gtf.gz &&
                zcat gtf.dir/rna.gtf.gz gtf.dir/ensembl.gtf.gz > %(outfile)s
                """

    P.run(statement)


@follows(mkdir("featurecounts.dir"))
@transform(map_with_bowtie2,
           regex("mapping.dir/(\S+).bam"),
           add_inputs(process_gtf),
           r"featurecounts.dir/\1/gene.tsv")
def count_features(infiles, outfile):
    """
    runs featurecounts to count reads over small RNA features
    """
    
    bamfile, gtf = infiles
    
    name = os.path.basename(bamfile)
    outfolder = name.replace(".bam","")
    intermediate = name.replace(".bam",".tsv")

    statement = """
                featureCounts -t exon -g gene_id -a %(gtf)s -o featurecounts.dir/%(outfolder)s/%(intermediate)s %(bamfile)s &&
                cut -f 1,7 featurecounts.dir/%(intermediate)s > featurecounts.dir/%(outfolder)s/gene.tsv 
                """

    P.run(statement)

################################################
# Perform mapping of tRNA's as set out in Hoffmann et al 2018
################################################

@follows(mkdir("tRNA-mapping.dir"))
@originate("tRNA-mapping.dir/tRNAscan.nuc.csv")
def trna_scan_nuc(outfile):
    """Scans genome using tRNAscanSE to identify nuclear tRNA"""

    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

    statement = "tRNAscan-SE -q -E  %(genome)s 2> tRNA-mapping.dir/tRNAscan.nuc.log | sed 1,3d > %(outfile)s"

# Need to modify if working with non eukaryotic organisms in pipeline.yml- -E to -U
    job_memory = "50G"

    P.run(statement)


@transform(trna_scan_nuc,
           regex("tRNA-mapping.dir/(\S+).nuc.csv"),
           r"tRNA-mapping.dir/\1.bed12")
def trna_scan_mito(infile, outfile):
    """Scans genome using tRNAscanSE to identify mitochrondrial tRNA then outputs
    a bed file of that."""

    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

    tmp_genome = P.get_temp_filename(".")

    statement = """
                cat %(genome)s | perl -lane 'BEGIN{$c=0;}if(m/^>chrM$/){$c=1}elsif(m/^>/){$c=0;}print if $c' > %(tmp_genome)s &&
                tRNAscan-SE -q -O  %(tmp_genome)s | sed 1,3d > tRNA-mapping.dir/tRNAscan.chrM.csv &&
                grep -v chrM %(infile)s > tRNA-mapping.dir/tRNAscan.nuc_mod.csv &&
                cat tRNA-mapping.dir/tRNAscan.nuc_mod.csv tRNA-mapping.dir/tRNAscan.chrM.csv > tRNA-mapping.dir/tRNAscan.csv &&
                perl %(cribbslab)s/perl/tRNAscan2bed12.pl tRNA-mapping.dir/tRNAscan.csv tRNA-mapping.dir/tRNAscan.bed12
                """
    # add onversion for csv to bed file
    P.run(statement)
    os.unlink(tmp_genome)

@transform(os.path.join(PARAMS["genome_dir"],
                            PARAMS["genome"] + ".fa"),
           regex("\S+/([\w_]+[\d_]+).fa"),
           add_inputs(trna_scan_mito),
           r"\1_masked.fa")
def mask_trna_genomic(infiles, outfile):
    """use sam tools to mask fasta ing bedtools """

    genome, bedfile = infiles
    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

    statement = """bedtools maskfasta -fi %(genome)s -fo %(outfile)s -mc N -bed %(bedfile)s"""

    P.run(statement)


@transform(mask_trna_genomic,
           suffix(".fa"),
           add_inputs(trna_scan_mito),
           "_artificial.fa")
def create_pre_trna(infiles, outfile):
    """create pre-tRNA library and then index genome and build segemehl indexes"""

    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")
    genome_name = PARAMS['genome']

    masked_genome, bedfile = infiles

    bedfile = bedfile.replace(".bed12", "")

    statement = """
                perl %(cribbslab)s/perl/modBed12.pl %(bedfile)s.bed12 %(bedfile)s_pre-tRNAs.bed12 &&
                bedtools getfasta -name -split -s -fi %(genome)s -bed %(bedfile)s_pre-tRNAs.bed12 -fo %(genome_name)s_pre-tRNAs.fa &&
                cat %(masked_genome)s %(genome_name)s_pre-tRNAs.fa > %(genome_name)s_artificial.fa &&
                samtools faidx %(genome_name)s_artificial.fa &&
                segemehl.x -x %(genome_name)s_artificial.idx -d %(genome_name)s_artificial.fa 2> segemehl.log
                """

    job_memory = "80G"
    P.run(statement)

@transform(os.path.join(PARAMS["genome_dir"],
                            PARAMS["genome"] + ".fa"),
           regex("\S+/([\w_]+[\d_]+).fa"),
           add_inputs(trna_scan_mito),
           r"\1.fa")
def create_mature_trna(infiles,outfile):
    """will create a library of mature tRNAs
    - remove introns and make fasta from bed12"""

    masked_genome, bedfile = infiles

    statement = """bedtools getfasta -name -split -s -fi %(masked_genome)s -bed %(bedfile)s -fo %(outfile)s"""

    P.run(statement)

@transform(create_mature_trna,
           suffix(".fa"),
           "_mature.fa")
def add_cca_tail(infile, outfile):
    """add CCA tail to the RNA chromosomes and remove pseudogenes"""

    statement = """perl %(cribbslab)s/perl/addCCA.pl %(infile)s %(outfile)s"""

    P.run(statement)

@transform(add_cca_tail,
           suffix("_mature.fa"),
           "_cluster.fa")
def mature_trna_cluster(infile, outfile):
    """mature tRNA clustering - only identical tRNAs are clustered"""

    cluster_info = outfile.replace("_cluster.fa","_clusterInfo.fa")

    statement = "perl %(cribbslab)s/perl/clustering.pl %(infile)s %(outfile)s %(cluster_info)s"

    P.run(statement)

@transform(mature_trna_cluster,
           suffix(".fa"),
           "")
def index_trna_cluster(infile, outfile):
    """index tRNA clusters"""

    cluster_index = infile.replace("fa","idx")
    cluster_dict = infile.replace("fa","dict")

    job_memory = "4G"
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    statement = """samtools faidx %(infile)s &&
                   segemehl.x -x %(cluster_index)s -d %(infile)s 2> segemehl_cluster.log &&
                   picard %(picard_opts)s CreateSequenceDictionary R=%(infile)s O=%(cluster_dict)s"""


    P.run(statement)


@transform()
def pre_mapping_artificial(infile, outfile):
    """pre-mapping of reads against the artificial genome"""

    statement = """ """

    P.run(statement)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
