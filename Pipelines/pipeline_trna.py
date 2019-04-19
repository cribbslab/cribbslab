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
 * a bowtie indexed genome
 * ensembl gtf: can be downloaded from 

segemehl is quite a slow mapper in comparrison to others. However it improves the quality of the tRNA alignment

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
import pandas as pd
import cgatcore.pipeline as P
import cgatcore.experiment as E
import ModuleTrna
import cgat.IndexedFasta as IndexedFasta


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

    statement = "fastqc -q -o fastqc_pre.dir/ %(infile)s"

    P.run(statement)

@follows(fastqc_pre)
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

        statement = "cp %(infile)s %(outfile)s"

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
                seqtk sample -s100 %(infile)s 500000 | gzip > %(outfile)s
                """

    P.run(statement)

@follows(mkdir("mapping.dir"))
@transform(downsample_fastq,
           regex("downsample.dir/(\S+).fastq.gz"),
           add_inputs(os.path.join(PARAMS["bowtie_genome_dir"],
                            PARAMS["bowtie_genome"] + ".fa")),
           r"mapping.dir/\1.bam")
def map_with_bowtie(infiles, outfile):
    """
    map reads with bowtie to get general alignment so features can be counted
    over RNA gene_biotypes
    """
    fastq, genome = infiles
    tmp_fastq = P.get_temp_filename(".")
    temp_file = P.get_temp_filename(".")
    genome = genome.replace(".fa", "")

    statement = """gzip -dc %(fastq)s > %(tmp_fastq)s && bowtie -k 10 -v 2 --best --strata --sam  %(genome)s  %(tmp_fastq)s 2> %(outfile)s_bowtie.log | samtools view -bS |
                   samtools sort -T %(temp_file)s -o %(outfile)s &&
                   samtools index %(outfile)s
                """

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
                zcat %(gtf_location)s | cgat gff2bed --set-name=gene_biotype | cgat bed2gff --as-gtf | gzip > gtf.dir/ensembl.gtf.gz &&
                zcat gtf.dir/rna.gtf.gz gtf.dir/ensembl.gtf.gz > %(outfile)s
                """

    P.run(statement)


@follows(mkdir("featurecounts.dir"))
@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           add_inputs(process_gtf),
           r"featurecounts.dir/\1/\1.feature_small.tsv")
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
                cut -f 1,7 featurecounts.dir/%(outfolder)s/%(intermediate)s > %(outfile)s 
                """

    P.run(statement)

@collate(count_features,
         regex("featurecounts.dir/(\S+)/(\S+).feature_small.tsv"),
         r"featurecounts.dir/merged_features.tsv")
def merge_features(infiles, outfile):
    """This function will merge all of the outputs from featurecounts and
    create a single tsv file for all samples"""


    features = ModuleTrna.merge_feature_data(infiles)
    
    features.to_csv(outfile, sep="\t", header=True, index=True)

###############################################
# Quality statistics for small RNA on genome mapping
###############################################

@follows(mkdir("genome_statistics.dir"))
@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           r"genome_statistics.dir/\1.strand")
def strand_specificity(infile, outfile):
    '''This function will determine the strand specificity of your library
    from the bam file'''

    statement = (
        "cgat bam2libtype "
        "--max-iterations 10000 "
        "< {infile} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


@follows(mkdir("genome_statistics.dir"))
@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           r"genome_statistics.dir/\1.nreads")
def count_reads(infile, outfile):
    '''Count number of reads in input files.'''

    statement = '''printf "nreads \\t" >> %(outfile)s'''

    P.run(statement)

    statement = '''samtools view %(infile)s | wc -l | xargs printf >> %(outfile)s'''

    P.run(statement)

# Add inputs adds all samples as infiles
@follows(mkdir("genome_statistics.dir"))
@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           add_inputs(count_reads, get_repeat_gff),
           r"genome_statistics.dir/\1.readstats")
def build_bam_stats(infiles, outfile):
    '''count number of reads mapped, duplicates, etc.
    Excludes regions overlapping repetitive RNA sequences
    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename with number of reads per sample
    outfile : str
       Output filename with read stats
    annotations_interface_rna_gtf : str
        :term:`PARMS`. :term:`gtf` format file with repetitive rna
    '''

    job_memory = "32G"
 
    # Only one sample 
    if len(infiles)==3:
        bamfile, readsfile, rna_file = infiles
    # If there are multiple samples, programme specifies which .nreads file to use, by matching name to bam file
    else:
        bamfile = infiles[0]
        rna_file = infiles[-1]
        # Split file name up into directory and file name(/), then further split up by file name and file type and take file name (.)
        bam_name = bamfile.split('/')[1].split('.')[0]
        for i in range(1,len(infiles)-1):
            nread_name = infiles[i].split('/')[1].split('.')[0]
            if bam_name ==  nread_name:
                readsfile = infiles[i]
                break
            else:
                continue    

    nreads = ModuleTrna.getNumReadsFromReadsFile(readsfile)
    track = P.snip(os.path.basename(readsfile),
                   ".nreads")

    # if a fastq file exists, submit for counting
    if os.path.exists(track + ".fastq.gz"):
        fastqfile = track + ".fastq.gz"
    elif os.path.exists(track + ".fastq.1.gz"):
        fastqfile = track + ".fastq.1.gz"
    else:
        fastqfile = None

    if fastqfile is not None:
        fastq_option = "--fastq-file=%s" % fastqfile
    else:
        fastq_option = ""

    statement = '''
    cgat bam2stats
         %(fastq_option)s
         --force-output
         --mask-bed-file=%(rna_file)s
         --ignore-masked-reads
         --num-reads=%(nreads)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(bamfile)s
    > %(outfile)s
    '''

    P.run(statement)

@follows(mkdir("genome_statistics.dir"))
@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           r"genome_statistics.dir/\1.idxstats")
def full_genome_idxstats(infile, outfile):
    """This will generate idxstats to count the number of mapped
       and unmapped reads per contig"""

    statement = "samtools idxstats %(infile)s > %(outfile)s"

    P.run(statement)

@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           r"genome_statistics.dir/\1.stats")
def build_samtools_stats(infile, outfile):
    '''gets stats for bam file so number of reads per chromosome can
    be plotted later'''

    statement = '''samtools stats %(infile)s > %(outfile)s'''

    P.run(statement)


@transform(map_with_bowtie,
           regex("mapping.dir/(\S+).bam"),
           add_inputs(os.path.join(PARAMS["bowtie_genome_dir"],
                            PARAMS["bowtie_genome"] + ".fa")),
           r"genome_statistics.dir/\1.genomecov")
def genome_coverage(infiles, outfile):
    """runs bedtoools genomecov to look at the coverage over all
       samples """

    infile, genome = infiles

    statement = """bedtools genomecov -ibam %(infile)s -g %(genome)s > %(outfile)s"""

# Maybe should use hg38_mature.fa instead, would  add input from add_cca_tail
# should use -d to look at every position
    P.run(statement)

################################################
# Perform mapping of tRNA's as set out in Hoffmann et al 2018
################################################

@active_if(PARAMS['trna_scan_load'])
@follows(mkdir("tRNAscan_load.dir"))
@originate("tRNA-mapping.dir/tRNAscan.nuc.csv")
def trna_scan_load(outfile):
    """ If already downloaded trna nuclear genome from http://gtrnadb.ucsc.edu/index.html """

    trna_folder = os.path.join(PARAMS['trna_scan_path'], PARAMS['trna_scan_folder'] + ".tar.gz")
    genome = PARAMS['trna_scan_folder']
    trna_file = "tRNAscan_load.dir/" + genome  + "-detailed.out"
    tmp_genome = P.get_temp_filename(".")

    statement = """ tar -xzf %(trna_folder)s -C tRNAscan_load.dir && 
    cut -f1-9,16 %(trna_file)s | sed 1,3d > %(tmp_genome)s && 
    awk -F"\\t" '{ $10 = ($10 == "pseudo" ? $10 : "") } 1' OFS=, %(tmp_genome)s | sed 's/,/\\t/g' > %(outfile)s   
    """

    P.run(statement)

#awk '{$10=$11=$12=$13=$14=$15=""; print $0}' %(tmp_genome)s > %(outfile)s
#awk '{gsub("high confidence set","hcs");gsub("secondary filtered","sf");gsub("isotype mismatch","isoMis");print}'

   


## Using if else statements
@follows(mkdir("tRNA-mapping.dir"))
@originate("tRNA-mapping.dir/tRNAscan.nuc.csv")
def trna_scan_nuc(outfile):
    """Scans genome using tRNAscanSE to identify nuclear tRNA"""
    if( PARAMS['index_run']):

        genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

        statement = "tRNAscan-SE -q %(genome)s 2> tRNA-mapping.dir/tRNAscan.nuc.log | sed 1,3d > %(outfile)s"

# Need to modify if working with non eukaryotic organisms in pipeline.yml- -E to -U
        job_memory = "50G"

        P.run(statement)
    else:
        indexed_genome = os.path.join(PARAMS['index_genome_dir'], PARAMS['index_genome'] + ".nuc.csv")
        os.symlink(indexed_genome, 'tRNA-mapping.dir/tRNAscan.nuc.csv')
# softlink to location of nuc.csv file
# Need option if downloaded from database


@transform(trna_scan_nuc,
           regex("tRNA-mapping.dir/(\S+).nuc.csv"),
           r"tRNA-mapping.dir/\1.bed12")
def trna_scan_mito(infile, outfile):
    """Scans genome using tRNAscanSE to identify mitochrondrial tRNA then cat the output of nuclear
       scan outputs a bed file of that."""

    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

    tmp_genome = P.get_temp_filename(".")

# For python script

    statement = """
                cat %(genome)s | perl -lane 'BEGIN{$c=0;}if(m/^>chrM$/){$c=1}elsif(m/^>/){$c=0;}print if $c' > %(tmp_genome)s &&
                tRNAscan-SE -q -O  %(tmp_genome)s |sed 1,3d > tRNA-mapping.dir/tRNAscan.chrM.csv &&
                grep -v chrM %(infile)s  > tRNA-mapping.dir/tRNAscan.nuc_mod.csv &&
                cat tRNA-mapping.dir/tRNAscan.nuc_mod.csv tRNA-mapping.dir/tRNAscan.chrM.csv > tRNA-mapping.dir/tRNAscan.csv &&
                python %(cribbslab)s/python/tRNAscan2bed12.py -I tRNA-mapping.dir/tRNAscan.csv -S %(outfile)s
                """
#  --info-file-out=%(trna_bed_info)s
    # add conversion for csv to bed file
# | cat | tr "\t" "," >
    P.run(statement)
    os.unlink(tmp_genome)

@transform(os.path.join(PARAMS["genome_dir"],
                            PARAMS["genome"] + ".fa"),
           regex("\S+/([\w_]+[\d_]+).fa"),
           add_inputs(trna_scan_mito),
           r"tRNA-mapping.dir/\1_masked.fa")
def mask_trna_genomic(infiles, outfile):
    """use sam tools to mask fasta ing bedtools """

    genome, bedfile = infiles
    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")

    statement = """bedtools maskfasta -fi %(genome)s -fo %(outfile)s -mc N -bed %(bedfile)s"""

    P.run(statement)


@transform(mask_trna_genomic,
           regex("tRNA-mapping.dir/(\S+)_masked.fa"),
           add_inputs(trna_scan_mito),
           r"tRNA-mapping.dir/\1_pre-tRNAs.fa")
def create_pre_trna(infiles, outfile):

    masked_genome, bedfile = infiles
    genome = os.path.join(PARAMS['genome_dir'], PARAMS['genome'] + ".fa")
    genome_name = PARAMS['genome']

    bedfile_name = bedfile.replace(".bed12","")

    statement = """
               python %(cribbslab)s/python/modBed12.py -I %(bedfile)s -S %(bedfile_name)s_pre-tRNAs.bed12 &&
                bedtools getfasta -name -split -s -fi %(genome)s -bed %(bedfile_name)s_pre-tRNAs.bed12 -fo %(outfile)s """

    P.run(statement)

@transform(create_pre_trna,
           regex("tRNA-mapping.dir/(\S+)_pre-tRNAs.fa"),
           add_inputs(mask_trna_genomic),
           r"tRNA-mapping.dir/\1_artificial.fa")
def create_artificial(infiles, outfile):
    """create pre-tRNA library and then index genome and build bowtie indexes"""

    genome_name = PARAMS['genome']

    pre_trna, masked_genome = infiles

    statement = """
                cat %(masked_genome)s tRNA-mapping.dir/%(genome_name)s_pre-tRNAs.fa > %(outfile)s &&
                samtools faidx %(outfile)s
                """

    P.run(statement)

@transform(create_artificial,
           regex("tRNA-mapping.dir/(\S+)_artificial.fa"),
           r"tRNA-mapping.dir/\1_artificial.1.ebwt")
def bowtie_index_artificial(infile, outfile):
    '''generate a bowtie index of the artificial genome 
       ================================================
       ================================================
       Generating a bowtie index can take a while..
       Please be patient, do something else.
       ================================================
       '''

    genome_name = PARAMS['genome']

    statement = """ bowtie-build %(infile)s tRNA-mapping.dir/%(genome_name)s_artificial 2> tRNA-mapping.dir/bowtie-build_artificial.log """

    P.run(statement)


@transform(os.path.join(PARAMS["genome_dir"],
                            PARAMS["genome"] + ".fa"),
           regex("\S+/([\w_]+[\d_]+).fa"),
           add_inputs(trna_scan_mito),
           r"tRNA-mapping.dir/\1.fa")
def create_mature_trna(infiles,outfile):
    """will create a library of mature tRNAs
    - remove introns and make fasta from bed12"""

    masked_genome, bedfile = infiles

    statement = """bedtools getfasta -name -split -s -fi %(masked_genome)s -bed %(bedfile)s -fo %(outfile)s"""

    P.run(statement)

@transform(create_mature_trna,
           regex("tRNA-mapping.dir/(\S+).fa"),
           r"tRNA-mapping.dir/\1_mature.fa")
def add_cca_tail(infile, outfile):
    """add CCA tail to the RNA chromosomes and remove pseudogenes"""

    statement = """python %(cribbslab)s/python/addCCA.py -I %(infile)s -S %(outfile)s"""

    P.run(statement)

@transform(add_cca_tail,
           regex("tRNA-mapping.dir/(\S+)_mature.fa"),
           r"tRNA-mapping.dir/\1_cluster.fa")
def mature_trna_cluster(infile, outfile):
    """mature tRNA clustering - only identical tRNAs are clustered"""

    cluster_info = outfile.replace("_cluster.fa","_clusterInfo.fa")

    statement = """python %(cribbslab)s/python/trna_cluster.py -I %(infile)s -S %(outfile)s --info-file-out=%(cluster_info)s"""

    P.run(statement)

@transform(mature_trna_cluster,
           regex("tRNA-mapping.dir/(\S+).fa"),
           r"tRNA-mapping.dir/\1.1.ebwt")
def index_trna_cluster(infile, outfile):
    """index tRNA clusters"""

    genome_name = PARAMS['genome']

    job_memory = "4G"

    statement = """samtools faidx %(infile)s &&
                   bowtie-build %(infile)s tRNA-mapping.dir/%(genome_name)s_cluster 2> bowtie_cluster.log
                """
    P.run(statement)


@follows(mkdir("pre_mapping_bams.dir"))
@transform(process_reads,
           regex("processed.dir/(\S+)_processed.fastq.gz"),
           add_inputs(bowtie_index_artificial),
           r"pre_mapping_bams.dir/\1.bam")
def pre_mapping_artificial(infiles, outfile):
    """pre-mapping of reads against the artificial genome"""

    fastq, bowtie_index_artificial = infiles

    index_name = bowtie_index_artificial.replace(".1.ebwt", "")
    fastq_name = fastq.replace(".fastq.gz","")
    fastq_name = fastq.replace("processed.dir/","")

    statement = """bowtie -n 3 -k 1 --best -e 800 --sam  %(index_name)s %(fastq)s 2> pre_mapping_bams.dir/%(fastq_name)s.log |
                   samtools view -b -o %(outfile)s
                   """


    job_memory = "20G"
    P.run(statement)


@transform(pre_mapping_artificial,
         regex("pre_mapping_bams.dir/(\S+).bam"),
           add_inputs(create_pre_trna),
         r"pre_mapping_bams.dir/\1_filtered.bam")
def remove_reads(infiles, outfile):
    """remove all of the reads mapping at least once to the genome"""

    infile, pre_trna_genome = infiles

    temp_file = P.get_temp_filename(".")
    temp_file1 = P.get_temp_filename(".")
    
    statement = """samtools view -h %(infile)s> %(temp_file)s && 
                   perl %(cribbslab)s/perl/removeGenomeMapper.pl %(pre_trna_genome)s %(temp_file)s %(temp_file1)s &&
                   samtools view -b %(temp_file1)s > %(outfile)s"""

    job_memory = "50G"
    P.run(statement)
    os.unlink(temp_file)
    os.unlink(temp_file1)

@transform(create_pre_trna,
           regex("tRNA-mapping.dir(\S+)_pre-tRNAs.fa"),
           r"tRNA-mapping.dir/\1_mature.bed")
def create_mature_bed(infile, outfile):
    """remove pre-tRNA regions and form a bed file of the mature tRNAs"""

    statement = """python %(cribbslab)s/python/trna_keep_mature.py -I %(infile)s -S %(outfile)s """

    P.run(statement)

@transform(remove_reads,
         regex("pre_mapping_bams.dir/(\S+)_filtered.bam"),
         add_inputs(create_mature_bed),
         r"tRNA-mapping.dir/\1_filtered.fastq.gz")
def keep_mature_trna(infiles, outfile):
    """remove pre-tRNA reads, keep only mature tRNA reads"""

    samfile, bedfile = infiles
    bedfile = bedfile.replace(".bed12", "")
    

    statement = """bedtools intersect -f 1 -wa -abam %(samfile)s -b %(bedfile)s |
                   cgat bam2fastq %(outfile)s"""

    P.run(statement)

#################################################
# Post processing of mapping data
#################################################

@follows(mkdir("post_mapping_bams.dir"))
@transform(keep_mature_trna,
           regex("tRNA-mapping.dir/(\S+)_filtered.fastq.gz"),
           add_inputs(mature_trna_cluster, index_trna_cluster),
           r"post_mapping_bams.dir/\1_trna.bam")
def post_mapping_cluster(infiles, outfile):
    """post mapping against the cluster tRNAs """

    fastqfile, database, bowtie_index_cluster = infiles

    genome_name = bowtie_index_cluster.replace(".1.ebwt","")

    temp_file = P.get_temp_filename(".")

    statement = """bowtie -n 3 -k 1 --best -e 800 --sam  %(genome_name)s %(fastqfile)s 2> tRNA-mapping.dir/cluster.log | samtools view -bS |
                   samtools sort -T %(temp_file)s -o %(outfile)s &&
                   samtools index %(outfile)s"""

    job_memory = "40G"
    P.run(statement)


@transform(mature_trna_cluster,
           regex("tRNA-mapping.dir/(\S+).fa"),
           r"\1_contig.tsv")
def get_contig_cluster(infile, outfile):
    """This will generate a contig file of the cluster genome"""

    fasta = IndexedFasta.IndexedFasta(infile)
    contigs = []

    for contig, size in fasta.getContigSizes(with_synonyms=False).items():
        contigs.append([contig, size])
    df_contig = pd.DataFrame(contigs, columns=['contigs', 'size'])
    df_contig.sort_values('contigs', inplace=True)
    df_contig.to_csv(outfile, sep="\t", header=False, index=False)
    

@follows(mkdir("variant_calling.dir/"))
@transform(post_mapping_cluster,
           regex("post_mapping_bams.dir/(\S+)_trna.bam"),
           add_inputs(mature_trna_cluster),
           r"variant_calling.dir/\1.var.raw.vcf")
def samtools_pileup(infiles, outfile):
    """use samtools mpileup and bcftools to call variants"""

    infile, cluster_fa = infiles

    statement = """samtools mpileup --no-BAQ --output-tags DP,AD -f %(cluster_fa)s --BCF %(infile)s | 
                   bcftools call --consensus-caller --variants-only --pval-threshold 0.05 -Ov  > %(outfile)s
                   """

    P.run(statement)


@transform(samtools_pileup,
           regex("variant_calling.dir/(\S+).var.raw.vcf"),
           add_inputs(mature_trna_cluster),
           r"variant_calling.dir/\1_variants.vcf")
def samtools_norm_indels(infiles, outfile):
    """use bcftools to normalise for indels"""

    infile, cluster_fa = infiles

    statement = """bcftools norm -Ou -m-any %(infile)s | bcftools norm -Ov --check-ref w -f %(cluster_fa)s > %(outfile)s
                   """

    P.run(statement)


@transform(samtools_norm_indels,
           regex("variant_calling.dir/(\S+)_variants.vcf"),
           r"variant_calling.dir/\1.var.flt.vcf")
def filter_vcf(infile, outfile):
    """filters a raw samtools mpileup vcf by depth"""

    statement = """bcftools view %(infile)s | vcfutils.pl varFilter -D100 > %(outfile)s """


    P.run(statement)


##############################################
# Identify tRNA fragment/full length position
##############################################

@follows(remove_reads)
@transform(post_mapping_cluster,
           suffix(".bam"),
           ".idxstats")
def idx_stats_post(infile, outfile):
    """perform samtools idxstats to determine percent expression of each cluster of tRNA
    can be computed - will eventually be passed into r to perform differential 
    expression analysis and calulate the percent computed"""

    statement = "samtools idxstats %(infile)s > %(outfile)s"

    P.run(statement)

@collate(idx_stats_post,
         regex("post_mapping_bams.dir/(\S+).idxstats"),
         r"merged_idxstats.txt.gz")
def merge_idx_stats(infiles, outfile):


    final_df = ModuleTrna.merge_counts_data(infiles)

    final_df.to_csv(outfile, sep="\t", compression="gzip")

####################################
# Collect information regarding the per base % and create a
# table of sample info for plotting
####################################

@transform(post_mapping_cluster,
           regex("post_mapping_bams.dir/(\S+)_trna.bam"),
           add_inputs(mature_trna_cluster),
           r"post_mapping_bams.dir/\1_pileup.tsv")
def create_coverage(infiles, outfile):
    '''use bam_pileup2tsv to compute the coverage over each base'''

    infile, fasta = infiles

    statement = '''cgat bam_pileup2tsv %(infile)s -m coverage-vcf -f %(fasta)s -L %(outfile)s.log > %(outfile)s '''

    P.run(statement)
# need to merge coverage?
'''
@follows(create_coverage,
         regex("post_mapping_bams.dir/(\S+)_pileup.tsv"),
         add_inputs(idx_stats_post),
         r"\1_top50_coverage.tsv")
def coverage_plot(infiles, outfile):
    ''' '''

    coverage, idx = infiles

    ModuleTrna.coverage(idx, coverage, outfile)
'''

'''
@follows(count_features)
@follows(mkdir("plots.dir"))
@merge(count_features, "plots.dir/feature_count.png")
def feature_count_plot(infiles, outfile):
    '''''' Create plot of proportion of gene type mapped for each sample''''''

    infiles = ":".join(infiles)
    statement = """Rscript %(cribbslab)s/R/feature_count_plots.r
                --input=%(infiles)s 
                --output=%(outfile)s
                """
    P.run(statement)

'''

@follows(strand_specificity, count_reads, count_features, build_bam_stats,
         full_genome_idxstats, build_samtools_stats, genome_coverage,
         bowtie_index_artificial, index_trna_cluster, remove_reads,
         keep_mature_trna, merge_idx_stats, create_coverage, filter_vcf)
def full():
    pass

@originate("multiqc_data")
def run_multiqc(outfile):
    ''' Run multiqc and overwrite any old reports '''

    statement = '''multiqc -f . '''

    P.run(statement)

@follows(run_multiqc)
@originate("QC_report.html")
def run_rmarkdown(outfile):

    cwd = os.getcwd()
    statement = '''cp %(cribbslab)s/R/QC_report.Rmd . | Rscript -e "rmarkdown::render('QC_report.Rmd', clean=TRUE)" '''

    P.run(statement)


@follows(run_multiqc, run_rmarkdown)
def build_report():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
