#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - 
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-



# Download a sample data from IGSR (https://www.internationalgenome.org/)

wget -P ~/Desktop/code/WGS-Variant-Calling/reads \
      ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P ~/Desktop/code/WGS-Variant-Calling/reads \
      ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


########################## Prep files (TO BE GENERATED ONLY ONCE) ##########################
echo "Run Prep files..."

# Download the human reference genome: hg38

wget -P ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/ \
      https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/Desktop/code/WGS-Variant-Callingo/supporting_files/hg38/hg38.fa.gz


# Index the reference genome to a .fai file using samtools
samtools faidx ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/hg38.fa


# Generate the reference dictionary (a .dict file) using gatk
gatk CreateSequenceDictionary \
      -R ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/hg38.fa \
      -O ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/hg38.dict


# Download known sites files for BQSR from GATK resource bundle
wget -P ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/ \
      https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/Desktop/code/WGS-Variant-Calling/supporting_files/hg38/ \
      https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx



########################## VARIANT CALLING STEPS ##########################


# Setup directories
my_path="/Users/myname/Desktop/code/WGS-Variant-Calling"
ref="${my_path}/supporting_files/hg38/hg38.fa"
known_sites="${my_path}/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="${my_path}/aligned_reads"
reads="${my_path}/reads"
results="${my_path}/results"
data="${my_path}/data"



# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, quality looks okay.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} \
      ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > \
      ${aligned_reads}/SRR062634.paired.sam


#### Look inside the paried mapped results from the terminal:
# samtools view ${aligned_reads}/SRR062634.paired.sam | less

#### counts the number of alignments for each FLAG type
# samtools flagstat ${aligned_reads}/SRR062634.paired.sam   


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam \
      -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam


### Check the number of duplicate reads from the terminal:
#samtools flagstat ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. Build the machine learning model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam \
      -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
# to genereate reads ready to use for variation calling
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} \
      --bqsr-recal-file ${data}/recal_data.table \
      -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} \
      I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam \
      O=${aligned_reads}/alignment_metrics.txt


gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam \
      OUTPUT=${aligned_reads}/insert_size_metrics.txt \
      HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


#### Generate a report from the terminal:
# multiqc .



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam \
      -O ${results}/raw_variants.vcf

# Elapsed time: 160.90 minutes.
# cat ${results}/raw_variants.vcf | less



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf \
      --select-type SNP -O ${results}/raw_snps.vcf

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf \
      --select-type INDEL -O ${results}/raw_indels.vcf

#: << 'END'
#END



