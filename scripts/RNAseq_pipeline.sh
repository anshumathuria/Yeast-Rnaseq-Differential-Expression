#!/bin/bash

###############################################################
# RNA-seq Preprocessing Pipeline
# Organism: Saccharomyces cerevisiae
# Author: Anshu Mathuria
#
# Steps:
# 1. Quality Control (FastQC)
# 2. MultiQC Summary Report
# 3. Reference Genome Indexing
# 4. Read Alignment using Bowtie2
# 5. (Optional) SAM to BAM Conversion
###############################################################

##############################
# Step 1: Quality Control
##############################

# Run FastQC on all FASTQ files
# Run FastQC on all compressed FASTQ files
fastqc conA_rep1.fq
fastqc conA_rep2.fq
fastqc conB_rep1.fq
fastqc conB_rep2.fq

##############################
# Step 2: Generate MultiQC Report
##############################

# Summarize all FastQC reports into a single HTML report
multiqc .

##############################
# Step 3: Build Bowtie2 Index
##############################

# Build reference genome index (only needs to be done once)
bowtie2-build GCF_000146045.2_R64_genomic.fna Sac_index

##############################
# Step 4: Align Reads
##############################

# Align each sample against the reference genome

bowtie2 -x Sac_index -U conA_rep1.fq -S conA_rep1.sam

bowtie2 -x Sac_index -U conA_rep2.fq -S conA_rep2.sam

bowtie2 -x Sac_index -U conB_rep1.fq -S conB_rep1.sam

bowtie2 -x Sac_index -U conB_rep2.fq -S conB_rep2.sam

##############################
# Step 5: Optional SAM → BAM Conversion
##############################

# Although this workflow continues using SAM files,
# BAM is the preferred format for most RNA-seq pipelines
# because it is compressed, indexed, and computationally efficient.

# Uncomment the following commands if BAM files are required.

# samtools view -bS conA_rep1.sam > conA_rep1.bam
# samtools sort conA_rep1.bam -o conA_rep1_sorted.bam
# samtools index conA_rep1_sorted.bam

# samtools view -bS conA_rep2.sam > conA_rep2.bam
# samtools sort conA_rep2.bam -o conA_rep2_sorted.bam
# samtools index conA_rep2_sorted.bam

# samtools view -bS conB_rep1.sam > conB_rep1.bam
# samtools sort conB_rep1.bam -o conB_rep1_sorted.bam
# samtools index conB_rep1_sorted.bam

# samtools view -bS conB_rep2.sam > conB_rep2.bam
# samtools sort conB_rep2.bam -o conB_rep2_sorted.bam
# samtools index conB_rep2_sorted.bam
