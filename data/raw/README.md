# Data

This directory is intended to store the input files required for the RNA-seq analysis pipeline.

## Directory Structure

```
data/
├── raw/          # Raw paired-end FASTQ files (not included)
└── reference/    # Reference genome and annotation files
```

## Raw Sequencing Data

The raw FASTQ files used for this project are **not included** in this repository due to their large size and potential data-sharing restrictions.

Place your paired-end FASTQ files inside the `raw/` directory before running the pipeline.

Example:

```
raw/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

## Reference Files

The `reference/` directory contains the reference genome and annotation files required for alignment and downstream analysis.

Typical files include:

- Genome sequence (`.fna` / `.fa`)
- Gene annotation (`.gff` / `.gtf`)
- STAR genome index (generated before alignment)

> **Note:** If the reference files are too large, they should be downloaded separately and are not typically stored in GitHub repositories.
