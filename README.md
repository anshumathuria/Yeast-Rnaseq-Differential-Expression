# RNA-seq Differential Gene Expression Analysis Pipeline

A reproducible RNA-seq analysis workflow for quality assessment, gene expression quantification, differential expression analysis, and visualization using Bash and R.

---

## Project Overview

This repository contains a complete RNA-seq analysis pipeline developed for differential gene expression (DGE) analysis. The workflow begins with raw paired-end RNA-seq reads, performs quality control and gene quantification, and concludes with statistical analysis and publication-quality visualizations.

The pipeline is designed to be modular, reproducible, and easy to adapt to other RNA-seq datasets.

---

## Workflow

```text
Raw RNA-seq FASTQ Files
          │
          ▼
Quality Assessment
(FastQC)
          │
          ▼
Read Alignment
          │
          ▼
Gene Quantification
(featureCounts)
          │
          ▼
Differential Expression Analysis
(DESeq2)
          │
          ▼
Visualization
(Heatmap • Volcano Plot)
```

---

## Repository Structure

```text
RNAseq-Analysis/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── data/
│   ├── raw/
│   ├── reference/
│   └── README.md
│
├── scripts/
│   ├── RNAseq_pipeline.sh
│   └── RNAseq_analysis.R
│
├── results/
│   ├── fastqc/
│   ├── figures/
│   │   ├── Heatmap.png
│   │   └── Volcano_plot.png
│   └── tables/
│
└── workflow/
    └── pipeline.png
```

---

## Tools Used

| Tool          | Purpose                                    |
| ------------- | ------------------------------------------ |
| FastQC        | Quality assessment of raw sequencing reads |
| featureCounts | Gene expression quantification             |
| DESeq2        | Differential gene expression analysis      |
| R             | Statistical analysis and visualization     |
| Bash          | Pipeline automation                        |

---

## Input

The pipeline requires:

* Paired-end RNA-seq FASTQ files
* Reference genome (FASTA)
* Gene annotation (GTF/GFF)

Raw sequencing data is **not included** in this repository. Place input FASTQ files inside the `data/raw/` directory before running the pipeline.

---

## Output

The pipeline generates:

* FastQC quality reports
* Gene count matrix
* Differentially expressed gene (DEG) table
* Heatmap of differentially expressed genes
* Volcano plot
* Summary tables

All outputs are stored in the `results/` directory.

---

## Running the Pipeline

Run the preprocessing pipeline:

```bash
bash scripts/RNAseq_pipeline.sh
```

Run the downstream differential expression analysis:

```bash
Rscript scripts/RNAseq_analysis.R
```

---

## Applications

This workflow can be adapted for:

* Differential gene expression analysis
* Bulk RNA-seq data analysis
* Comparative transcriptomics
* Functional genomics studies

---

## Author

**Anshu Mathuria**

M.Sc. Biotechnology & Bioinformatics
Institute of Bioinformatics and Applied Biotechnology (IBAB), Bengaluru

---

## License

This project is distributed under the MIT License.
