![migex icon](img/migex_icon.png)

# MiGEx — Microbial Genome Explorer

MiGEx (Microbial Genome Explorer) is a lightweight and automated bacterial genome analysis pipeline designed for rapid genome assembly, annotation, antimicrobial resistance detection, taxonomic confirmation, and report generation.

The pipeline focuses on:

* **Fast execution**
* **Simple deployment**
* **Readable outputs for microbiologists and bioinformaticians**
* **Automated PDF report generation**
* **Reproducible analyses using containers**

MiGEx is built with **Snakemake** and combines multiple state-of-the-art microbial genomics tools into a single workflow.

---

# Features

MiGEx performs:

## Quality Control

* Raw read quality assessment with FastQC
* Read filtering and trimming with fastp
* Post-filtering QC validation
* QC summary table generation

## Genome Assembly

* De novo assembly using Unicycler
* Assembly quality assessment with QUAST

## Functional and Genomic Analysis

* Genome annotation with Bakta
* Antimicrobial resistance detection using AMRFinderPlus
* Virulence gene detection with Abricate
* Plasmid detection with Abricate
* MLST sequence typing
* Taxonomic confirmation with Mash
* Assembly coverage analysis using Bowtie2 + Samtools

## Reporting

* Automated PDF report generation using RMarkdown
* Consolidated microbiology-oriented interpretation
* Publication-ready figures and summary tables

---

# Dependencies

MiGEx relies on the following major tools:

| Tool          | Purpose                     |
| ------------- | --------------------------- |
| FastQC        | Read quality control        |
| fastp         | Read trimming and filtering |
| Unicycler     | Genome assembly             |
| QUAST         | Assembly evaluation         |
| Bakta         | Genome annotation           |
| AMRFinderPlus | AMR detection               |
| Abricate      | Virulence/plasmid detection |
| MLST          | Sequence typing             |
| Mash          | Taxonomic identification    |
| Bowtie2       | Read mapping                |
| Samtools      | Coverage statistics         |
| RMarkdown     | PDF report generation       |

Most tools are executed using container images for reproducibility.

---

# Installation

## Clone the repository

```bash
git clone https://github.com/fmauffrey/MiGEx.git
cd MiGEx
```

---

## Create the Conda/Mamba environment

Using Mamba is strongly recommended.

```bash
mamba env create -n migex -f environment.yaml
conda activate migex
```

---

# Database Installation

Coming...

---

# Input Requirements

MiGEx expects paired-end Illumina reads compressed as:

```text
sample_1.fastq.gz
sample_2.fastq.gz
```

All reads should be placed in the same input directory.

---

# Initiate an analysis

The first thing to do to run the pipeline is to run
```text
migex init -i [folder] -o [run_folder] --organism [name]
```

The list of organisms available is available here. This is meant to facilitate the completion of 
the config.file. Only most common organisms are listed. Other Options are available for mlst. See 
the mlst documentation for available values. This can then be changed manually in the config file.

# Configuration

The pipeline uses a `config.yaml` file.

Example:

```yaml
path: "/path/to/fastq"

samples:
  - sampleA
  - sampleB

fastp:
  qualified_quality_phred: 20
  unqualified_percent_limit: 30
  average_quality: 20
  length_limit: 50

unicycler:
  mode: normal

amrfinder:
  organism: ""
  ident_min: 0.8
  coverage_min: 0.8

bakta:
  database: "/path/to/bakta_db"
```

---

# Running MiGEx

## Run the complete workflow
 
 Coming...
---

# Main Outputs

## Quality Control

```text
1-QC/
```

Contains:

* FastQC reports
* fastp HTML reports
* QC summary CSV

---

## Genome Assembly

```text
2-Assembly/
```

Contains:

* Unicycler assemblies
* QUAST metrics
* Assembly logs

---

## Functional Analysis

```text
3-Analysis/
```

Contains:

* AMR predictions
* Genome annotations
* Virulence genes
* Plasmid markers
* MLST typing
* Mash taxonomic results
* Coverage statistics

---

## AMR

Put the list for organism choice. Specify that it is needed to have point mutations.

## MLST

Leave blank if you don't know the species sequenced. The auto mode will select the most appropriate organism. If the result of the identication does not match with the automatically selected organism, you can modify the config file with the appropriate organism, delete the report and the mlst folder and rerun the analysis. Only the mlst analysis will be performed.

## Final Reports

```text
4-Reports/
```

Each sample receives an individual PDF report.

The report summarizes:

* Sequencing quality
* Assembly metrics
* Species identification
* MLST typing
* Resistance genes
* Virulence factors
* Plasmid content
* Coverage information

---

# Current Status

MiGEx is currently under active development.