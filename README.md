# Genomic_pipeline

A reproducible bioinformatics pipeline for processing genomic annotation and sequence data.

## Overview  
This pipeline provides scripts to clean, annotate, and separate genomic-data files by chromosome. The main goal is to maintain a clear folder structure and modular scripts, while keeping large raw/clean datasets **outside** the version control repository.

## Repository Structure  

Genomic_pipeline/
├── 00_main_pipeline.py              # Main driver script
├── 01_main_chr_clean.py             # Chromosome cleaning script
├── 02_individual_chr_gff3.py           # Script for per-chromosome GFF3 processing
├── 03_single_annotation.py                 # Single annotation script
├── 03_single_annotation_unique.py          # Variation of annotation script
├── CG_cont.py              # Script for CG content analysis
├── chr1_annotation.py              # Example script for chromosome 1
├── chromosome_clean.py                 # Auxiliary cleaning script
├── chromosome_sep.py                # Chromosome separation script
├── data/
│ ├── raw/ # (Empty placeholder) raw large datasets
│ │ └── .gitkeep
│ └── clean/ # (Empty placeholder) processed data
│ └── .gitkeep
├── .gitignore # Ignore large files, virtual-env etc.
└── README.md # This file


## How It Works  
1. You place large reference and annotation files in `data/raw/` (e.g., genome files, GFF3 annotation).  
2. The script `01_main_chr_clean.py` and others operate on those raw data to produce cleaned files in `data/clean/`.  
3. Auxiliary scripts further separate by chromosome (`chromosome_sep.py`), compute annotation per chromosome (`02_individual_chr_gff3.py`), and refine into unique annotations (`03_single_annotation_unique.py`).  
4. The main driver `00_main_pipeline.py` ties the steps together.  
5. You should **not** commit large data files into the repo; this repository only tracks the code and folder structure for reproducibility.

## Prerequisites  
- Python 3.x  
- Required Python packages: (list packages here, e.g. pandas, gffutils, etc.)  
  E.g.:  
  ```bash
  pip install pandas gffutils numpy

Usage

# Activate your environment
source bioenv/bin/activate

# Ensure your raw data is copied into data/raw/
cp /path/to/hg38.fa data/raw/
cp /path/to/gencode.v49.basic.annotation.gff3 data/raw/

# Run the main pipeline
python 00_main_pipeline.py