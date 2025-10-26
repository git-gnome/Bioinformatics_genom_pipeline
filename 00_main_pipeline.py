#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
00_main_pipeline.py
Main pipeline to extract chromosome FASTA, GFF3, and annotate it.
"""

import os
import subprocess
import sys

# ===== USER SET VARIABLES =====
chrom = "chr1"                             # Chromosome to process
raw_fasta = "data/raw/hg38.fa"             # Full genome FASTA
gff3_full = "data/raw/gencode.v49.basic.annotation.gff3"  # Full GFF3
clean_dir = "data/clean"                   # Clean folder
# =============================

script_dir = os.getcwd()

# 1. Extract chromosome FASTA
print(f"[Step 1] Extracting {chrom} FASTA...")
subprocess.run([
    sys.executable,
    os.path.join(script_dir, "01_main_chr_clean.py"),
    chrom, raw_fasta, clean_dir
])

# 2. Extract chromosome GFF3
print(f"[Step 2] Extracting {chrom} GFF3...")
subprocess.run([
    sys.executable,
    os.path.join(script_dir, "02_individual_chr_gff3.py"),
    chrom, gff3_full, clean_dir
])

# 3. Annotate chromosome
print(f"[Step 3] Annotating {chrom}...")
fasta_chr = os.path.join(clean_dir, f"{chrom}.fasta")
gff_chr = os.path.join(clean_dir, f"{chrom}.gff3")

subprocess.run([
    sys.executable,
    os.path.join(script_dir, "03_single_annotation_unique.py"),
    chrom, fasta_chr, gff_chr, clean_dir
])

print("Pipeline completed successfully!")
