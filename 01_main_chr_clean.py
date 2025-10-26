#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
01_main_chr_clean.py
Extracts a single chromosome from a raw genome FASTA and saves it individually in clean/ folder.
"""

from Bio import SeqIO
import os
import sys

# ===== USER SET VARIABLES =====
chrom = sys.argv[1]        # Chromosome to extract (e.g., chr1)
raw_fasta = sys.argv[2]    # Full genome FASTA (e.g., data/raw/hg38.fa)
clean_dir = sys.argv[3]    # Output folder (e.g., data/clean)
# =============================

# Ensure output folder exists
os.makedirs(clean_dir, exist_ok=True)

out_path = os.path.join(clean_dir, f"{chrom}.fasta")

# Skip if already exists
if os.path.exists(out_path):
    print(f"{chrom}.fasta already exists, skipping extraction.")
else:
    print(f"Processing raw genome for {chrom}...")
    for record in SeqIO.parse(raw_fasta, "fasta"):
        if record.id == chrom:
            SeqIO.write(record, out_path, "fasta")
            print(f"Saved {chrom} -> {out_path}")
            break
