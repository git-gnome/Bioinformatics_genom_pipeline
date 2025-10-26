#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
02_individual_chr_gff3.py
Extracts a single chromosome from a GFF3 file and saves it as a separate GFF3.
"""

import os
import sys

# ===== USER SET VARIABLES =====
chrom = sys.argv[1]           # Chromosome name (e.g., chr1)
gff3_full = sys.argv[2]       # Full GFF3 path (e.g., data/raw/gencode.v49.basic.annotation.gff3)
clean_dir = sys.argv[3]       # Output folder (e.g., data/clean)
# =============================

out_gff3 = os.path.join(clean_dir, f"{chrom}.gff3")

# Ensure output folder exists
os.makedirs(clean_dir, exist_ok=True)

# Skip if already exists
if os.path.exists(out_gff3):
    print(f"{chrom}.gff3 already exists, skipping extraction.")
else:
    print(f"Extracting {chrom} from full GFF3...")
    with open(gff3_full, "r") as infile, open(out_gff3, "w") as outfile:
        for line in infile:
            if line.startswith("#") or line.startswith(f"{chrom}\t"):
                outfile.write(line)
    print(f"Saved {chrom} GFF3 -> {out_gff3}")
