#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
03_single_annotation.py
Annotates chr1 using a chr1-only GFF3.
Outputs a tab-separated summary of gene features with unique CDS.
"""

from Bio import SeqIO
import gffutils
import os

# Paths
fasta_path = "data/clean/chr1.fasta"       # cleaned chr1
output_path = "data/annotated_chr1.tsv"

gff3_chr1 = "data/clean/chr1.gff3"

db_chr1 = "data/clean/chr1.db"                   # SQLite database for gffutils

# Create GFF database (if it doesn't exist)
if not os.path.exists(db_chr1):
    print("Creating chr1-only GFF database (fast, ignoring duplicates)...")
    # Allow duplicates by using merge strategy 'merge'
    gffutils.create_db(
        gff3_chr1,
        dbfn=db_chr1,
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )

# Load the database
db = gffutils.FeatureDB(db_chr1, keep_order=True)

# Load chromosome sequence
chr_record = SeqIO.read(fasta_path, "fasta")

# Prepare output
with open(output_path, "w") as out_file:
    out_file.write("gene_id\tgene_name\tstart\tend\tstrand\tfeature_type\n")
    
    # Iterate over all genes on chr1
    for feature in db.features_of_type("gene"):
        # Skip if gene is not on chr1
        if feature.seqid != chr_record.id:
            continue

        gene_id = feature.attributes.get("gene_id", ["NA"])[0]
        gene_name = feature.attributes.get("gene_name", ["NA"])[0]

        # Write gene
        out_file.write(f"{gene_id}\t{gene_name}\t{feature.start}\t{feature.end}\t{feature.strand}\tgene\n")

        # Track unique CDS per gene
        seen_cds = set()
        for cds in db.children(feature, featuretype="CDS", order_by="start"):
            key = (cds.start, cds.end, cds.strand)
            if key not in seen_cds:
                out_file.write(f"{gene_id}\t{gene_name}\t{cds.start}\t{cds.end}\t{cds.strand}\tCDS\n")
                seen_cds.add(key)

print(f"Annotation completed! Results saved to {output_path}")
