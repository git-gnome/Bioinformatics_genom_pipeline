#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
03_single_annotation.py
Annotates chr1 using Gencode v49 basic GFF3.
Duplicates are ignored to speed up DB creation.
"""

from Bio import SeqIO
import gffutils
import os

# Paths
fasta_path = "data/clean/chr1.fasta"
output_path = "data/clean/annotated_chr1.tsv"
gff3_chr1 = "data/clean/chr1.gff3"  # extracted chr1 GFF3
db_chr1 = "data/clean/chr1.db"

# Create GFF database if it doesn't exist
if not os.path.exists(db_chr1):
    print("Creating chr1-only GFF database (fast, ignoring duplicates)...")
    gffutils.create_db(
        gff3_chr1,
        dbfn=db_chr1,
        force=True,
        keep_order=True,
        merge_strategy='create_unique',  # ignore duplicate IDs
        sort_attribute_values=True
    )

# Load the database
db = gffutils.FeatureDB(db_chr1, keep_order=True)

# Load chromosome sequence
chr_record = SeqIO.read(fasta_path, "fasta")

# Write annotations to TSV
with open(output_path, "w") as out_file:
    out_file.write("gene_id\tgene_name\tstart\tend\tstrand\tfeature_type\n")

    # Iterate over all genes on this chromosome
    for feature in db.region(seqid=chr_record.id, featuretype="gene"):
        gene_id = feature.attributes.get("gene_id", ["NA"])[0]
        gene_name = feature.attributes.get("gene_name", ["NA"])[0]
        out_file.write(f"{gene_id}\t{gene_name}\t{feature.start}\t{feature.end}\t{feature.strand}\tgene\n")

        # List CDS or exons
        for cds in db.children(feature, featuretype="CDS", order_by="start"):
            out_file.write(f"{gene_id}\t{gene_name}\t{cds.start}\t{cds.end}\t{cds.strand}\tCDS\n")

print(f"Annotation completed! Results saved to {output_path}")
