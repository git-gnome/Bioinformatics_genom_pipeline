#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
chr1_annotation.py
Annotates chr1 using a chr1-only GFF3 file (Gencode v49 basic).
Outputs a tab-separated summary of gene features with progress updates.
"""

from Bio import SeqIO
import gffutils
import os

# =========================
# Paths
# =========================
fasta_path = "data/clean/chr1.fasta"          # cleaned chr1 fasta
output_path = "data/annotated_chr1.tsv"       # output TSV
gff3_chr1 = "data/chr1.gff3"                 # chr1-only GFF3
db_chr1 = "data/chr1.db"                     # gffutils DB for chr1

# =========================
# Create or load GFF database
# =========================
if not os.path.exists(db_chr1):
    print("Creating chr1-only GFF database (fast)...")
    gffutils.create_db(
        gff3_chr1,
        dbfn=db_chr1,
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )
    print("GFF database created successfully!")

db = gffutils.FeatureDB(db_chr1, keep_order=True)

# =========================
# Load chromosome sequence
# =========================
chr_record = SeqIO.read(fasta_path, "fasta")

# =========================
# Annotate genes
# =========================
with open(output_path, "w") as out_file:
    out_file.write("gene_id\tgene_name\tstart\tend\tstrand\tfeature_type\n")
    
    gene_count = 0
    for gene in db.features_of_type("gene", seqid=chr_record.id):
        gene_id = gene.attributes.get("gene_id", ["NA"])[0]
        gene_name = gene.attributes.get("gene_name", ["NA"])[0]
        out_file.write(f"{gene_id}\t{gene_name}\t{gene.start}\t{gene.end}\t{gene.strand}\tgene\n")
        
        # Optional: annotate CDS or exons
        for cds in db.children(gene, featuretype="CDS", order_by="start"):
            out_file.write(f"{gene_id}\t{gene_name}\t{cds.start}\t{cds.end}\t{cds.strand}\tCDS\n")
        
        gene_count += 1
        # Print progress every 100 genes
        if gene_count % 10 == 0:
            print(f"Annotated {gene_count} genes...")

print(f"Annotation completed! {gene_count} genes annotated.")
print(f"Results saved to {output_path}")
