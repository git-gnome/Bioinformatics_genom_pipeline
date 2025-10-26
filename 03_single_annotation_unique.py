#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
03_single_annotation_unique.py
Annotates a single chromosome using a GFF3 file.
Removes duplicate CDS entries for a cleaner summary.
"""

from Bio import SeqIO
import gffutils
import os
import sys

# ===== USER SET VARIABLES =====
chrom = sys.argv[1]           # Chromosome name (e.g., chr1)
fasta_path = sys.argv[2]      # Chromosome FASTA path (e.g., data/clean/chr1.fasta)
gff3_chr = sys.argv[3]        # Chromosome GFF3 path (e.g., data/clean/chr1.gff3)
clean_dir = sys.argv[4]       # Output folder (e.g., data/clean)
# =============================

db_chr = os.path.join(clean_dir, f"{chrom}_unique.db")
output_path = os.path.join(clean_dir, f"annotated_{chrom}_unique.tsv")

# Create GFF database (duplicates handled by 'create_unique')
if not os.path.exists(db_chr):
    print(f"Creating {chrom}-only GFF database (fast, ignoring duplicates)...")
    gffutils.create_db(
        gff3_chr,
        dbfn=db_chr,
        force=True,
        keep_order=True,
        merge_strategy='create_unique',  # avoids duplicate IDs
        sort_attribute_values=True
    )

db = gffutils.FeatureDB(db_chr, keep_order=True)

# Load chromosome sequence
chr_record = SeqIO.read(fasta_path, "fasta")

# Keep track of already written CDS to remove duplicates
seen_cds = set()

# Write annotations
with open(output_path, "w") as out_file:
    out_file.write("gene_id\tgene_name\tstart\tend\tstrand\tfeature_type\n")
    
    # Iterate over all genes on the chromosome
    for gene in db.region(seqid=chr_record.id, featuretype="gene"):
        gene_id = gene.attributes.get("gene_id", ["NA"])[0]
        gene_name = gene.attributes.get("gene_name", ["NA"])[0]
        out_file.write(f"{gene_id}\t{gene_name}\t{gene.start}\t{gene.end}\t{gene.strand}\tgene\n")

        # List CDS or exons, skip duplicates
        for cds in db.children(gene, featuretype="CDS", order_by="start"):
            cds_key = (gene_id, cds.start, cds.end, cds.strand)
            if cds_key not in seen_cds:
                out_file.write(f"{gene_id}\t{gene_name}\t{cds.start}\t{cds.end}\t{cds.strand}\tCDS\n")
                seen_cds.add(cds_key)

print(f"Annotation completed! Results saved to {output_path}")
