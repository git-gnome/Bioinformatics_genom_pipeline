

import os
from collections import Counter
from Bio import SeqIO
from pyfaidx import Fasta
import shutil

# Paths
fasta_path = "/home/gnome_/Work/Bioinformatics/data/raw/hg38.fa"
data_dir = "/home/gnome_/Work/Bioinformatics/data"
raw_dir = os.path.join(data_dir, "raw")
clean_dir = os.path.join(data_dir, "clean")

# Create subfolders
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(clean_dir, exist_ok=True)

# Move the raw genome to the raw folder
shutil.copy(fasta_path, raw_dir)

# Load genome with pyfaidx
genome = Fasta(fasta_path)
total_length = 0
base_counts = Counter()

for chrom in genome.keys():
    seq = genome[chrom][:].seq.upper()
    total_length += len(seq)
    base_counts.update(seq)

print(f"Total genome length: {total_length:,} bp")
print("Base composition:")
for base in "ACGTN":
    print(f"{base}: {base_counts.get(base, 0):,}")

# Primary chromosomes
primary_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
primary_seqs = []

# Process genome and save primary chromosomes
for record in SeqIO.parse(fasta_path, "fasta"):
    seq = record.seq.upper()
    counts = Counter(seq)
    print(f"{record.id}: {len(seq):,} bp, N={counts.get('N',0):,}")

    if record.id in primary_chroms:
        # Save individual chromosome in clean folder
        chrom_file = os.path.join(clean_dir, f"{record.id}.fasta")
        SeqIO.write(record, chrom_file, "fasta")
        primary_seqs.append(record)

# Save combined clean genome
clean_genome_file = os.path.join(clean_dir, "hg38_primary.fasta")
SeqIO.write(primary_seqs, clean_genome_file, "fasta")

print(f"Clean genome saved in: {clean_genome_file}")
print(f"Individual chromosomes saved in: {clean_dir}")
