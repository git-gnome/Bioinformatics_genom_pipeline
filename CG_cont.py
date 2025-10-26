
# source /home/gnome_/Work/Bioinformatics/bioenv/bin/activate
# source bioenv/bin/activate
# python CG_cont.py

# cg_content_chr1.py
from Bio import SeqIO
from collections import Counter

# Path to your cleaned chromosome 1 fasta
fasta_path = "/home/gnome_/Work/Bioinformatics/data/clean/chr1.fasta"

# Parse the fasta
for record in SeqIO.parse(fasta_path, "fasta"):
    seq = record.seq.upper()
    counts = Counter(seq)

    total_length = len(seq)
    gc_count = counts.get("G", 0) + counts.get("C", 0)
    cg_content = gc_count / total_length * 100

    print(f"Chromosome: {record.id}")
    print(f"Total length: {total_length:,} bp")
    print("Base composition:")
    for base in "ACGTN":
        print(f"  {base}: {counts.get(base,0):,}")
    print(f"GC content: {cg_content:.2f}%")
