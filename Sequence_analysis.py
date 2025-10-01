from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import gc_fraction
from collections import Counter
import matplotlib.pyplot as plt

#---------------------A. DNA Sequence Analysis-------------------
# --------------------1. Reading a FASTA File--------------------
print("Reading FASTA File")
records = list(SeqIO.parse("Fasta_Seq.fasta", "fasta"))
for rec in records:
    print(f"ID: {rec.id}, Length of Sequence: {len(rec.seq)}, Sequence: {rec.seq}")

# -------------------2. Computing GC Content ----------------------
print("The GC content of FASTA File is:\n")
for recs in records:
    print(f"ID: {recs.id} | GC: {gc_fraction(recs.seq) * 100:.2f}%")

# Visualization
print("The plot for the GC content is:\n")
gc_contents = [gc_fraction(recs.seq) * 100 for recs in records]

plt.hist(gc_contents, bins=20, color="skyblue", edgecolor="black")
plt.title("GC-Content Distribution")
plt.xlabel("GC%")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()

# --------------------3. Motif Search ------------------------------
motifs = ["ATG", "TATA", "CGCG"]
m_counts = {m: 0 for m in motifs}
seq_num = 1
print("The Occurence of Motifs: ATG, TATA, CGCG, in the sequences are:")
for rec in records:
    seq = str(rec.seq)
    print(f"Sequence Number: {seq_num}")
    for m in motifs:
        count = seq.count(m)
        m_counts[m] += count
        print(f" {m} : {count}")
    seq_num += 1


#Visualization
print("The plot for Motif Frequency is:")
plt.bar(m_counts.keys(), m_counts.values(), color= ['orange', 'pink', 'green'])
plt.title("Motif Occurences in DNA Sequences")
plt.xlabel("Motif")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()

# ------------------------4. DNA Blast---------------------------
Sequence = records[0].seq
results = NCBIWWW.qblast("blastn", "nt", Sequence)
with open("blast_results.xml", "w") as file:
    file.write(results.read())
results.close()
with open("blast_results.xml") as results_handle:
    blast_records = NCBIXML.parse(results_handle)
    blast_record = list(blast_records)

# Number of Good Hits
E_VALUE_THRESH = 0.00000000001 
count = 0
for record in blast_record:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                count += 1
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                print()
                
print(f"There are {count} similar sequences in Blast output")

# -----------------B. Protein Analysis----------------------
# -----------------1. Transcription & Translation -------------------
output_file = "Protein_Seq.fasta"

with open(output_file, "w") as file:
    for rec in records:
        rna_seq = rec.seq.transcribe()
        protein_seq = rna_seq.translate()

        file.write(f">{rec.id}\n{protein_seq}\n")
print(f"Translated Sequences are saved in {output_file}")

# ------------------2. Common A.A Frequency ----------------------------

protein_records = list(SeqIO.parse("Protein_Seq.fasta", "fasta"))
all_proteins = "".join([str(rec.seq) for rec in protein_records])
amino_acid_count = Counter(all_proteins)
amino_acid_count.pop('*', None)
amino_acid_count.most_common(10)
amino_list = dict(sorted(amino_acid_count.items()))

#Visualization
plt.bar(amino_list.keys(), amino_list.values(), width= 0.5, color= ['b', 'r', 'm', 'c'])
plt.title("Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()

# ---------------------3. Protein Blast --------------------------------------

Sequence1 = protein_records[0].seq
results = NCBIWWW.qblast("blastp", "pdb", Sequence1)
with open("pblast_results.xml", "w") as file:
    file.write(results.read())
results.close()
with open("pblast_results.xml") as results_handle:
    blast_recs = NCBIXML.parse(results_handle)
    blast_rec = list(blast_recs)

# Number of Good Hits
E_VALUE_THRESH = 0.00000000001 
count = 0
for record in blast_rec:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                count += 1
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                print()
                
print(f"There are {count} similar sequences in Blast output")

