

# Import the Biopython functions we need.
# Remember that you can look up BioPython functions here - https://biopython.org/docs/latest/api/Bio.html
from Bio import SeqIO, pairwise2
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

# PrettyTable is a nice package for formatting tables - https://pypi.org/project/prettytable/
import prettytable
from Bio.SeqRecord import SeqRecord
from prettytable import PrettyTable

import datetime
import random
import statistics

from Bio.Seq import Seq
from Bio import pairwise2 as pw
from Bio import AlignIO

# PLEASE DONT FORGET TO SET YOUR EMAIL THIS IS GOOD ETIQUETTE FOR A SERVICE USER
EMAIL = 'jamesalastairmaxwell@googlemail.com'


def Part1(ids):
	# Initialize PrettyTable
	task1 = PrettyTable()
	task2 = PrettyTable()
	task3 = PrettyTable()
	part2 = PrettyTable()

	# Define the column names for your table
	task1.field_names = ["Gene ID", "Accession", "Nucleotides Length", "3-UTR Length"]
	task2.field_names = ["Gene ID", "Accession", "Nucleotides Length", "% A", "% C", "% G", "% T"]
	task3.field_names = ["Gene ID", "Accession", "Nucleotides Length", "Protein Length", "Most Frequent Amino Acid"]
	part2.field_names = ["Gene ID", "Accession", "Nucleotides Length", "Exons"]

	# Iterate over each gene ID and process the information
	for gene_id in ids:
		Entrez.email = EMAIL
		handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
		gene = SeqIO.read(handle, "genbank")
		handle.close()

		nucleotides_len = len(gene.seq)

		a = c = g = t = 0

		for n in gene.seq:
			if n == 'A':
				a += 1
			if n == 'C':
				c += 1
			if n == 'G':
				g += 1
			if n == 'T':
				t += 1

		three_utr_len = None  # Default, in case we don't find a 'CDS'
		protein_sequence = None
		exons = 0

		# Find the coding sequence location again
		for f in gene.features:
			if f.type == 'CDS':
				# Find the location of the coding sequence
				coding_location = f.location
				# Extract the 3'-UTR by taking all the sequence after the coding sequence
				seq3utr = gene.seq[coding_location.end:]
				three_utr_len = len(seq3utr)

				protein_sequence = f.translate(gene).seq

			if f.type == "exon":
				exons += 1

		# Add a row with the data for this gene
		task1.add_row([
			gene_id,
			gene.annotations['accessions'][0] if 'accessions' in gene.annotations else 'N/A',
			# Assuming we need the first accession ID
			nucleotides_len,
			three_utr_len if three_utr_len is not None else 'N/A'
		])

		task2.add_row([
			gene_id,
			gene.annotations['accessions'][0] if 'accessions' in gene.annotations else 'N/A',
			# Assuming we need the first accession ID
			nucleotides_len,
			str(round(100 * a / nucleotides_len)) + '%',
			str(round(100 * c / nucleotides_len)) + '%',
			str(round(100 * g / nucleotides_len)) + '%',
			str(round(100 * t / nucleotides_len)) + '%'
		])

		task3.add_row([
			gene_id,
			gene.annotations['accessions'][0] if 'accessions' in gene.annotations else 'N/A',
			# Assuming we need the first accession ID
			nucleotides_len,
			len(protein_sequence),
			max(set(protein_sequence), key=protein_sequence.count)
		])

		part2.add_row([
			gene_id,
			gene.annotations['accessions'][0] if 'accessions' in gene.annotations else 'N/A',
			# Assuming we need the first accession ID
			nucleotides_len,
			exons
		])

	# Print the table to the console
	print(task1)
	print(rn)
	print("\n")
	print(task2)
	print(rn)
	print("(Due to rounding, the displayed percentages may not total exactly 100%)\n\n")
	print(task3)
	print(rn)

	return part2


def Part2(ids):
	Entrez.email = EMAIL
	handle = Entrez.efetch(db="nucleotide", id=1676318229, rettype="gb", retmode="text")
	gene = SeqIO.read(handle, "genbank")
	handle.close()
	shortest = gene.seq
	shortest_protein_sequence = None
	shortest_exons = []
	# Find the coding sequence location again
	for f in gene.features:
		if f.type == 'CDS':
			shortest_protein_sequence = f.translate(gene).seq
		if f.type == "exon":
			shortest_exons.append(gene.seq[f.location.start:f.location.end])

	handle = Entrez.efetch(db="nucleotide", id=1676317140, rettype="gb", retmode="text")
	gene = SeqIO.read(handle, "genbank")
	handle.close()
	longest = gene.seq
	longest_protein_sequence = None
	longest_exons = []
	# Find the coding sequence location again
	for f in gene.features:
		if f.type == 'CDS':
			longest_protein_sequence = f.translate(gene).seq
		if f.type == "exon":
			longest_exons.append(gene.seq[f.location.start:f.location.end])

	alignments = pw.align.globalxx(shortest, longest)
	protein_alignments = pw.align.globalxx(shortest_protein_sequence, longest_protein_sequence)

	# distribution = generate_score_distribution(shortest, longest, "nucleotide_score_distro.txt")
	# protein_distribution = generate_score_distribution(shortest_protein_sequence, longest_protein_sequence, "protein_score_distro.txt")
	distribution = read_score_distribution("nucleotide_score_distro.txt")
	protein_distribution = read_score_distribution("protein_score_distro.txt")

	# the alignment score
	print("Nucleotide alignment score =", alignments[0][2])
	print("Nucleotide alignment Z-score =",
		  (alignments[0][2] - statistics.stdev(distribution)) / statistics.mean(distribution))
	print("Protein isoform alignment score =", protein_alignments[0][2])
	print("Protein isoform alignment Z-score =",
		  (protein_alignments[0][2] - statistics.stdev(protein_distribution)) / statistics.mean(protein_distribution))

	exons = {"declare"}
	exons.remove("declare")
	for i in ids:
		handle = Entrez.efetch(db="nucleotide", id=i, rettype="gb", retmode="text")
		gene = SeqIO.read(handle, "genbank")
		for f in gene.features:
			if f.type == "exon":
				exons.add(gene.seq[f.location.start:f.location.end])

	print("Number of exons:", len(exons))
	print("Exon lengths:")
	print([len(e) for e in exons])

	extras = []
	for exon in longest_exons:
		if exon not in shortest_exons:
			extras.append(exon)
	print("Exons in the longest transcript but not the shortest:")
	print(extras)

	generate_blastx(extras)


def Part3():
	Entrez.email = EMAIL
	cadherin_proteins = {"declare"}
	cadherin_proteins.remove("declare")
	for i in range(1, 50):
		query = "\"Homo sapiens\"[Organism] AND CDH" + str(i) + "[Gene] AND swissprot[filter]"
		handle = Entrez.esearch(db="protein", term=query, retmax=100)
		record = Entrez.read(handle)
		handle.close()

		# print(query)
		# print(record["IdList"])

		cadherin_proteins.update(record["IdList"])

	print(cadherin_proteins)  # maybe get accessions from this

	CDH7_id = 296434420  # id for cadherin-7
	cadherin_proteins.remove(str(CDH7_id))
	handle = Entrez.efetch(db="protein", id=CDH7_id, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()
	base_seq = record.seq
	cadherin_alignments = []

	for p in cadherin_proteins:
		handle = Entrez.efetch(db="protein", id=p, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		handle.close()
		comp_seq = record.seq

		cadherin_alignments.append((p, pw.align.globalms(base_seq, comp_seq, 1, -1, -1, -1)[0].score))

	# print(cadherin_alignments)

	cadherin_alignments.sort(key=lambda x: x[1])
	for p in cadherin_alignments:
		print(p)

	closest_cadherins = [('190359622', 'CDH20'), ('116241276', 'CDH10'), ('3023435', 'CDH18')]#
	print("The Cadherins most similar to Cadherin-7 are:")
	print(closest_cadherins)
	closest_cadherin_isoform = '37999814'

	handle = Entrez.efetch(db="protein", id=p, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()
	comp_seq = record.seq


# # Write the sequence record to a FASTA file
# with open("output.fasta", "w") as output_handle:
# 	SeqIO.write(record, output_handle, "fasta")


def generate_blastx(extras):
	i = 0
	print(i)
	for e in extras:
		result_handle = NCBIWWW.qblast("blastx", 'nr', e)
		blast_record = NCBIXML.read(result_handle)

		result_handle_table(blast_record, len(e))
		i += 1


def generate_blastp(proteins):
	i = 0
	# print(i)
	for p in proteins:
		Entrez.email = EMAIL
		handle = Entrez.efetch(db="protein", id=p, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		handle.close()

		result_handle = NCBIWWW.qblast("blastp", 'nr', record.seq)
		blast_record = NCBIXML.read(result_handle)

		result_handle_table(blast_record, len(p))
		i += 1


def shuffle_sequence(seq):
	# Shuffles a given sequence.
	seq_list = list(seq)
	random.shuffle(seq_list)
	return ''.join(seq_list)


def generate_score_distribution(seq1, seq2, filename, num_shuffles=1000):
	# Generates a distribution of alignment scores.
	scores = []
	f = open(filename, "w")
	for i in range(0, num_shuffles):
		shuffled_seq2 = shuffle_sequence(seq2)
		alignment = pw.align.globalxx(seq1, shuffled_seq2, one_alignment_only=True)
		score = alignment[0].score
		scores.append(score)
		print(i)
		f.write(str(score))
		f.write(",")
		print(datetime.datetime.now())
	f.close()
	return scores


def read_score_distribution(filename):
	f = open(filename, "r")
	raw_text = f.read()
	# print(raw_text)
	distribution = raw_text.split(",")
	distribution.pop()
	scores = []
	for s in distribution:
		scores.append(float(s))
	return scores


def result_handle_table(blast_record, query_len):
	t = PrettyTable(
		['Description', 'Max Score', 'Total Score', 'Query Cover', 'E Value', 'Per Ident', 'Acc Len', 'Accession'])
	t._max_width = {"Description": 30}
	for alignment in blast_record.alignments:
		score = 0
		max_score = 0
		query_cover = 0
		perc_ident = 100
		for hsp in alignment.hsps:
			score = score + hsp.bits
			max_score = max(max_score, hsp.bits)
			query_cover += ((hsp.query_end - hsp.query_start) / query_len)
			perc_ident = min(perc_ident, hsp.identities / hsp.align_length)
		if hsp.expect < 0.05:
			t.add_row(
				[alignment.hit_def.split('>')[0], round(max_score), round(score), '{:3.0f}%'.format(query_cover * 100),
				 '{:1.2e}'.format(hsp.expect), '{:3.2f}%'.format(perc_ident * 100), alignment.length,
				 alignment.accession])

	return print(t[0:10])


transcript_ids = [1519244509, 1890333824, 1676318229, 1676317140]
rn = datetime.datetime.now

exonTable = Part1(transcript_ids)
exonTable.sortby = "Nucleotides Length"
print(exonTable)
print(rn)

Part2(transcript_ids)

Part3()
