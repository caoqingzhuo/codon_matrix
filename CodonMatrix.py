import torch
from Bio import SeqIO

Codons = [a+b+c for a in 'ATCG' for b in 'ATCG' for c in 'ATCG']
CodonDict = {codon: index for index, codon in enumerate(Codons)}
print(CodonDict)

class CodonMatrix:
	def __init__(self, transcript, CodonDict = CodonDict):
		self.CodonDict = CodonDict
		self.codon_count = torch.zeros(64, 64, dtype=torch.int32)  # Integer matrix for codon counts
		self.transcript = transcript
		self.CalCodonCount

	def CalCodonCount(self):
        # Calculate codon pair counts
		for i in range(0, len(self.transcript) - 6, 3):
			first_codon = self.transcript[i:i + 3]
			second_codon = self.transcript[i + 3:i + 6]
			if first_codon in self.CodonDict and second_codon in self.CodonDict:
				first_index = self.CodonDict[first_codon]
				second_index = self.CodonDict[second_codon]
				self.codon_count[first_index][second_index] += 1
		return self.codon_count


class CodonMatrixs:
	def __init__(self, mRNAsDict, CodonDict):
		self.codon_count = torch.zeros(64, 64, dtype=torch.int32)  # Integer matrix for aggregated codon counts
		self.transition_probabilities = torch.zeros(64, 64, dtype=torch.float32)
		for key,val in mRNAsDict.items():
			transcript = str(val.seq)
			codon_matrix = CodonMatrix(transcript).CalCodonCount()
			self.update_codon_count(codon_matrix)
		self.calculate_transition_probabilities()

	def update_codon_count(self, codon_count):
		self.codon_count = self.codon_count + codon_count

	def calculate_transition_probabilities(self):
		row_sums = self.codon_count.sum(dim=1, keepdim=True).float()
		self.transition_probabilities = torch.where(row_sums != 0, self.codon_count.float() / row_sums, torch.zeros_like(self.codon_count).float())
		return self.transition_probabilities

fasta_file = "./c_elegans.PRJNA13758.WS285.CDS_transcripts.fa"
mRNA = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

CodonMatrixTransitionP = CodonMatrixs(mRNA, CodonDict).calculate_transition_probabilities()




#mRNAsCondon = CodonMatrix()
