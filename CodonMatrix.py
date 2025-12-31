import torch
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

Codons = [a+b+c for a in 'ATCG' for b in 'ATCG' for c in 'ATCG']
CodonDict = {codon: index for index, codon in enumerate(Codons)}
print(CodonDict)

class CodonMatrix:
    def __init__(self, transcript, CodonDict = CodonDict):
        self.CodonDict = CodonDict
        self.codon_count = torch.zeros(64, 64, dtype=torch.int32)  # Integer matrix for codon counts
        self.transcript = transcript
        self.CalCodonCount()

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

fasta_file = "./p_pacificus_clean_CDS.fasta"
mRNA = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

CodonMatrixTransitionP = CodonMatrixs(mRNA, CodonDict).calculate_transition_probabilities()

#mRNAsCondon = CodonMatrix()
print(CodonMatrixTransitionP)
print(CodonMatrixTransitionP.shape)
row_sums = CodonMatrixTransitionP.sum(dim=1)
print(row_sums.min().item(), row_sums.max().item())
zero_rows = (CodonMatrixTransitionP.sum(dim=1) == 0).sum().item()
print("zero rows:", zero_rows, "/ 64")

M = CodonMatrixTransitionP.detach().cpu().numpy()

# 把 0 mask 掉（不参与颜色）
Mm = np.ma.masked_where(M == 0, M)

# vmin 用非零最小值（但给个下限避免太极端）
vmin = max(Mm.min(), 1e-5)
vmax = Mm.max()

plt.figure(figsize=(7, 6))
plt.imshow(Mm, aspect="equal", interpolation="nearest",
           norm=LogNorm(vmin=vmin, vmax=vmax))
plt.colorbar(label="Transition probability (log, zeros masked)")
plt.title("Codon Transition Matrix (64×64)")
plt.xlabel("Next codon index")
plt.ylabel("Current codon index")
plt.tight_layout()
plt.savefig("heatmap_log_mask0.png", dpi=300)
plt.show()