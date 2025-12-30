import torch
from Bio import SeqIO
from Bio.Data import CodonTable

Codons = [a+b+c for a in 'ATCG' for b in 'ATCG' for c in 'ATCG']
CodonDict = {codon: index for index, codon in enumerate(Codons)}
print(CodonDict)

class CodonMatrix:
    def __init__(self, transcript, CodonDict = CodonDict):
        self.CodonDict = CodonDict
        self.codon_count = torch.zeros(64, 64, dtype=torch.int32)  # Integer matrix for codon counts
        cleaned = CodonMatrix.clean_codons(transcript, table_id=5)
        self.transcript = cleaned if cleaned is not None else ""
        self.CalCodonCount()

    @staticmethod
    def get_start_stop_codons(table_id):
        table = CodonTable.unambiguous_dna_by_id[table_id]
        stop_codons = set(table.stop_codons)
        start_codons = set(table.start_codons)
        return start_codons, stop_codons
    
    @staticmethod
    def clean_codons(seq, table_id):
        start_codons, stop_codons = CodonMatrix.get_start_stop_codons(table_id)
        start_pos = None
        for i in range(0, len(seq) - 2, 3):
            if seq[i:i+3] in start_codons:
                start_pos = i
                break

        if start_pos is None:
            return None

        for i in range(start_pos, len(seq) - 2, 3):
            if seq[i:i+3] in stop_codons:
                cds = seq[start_pos:i+3]
                return cds if len(cds) % 3 == 0 else None
        return None


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

fasta_file = "./mito_CDS.fasta"
mRNA = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

CodonMatrixTransitionP = CodonMatrixs(mRNA, CodonDict).calculate_transition_probabilities()

#mRNAsCondon = CodonMatrix()
print(CodonMatrixTransitionP)
print(CodonMatrixTransitionP.shape)
row_sums = CodonMatrixTransitionP.sum(dim=1)
print(row_sums.min().item(), row_sums.max().item())
zero_rows = (CodonMatrixTransitionP.sum(dim=1) == 0).sum().item()
print("zero rows:", zero_rows, "/ 64")

