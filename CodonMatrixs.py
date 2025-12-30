from Bio import SeqIO
import torch
import numpy as np

class CodonMatrixs:
    def __init__(self):
        self.codon_matrices = []

    def add_codon_matrix(self, codon_matrix):
        self.codon_matrices.append(codon_matrix)

    def recalculate_transition_probabilities(self):
        if not self.codon_matrices:
            return None

        # Initialize the combined transition counts matrix
        combined_transition_counts = torch.zeros((64, 64))

        # Sum the transition matrices
        for codon_matrix in self.codon_matrices:
            combined_transition_counts += codon_matrix.get_transition_matrix()

        # Calculate the average transition probabilities
        row_sums = combined_transition_counts.sum(dim=1, keepdim=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        combined_transition_probabilities = combined_transition_counts / row_sums

        return combined_transition_probabilities
