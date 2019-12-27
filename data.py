import random
from typing import List

class Gene():
    def __init__(self, position: int = None, variants: str = None,
                 cros_prob: float = None, mut_prob: float = None) -> None:
        self.position = position
        self.variants = variants
        self.cros_prob = cros_prob
        self.mut_prob = mut_prob

    def select_aminoacid(self) -> str:
        x = random.choice(self.variants)
        return x

class Protein():
    def __init__(self, genes) -> None:
        self.variance = dict()
        self.mlambda = None
        for gene in genes:
            self.variance[gene.position] = gene.select_aminoacid()

    def update_variance_from_sequence(self, genes: List[Gene], sequence: str) -> None:
        self.variance = dict()
        self.mlambda = None
        for gene in genes:
            pos = gene.position
            self.variance[pos] = sequence[pos - 1]

    def get_sequence(self, pattern_sequence: str) -> str:
        pattern_sequence = list(pattern_sequence)
        for (position, aminoacid) in self.variance.items():
            pattern_sequence[position - 1] = aminoacid
        return ''.join(pattern_sequence)

    def update_variance(self, position: int, aminoacid: str) -> None:
        self.variance[position] = aminoacid
        self.mlambda = None

    def get_variance(self) -> str:
        return ''.join(self.variance.values())

    def set_variance(self, aminoacids: str) -> None:
        positions = self.variance.keys()
        for position, aminoacid in zip(positions, aminoacids):
            self.update_variance(position, aminoacid)