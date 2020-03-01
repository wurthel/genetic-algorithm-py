import random
from typing import List


class Gene:
    def __init__(self, position: int = None, variants: str = None,
                 cros_prob: float = None, mut_prob: float = None) -> None:
        self._position = position
        self._variants = variants
        self._cros_prob = cros_prob
        self._mut_prob = mut_prob

    def select_aminoacid(self) -> str:
        x = random.choice(self._variants)
        return x

    @property
    def variants(self):
        return self._variants

    @property
    def cros_prob(self):
        return self._cros_prob

    @property
    def mut_prob(self):
        return self._mut_prob

    @property
    def position(self):
        return self._position


class Protein:
    def __init__(self, sequence, genes) -> None:
        self.__original_sequence = sequence

        self._variance = dict()
        self._value = None

        for gene in genes:
            self.variance[gene.position] = gene.select_aminoacid()

    @property
    def variance(self):
        return self._variance

    @property
    def value(self):
        return self._value

    @property
    def sequence(self):
        sequence = list(self.__original_sequence)
        for k, v in self._variance.items():
            sequence[k - 1] = v
        return ''.join(sequence)

    def update_value(self, value):
        self._value = value

    def update_variance(self, position: int, aminoacid: str) -> None:
        self._variance[position] = aminoacid
        self._value = None

    def update_variance_from_sequence(self, genes: List[Gene], sequence: str) -> None:
        for gene in genes:
            k = gene.position
            v = sequence[k - 1]
            self.update_variance(k, v)

    def set_variance(self, positions: List[int], aminoacids: List[str]) -> None:
        for k, v in zip(positions, aminoacids):
            self.update_variance(k, v)
