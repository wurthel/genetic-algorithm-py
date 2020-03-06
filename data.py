import random
import numpy as np
from typing import List

CHARGED = "RHKDE"
NON_CHARGED = "STNQCGPAVILMFYW"


class Gene:
    def __init__(self, value: str = "x", position: int = 0,
                 cros_prob: float = 0.0, mut_prob: float = 0.0,
                 coordinates: np.ndarray = np.array([0.0, 0.0, 0.0])) -> None:
        self._value = value
        self._position = position
        self._cros_prob = cros_prob
        self._mut_prob = mut_prob
        self._coordinates = coordinates

    @property
    def type(self):
        if self.value in CHARGED:
            return "CHARGED"
        if self.value in NON_CHARGED:
            return "NONCHARGED"
        assert False

    @property
    def value(self):
        return self.value

    @value.setter
    def value(self, x: str):
        self._value = x

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, v: np.ndarray):
        self._coordinates = np.array(v)

    @property
    def cros_prob(self):
        return self._cros_prob

    @cros_prob.setter
    def cros_prob(self, x: float):
        self._cros_prob = x

    @property
    def mut_prob(self):
        return self._mut_prob

    @mut_prob.setter
    def mut_prob(self, x: float):
        self._mut_prob = x

    @property
    def position(self):
        return self._position

class Protein:
    def __init__(self, sequence) -> None:
        self.__original_sequence = sequence

        self._genes = []
        self._value = None

        for n, x in enumerate(sequence):
            gene = Gene(value=x, position=n)
            self._genes.append(gene)

    def __getitem__(self, item) -> Gene:
        return self.genes[item]

    @property
    def sequence(self):
        sequence = ""
        for x in self.genes:
            sequence += x.value
        return sequence

    @property
    def genes(self):
        return sorted(self._genes, key=lambda x: x.position)

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, x):
        self._value = x

    def update_genes_from_sequence(self, sequence: str) -> None:
        self.value = None
        for x, y in zip(self.genes, sequence):
            x.value = y
