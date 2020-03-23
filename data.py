import random
import numpy as np
from utils import read_pdb
from typing import List
import copy
from abc import abstractmethod

AMINOACIDS = "RHKDESTNQCGPAVILMFYW"
CHARGED = "RHKDE"
NON_CHARGED = "STNQCGPAVILMFYW"
POSITIVE_CHARGE = "RHK"
NEGATIVE_CHARGE = "DE"


#
# class Chromosome:
#     def __init__(self, genes: List[Gene]):
#         pass
#
#     @property
#     def(self):
#         return


class Gene:
    def __init__(self, value: str = "x",
                 coordinates: np.ndarray = np.array([0.0, 0.0, 0.0])) -> None:
        self._value: str = value
        self._coordinates: np.ndarray = coordinates

    @property
    def charged(self):
        if self.value in CHARGED:
            return True
        return False

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, x: str):
        self._value = x

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, v: np.ndarray):
        self._coordinates = v


class Protein:
    def __init__(self, sequence, value=None) -> None:
        self._genes = dict()
        self._value = value

        for n, x in enumerate(sequence):
            gene = Gene(value=x)
            self._genes[n] = gene

    def __getitem__(self, item):
        return self._genes[item]

    @property
    def charge(self):
        charge = 0
        for g in self.genes.values():
            if g.value in POSITIVE_CHARGE:
                charge += 1
            if g.value in NEGATIVE_CHARGE:
                charge -= 1
        return charge

    @property
    def sequence(self):
        sequence = ""
        for x in self.genes.values():
            sequence += x.value
        return sequence

    @property
    def genes(self):
        return self._genes

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, x):
        self._value = x

    def init_coordinates(self, coords):
        for gene in enumerate(self._genes):
            gene.coordinates = np.array(coords[i])

    def __copy__(self):
        new_protein = Protein(self.sequence, self.value)

        coords = []
        for _, gene in self.genes:
            coords.append(gene.coordinates)
        new_protein.init_coordinates(coords)

        return new_protein