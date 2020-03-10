import random
import numpy as np
from utils import read_pdb
from typing import List
import copy

CHARGED = "RHKDE"
NON_CHARGED = "STNQCGPAVILMFYW"


class Gene:
    def __init__(self, value: str = "x", position: int = 0,
                 cros_prob: float = 0.0, mut_prob: float = 0.0,
                 coordinates: np.ndarray = np.array([0.0, 0.0, 0.0])) -> None:
        self._value: str = value
        self._cros_prob: float = cros_prob
        self._mut_prob: float = mut_prob
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


class Protein:
    def __init__(self, sequence) -> None:
        self.__original_sequence = sequence

        self._genes = dict()
        self._value = None

        for n, x in enumerate(sequence):
            gene = Gene(value=x)
            self._genes[n] = gene

    def __getitem__(self, item):
        return self._genes[item]

    @property
    def charge(self):
        charge = 0
        for g in self.genes.values():
            if g.value in "KRH":
                charge += 1
            if g.value in "DE":
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

    def update_genes_from_sequence(self, sequence: str) -> None:
        assert len(self.genes) == len(sequence)
        self.value = None
        for x, y in zip(self.genes, sequence):
            x.value = y

    def read_coordinates(self, fname):
        """
        Считывает координаты аминокислот
        :param fname:
        :return:
        """
        molecule = read_pdb(fname)
        for a in molecule.values():
            if a.Name == "CA":
                resseq = a.ResSeq
                coords = a.Coordin
                self.genes[resseq - 1].coordinates = coords


class IsStableStruct:
    def __init__(self, max_charge=10, max_charged=30, min_distance=5.0):
        self.__max_charge = max_charge
        self.__max_charged = max_charged
        self.__min_distance = min_distance

    def __call__(self, p: Protein):
        try:
            self.max_charge(p)
            self.max_n_charged(p)
            self.distances(p)
            return True
        except:
            return False

    def max_charge(self, p: Protein):
        if abs(p.charge) > self.__max_charge:
            raise Exception

    def max_n_charged(self, p: Protein):
        n = len(p.genes)
        c = 0
        for i in range(0, n):
            if p[i].charged:
                c += 1
        if c > self.__max_charged:
            raise Exception

    def distances(self, p: Protein):
        n = len(p.genes)
        for i in range(0, n):
            for j in range(i + 1, n):
                if p[i].charged and p[j].charged:
                    dist = np.linalg.norm(p[i].coordinates - p[j].coordinates)
                    if dist < self.__min_distance:
                        raise Exception
