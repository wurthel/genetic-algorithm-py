import random
import copy
import time
import os
import shutil
from typing import List, Dict
from data import Gene, Protein, CHARGED, NON_CHARGED
from utils import read_sequence, read_coordinates
from abc import abstractmethod
import numpy as np


# def read_genes(config) -> List[Gene]:
#     genes = list()
#     for (position, params) in config.items():
#         params = params.split(' ')
#         position = int(position)
#         variants = params[0]
#         cros_prob = float(params[1])
#         mut_prob = float(params[2])
#         gene = Gene(position=position, variants=variants, cros_prob=cros_prob, mut_prob=mut_prob)
#         genes.append(gene)
#     return genes

class Evolution:
    def __init__(self, population, mut_prob, cros_prob, checker=None):
        self._population = population
        self._mut_prob = mut_prob
        self._cros_prob = cros_prob
        self._checker = checker

    @property
    def population(self):
        return self._population

    @abstractmethod
    def mutation(self):
        pass

    @abstractmethod
    def crossover(self):
        pass

    @abstractmethod
    def selection(self):
        pass


class BaseFunction:
    def __init__(self, input_file, output_file):
        self._input_file = input_file
        self._output_file = output_file
        self._computed = None

    @abstractmethod
    def compute(self):
        pass

    @abstractmethod
    def save_to_file(self, path):
        pass

    @abstractmethod
    def load_from_file(self, path):
        pass


class ProteinEvolution(Evolution, BaseFunction):
    def __init__(self, *args, **kwargs):
        super().__init__(args, kwargs)

        for protein in self._population:
            if protein.value is not None:
                self._computed[protein.sequence] = protein.value


def mutation(self, attempts=1):
    new_population = []
    for protein in self._population:
        r = random.random()
        if r < self._mut_prob:
            while attempts > 0:
                new_protein = copy.copy(protein)
                for gene in new_protein.genes():
                    r = random.random()
                    if r < gene.mut_prob:
                        gene.value = random.choice(NON_CHARGED) if gene.charged else random.choice(CHARGED)
                        new_protein.value = None
                if self._checker(new_protein):
                    new_population.append(new_protein)
                    break
                elif attempts == 1:
                    new_protein = copy.copy(protein)
                    new_population.append(new_protein)
                attempts -= 1
    self._population = new_population


def crossover(self, attempts=1):
    new_population = []
    for_cross = []

    for x in self._population:
        r = random.random()
        if r < self._cros_prob:
            for_cross.append(x)
        else:
            new_population.append(x)

    if len(for_cross) % 2 == 1:
        new_population.append(for_cross.pop())

    random.shuffle(for_cross)

    for protein1, protein2 in zip(for_cross[0:-1:2], for_cross[1::2]):
        r = random.random()
        if r < self._cros_prob:
            while attempts > 0:
                new_protein1, new_protein2 = protein1, protein2
                for gene1, gene2 in zip(new_protein1.genes, new_protein2.genes):
                    r = random.random()
                    assert gene1.cros_prob == gene2.cros_prob
                    if r < gene1.cros_prob and gene1.charged != gene2.charged:
                        gene1.value, gene2.value = gene2.value, gene1.value
                        new_protein1.value = new_protein2.value = None
                if self.is_stable_protein(new_protein1) and self.is_stable_protein(new_protein2):
                    new_population.append(new_protein1)
                    new_population.append(new_protein2)
                elif attempts == 1:
                    new_protein1, new_protein2 = protein1, protein2
                    new_population.append(new_protein1)
                    new_population.append(new_protein2)
                attempts -= 1
    self._population = new_population


# TODO: расширить возможности этапа селекции
def selection(self, eval_param=0.05):
    def distribution(p, m):
        def evaluate(p: float, n: int) -> float:
            return p * pow(1 - p, n - 1)

        vs = []
        for i in range(1, m + 1):
            v = 0
            for j in range(1, i + 1):
                v += evaluate(p, j)
            vs.append(v)
        return vs

    npopulation = []
    pop_size = len(self._population)
    n_bests = 3

    population = sorted(self._population, key=lambda x: x.value, reverse=True)

    bests = []
    for i in range(n_bests):
        x = population[i].copy()
        bests.append(x)

    q = distribution(eval_param, pop_size)
    for _ in range(pop_size):
        n, r = 0, random.uniform(0, q[-1])
        while r > q[n]:
            n += 1
        x = population[n].copy()
        npopulation.append(x)

    npopulation = sorted(npopulation, key=lambda x: x.value, reverse=True)[0:-n_bests]
    npopulation = bests + npopulation

    random.shuffle(npopulation)

    return npopulation


def compute(self):
    for_computing = []
    for protein in self._population:
        sequence = protein.sequence
        if sequence in self._computed:
            protein.value = self._computed[sequence]
        else:
            for_computing.append(protein)

    if for_computing:
        with open(".tempfile", "w") as ouf:
            for protein in for_computing:
                sequence = protein.sequence
                ouf.write(sequence + "\n")
        os.rename(".tempfile", self._output_file)

    while not os.path.exists(self._input_file):
        time.sleep(5)

    with open(self._input_file) as inf:
        for protein in for_computing:
            protein.value = float(inf.readline())
            sequence, value = protein.sequence, protein.value
            self._computed[sequence] = value

    os.remove(self._output_file)
    os.remove(self._input_file)


def save_to_file(self, path):
    with open(path, "a") as ouf:
        for sequence, value in self._computed.items():
            line = f'{sequence}\t{value}\n'
            ouf.write(line)


def load_from_file(self, path):
    with open(path, "r") as inf:
        for line in inf.readlines():
            sequence, value = line.split()
            self._computed[sequence] = value


def is_stable_protein(self, protein):
    if self._checker is not None:
        return self._checker.check(protein)
    return True


def generate_population(pdb_file: str, pop_size: int, computed_path: str) -> List[Protein]:
    population = []

    coordinates = read_coordinates(pdb_file)
    if os.path.exists(computed_path):
        with open(computed_path, "r") as inf:
            line = inf.readline()
            while line and len(population) <= pop_size:
                sequence, value = line.split()
                protein = Protein(sequence)
                protein.value = value
                protein.init_coordinates(coordinates)
                population.append(protein)

    sequence = read_sequence(pdb_file)
    while len(population) <= pop_size:
        protein = Protein(sequence)
        protein.init_coordinates(coordinates)
        population.append(protein)

    return population


def get_best_protein(population: List[Protein]) -> Protein:
    best_protein = max(population, key=lambda x: x.value)
    return best_protein
