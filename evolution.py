import random
from copy import copy
import time
import os
from typing import List
from data import Protein
from utils import read_sequence, read_coordinates
from abc import abstractmethod, ABC

PullAPlus = "STNQCWYEDH"
PullBPlus = "PGAVILMFEDH"
PullA = "STNQCWY"
PullB = "PGAVILMF"
PullC = "STNQCGPAVILMFYW"
PullD = "RKHED"

ResiduesSet1 = "STNQCGPAVILMFYWEDH"
ResiduesSet2 = "STNQCGPAVILMFYW"
ResiduesSet3 = "RHKDESTNQCGPAVILMFYW"

PositionsSet1 = [23, 56, 88, 89, 92, 93, 96, 188, 215, 218]
PositionsSet2 = [52, 60, 86, 121, 122, 124, 125, 126, 128, 129, 140, 141, 142, 144, 145, 146, 147, 148, 181, 184, 185,
                 186, 189, 191, 192, 193, 211, 218]
PositionsSet3 = [18, 19, 20, 22, 24, 26, 27, 30, 48, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 80, 81, 82, 84, 85,
                 87, 90, 91, 94, 95, 97, 98, 99, 100, 114, 115, 116, 117, 118, 119, 120, 127, 130, 131, 134, 136, 137,
                 138, 139, 143, 149, 150, 151, 155, 177, 178, 179, 180, 182, 183, 187, 190, 194, 195, 197, 198, 199,
                 207, 208, 209, 210, 212, 213, 214, 216, 217, 220, 221, 222, 224, 225, 226]


class Evolution(ABC):
    def __init__(self):
        self._population = None
        self._mut_prob = None
        self._cros_prob = None

    @property
    def population(self):
        return self._population

    @property
    def mut_prob(self):
        return self._mut_prob

    @property
    def cros_prob(self):
        return self._cros_prob

    @abstractmethod
    def mutation(self):
        pass

    @abstractmethod
    def crossover(self):
        pass

    @abstractmethod
    def selection(self):
        pass


class BaseFunction(ABC):
    def __init__(self):
        self._input_file = None
        self._output_file = None
        self._save_file = None
        self._computed = None

    @abstractmethod
    def compute(self):
        pass

    @abstractmethod
    def save_to_file(self):
        pass

    @abstractmethod
    def load_from_file(self, path):
        pass


class ProteinEvolution(Evolution, BaseFunction):
    def __init__(self, population, mut_prob, cros_prob, input_file, output_file, save_file, checker=None):
        super().__init__()
        self._population = population
        self._mut_prob = mut_prob
        self._cros_prob = cros_prob
        self._input_file = input_file
        self._output_file = output_file
        self._save_file = save_file
        self._checker = checker
        self._computed = dict()

        for protein in self._population:
            if protein.value is not None:
                self._computed[protein.sequence] = protein.value

    def mutation(self, attempts=1):
        def _mutation(protein: Protein):
            for n, gene in enumerate(protein.genes):
                new_value = gene.value
                if n + 1 in PositionsSet1:
                    if gene.polared:
                        if random.random() < 0.2:
                            new_value = random.choice(PullBPlus)
                        elif random.random() < 0.05:
                            new_value = random.choice(PullAPlus)
                    else:
                        if random.random() < 0.2:
                            new_value = random.choice(PullAPlus)
                        elif random.random() < 0.05:
                            new_value = random.choice(PullBPlus)
                    gene.value = new_value
                elif n + 1 in PositionsSet2:
                    if gene.polared:
                        if random.random() < 0.2:
                            new_value = random.choice(PullBPlus)
                        elif random.random() < 0.05:
                            new_value = random.choice(PullA)
                    else:
                        if random.random() < 0.2:
                            new_value = random.choice(PullAPlus)
                        elif random.random() < 0.05:
                            new_value = random.choice(PullB)
                elif n + 1 in PositionsSet3:
                    if gene.charged:
                        if random.random() < 0.2:
                            new_value = random.choice(PullC)
                        elif random.random() < 0.1:
                            new_value = random.choice(PullD)
                    else:
                        if random.random() < 0.2:
                            new_value = random.choice(PullD)
                        elif random.random() < 0.1:
                            new_value = random.choice(PullC)
                gene.value = new_value

        new_population = []
        for protein in self._population:
            if random.random() < self._mut_prob:
                while attempts > 0:
                    new_protein = copy(protein)
                    _mutation(new_protein)
                    if self.is_stable_protein(new_protein):
                        new_population.append(new_protein)
                        break
                    elif attempts == 1:
                        new_protein = copy(protein)
                        new_population.append(new_protein)
                    attempts -= 1

        self._population = new_population

    def crossover(self, attempts=1):
        new_population = []
        for_cross = []

        for x in self._population:
            if random.random() < self._cros_prob:
                for_cross.append(x)
            else:
                new_population.append(x)

        if len(for_cross) % 2 == 1:
            new_population.append(for_cross.pop())

        random.shuffle(for_cross)

        for protein1, protein2 in zip(for_cross[0:-1:2], for_cross[1::2]):
            while attempts > 0:
                new_protein1, new_protein2 = copy(protein1), copy(protein2)
                for gene1, gene2 in zip(new_protein1.genes, new_protein2.genes):
                    if random.random() < self._cros_prob:
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
    def selection(self, eval_param=0.05, save_n_best=3):
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

        population = sorted(self._population, key=lambda x: x.value, reverse=True)

        bests = []
        for i in range(save_n_best):
            x = population[i].copy()
            bests.append(x)

        q = distribution(eval_param, pop_size)
        for _ in range(pop_size):
            n, r = 0, random.uniform(0, q[-1])
            while r > q[n]:
                n += 1
            x = population[n].copy()
            npopulation.append(x)

        npopulation = sorted(npopulation, key=lambda x: x.value, reverse=True)[0:-save_n_best]
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

    def save_to_file(self):
        with open(self._save_file, "a") as ouf:
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

    def get_best_protein(self) -> Protein:
        best_protein = max(self.population, key=lambda x: x.value)
        return best_protein


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
