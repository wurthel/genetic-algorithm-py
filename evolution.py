import random
from copy import copy
import time
import os
from typing import List
from data import Protein
from utils import read_sequence, read_coordinates
from abc import abstractmethod, ABC
from data import Gene

PullAPlus = "STNQCWYEDH"
PullBPlus = "PGAVILMFEDH"
PullA = "STNQCWY"
PullB = "PGAVILMF"
PullC = "STNQCGPAVILMFYW"
PullD = "RKHED"

ResiduesSet1 = "STNQCGPAVILMFYWEDH"
ResiduesSet2 = "STNQCGPAVILMFYW"
ResiduesSet3 = "RHKDESTNQCGPAVILMFYW"

PositionsSet1 = {23, 56, 88, 89, 92, 93, 96, 188, 215}
PositionsSet2 = {52, 60, 86, 121, 122, 124, 125, 126, 128, 129, 140, 141, 142, 144, 145, 146, 147, 148, 181, 184, 185,
                 186, 189, 191, 192, 193, 211, 218}
PositionsSet3 = {18, 19, 20, 22, 24, 26, 27, 30, 48, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 80, 81, 82, 84, 85,
                 87, 90, 91, 94, 95, 97, 98, 99, 100, 114, 115, 116, 117, 118, 119, 120, 127, 130, 131, 134, 136, 137,
                 138, 139, 143, 149, 150, 151, 155, 177, 178, 179, 180, 182, 183, 187, 190, 194, 195, 197, 198, 199,
                 207, 208, 209, 210, 212, 213, 214, 216, 217, 220, 221, 222, 224, 225, 226}
PositionsSetUnion = set.union(PositionsSet1, PositionsSet2, PositionsSet3)


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
        self._saved = None

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
        self._saved = set()

    def mutation(self, attempts=1):
        """

        :param attempts: число попыток инциниализации на один protein
        :return:
        """

        def _mutation(protein: Protein):
            for i, gene in enumerate(protein.genes):
                new_value = gene.value
                if i + 1 in PositionsSet1:
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
                elif i + 1 in PositionsSet2:
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
                elif i + 1 in PositionsSet3:
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
                new_gene = Gene(value=new_value)
                protein.update_gene(i, new_gene)

        # must_be = 0
        # real = 0

        new_population = []
        for protein in self._population:
            new_protein = copy(protein)
            if random.random() < self._mut_prob:
                # must_be += 1
                _attempts = attempts
                while _attempts > 0:
                    attempt_protein = copy(new_protein)
                    _mutation(attempt_protein)
                    if self.is_stable_protein(attempt_protein):
                        new_protein = attempt_protein
                        # real += 1
                        break
                    _attempts -= 1
            new_population.append(new_protein)

        # print(f"Mutation: must_be/real {must_be}/{real}")

        self._population = new_population

    def crossover(self, attempts=1):
        pair_cros_prob = 0.5  # вероятность обмена аминокислотами
        new_population = []
        for_cross = []

        for protein in self._population:
            new_protein = copy(protein)
            if random.random() < self._cros_prob:
                for_cross.append(new_protein)
            else:
                new_population.append(new_protein)

        if len(for_cross) % 2 == 1:
            new_population.append(for_cross.pop())

        random.shuffle(for_cross)

        # must_be = 0
        # real = 0

        for protein1, protein2 in zip(for_cross[0:-1:2], for_cross[1::2]):
            # must_be += 2
            _attempts = attempts
            new_protein1, new_protein2 = protein1, protein2
            while _attempts > 0:
                attempt_protein1, attempt_protein2 = copy(new_protein1), copy(new_protein2)
                for i, (gene1, gene2) in enumerate(zip(attempt_protein1.genes, attempt_protein2.genes)):
                    if random.random() < pair_cros_prob:
                        new_gene1 = Gene(value=gene2.value)
                        new_gene2 = Gene(value=gene1.value)
                        attempt_protein1.update_gene(i, new_gene1)
                        attempt_protein2.update_gene(i, new_gene2)
                if self.is_stable_protein(attempt_protein1) and self.is_stable_protein(attempt_protein2):
                    new_protein1 = attempt_protein1
                    new_protein2 = attempt_protein2
                    # real += 2
                    break
                _attempts -= 1
            new_population.append(new_protein1)
            new_population.append(new_protein2)

        # print(f"Crossover: must_be/real {must_be}/{real}")

        self._population = new_population

    # TODO: расширить возможности этапа селекции
    def selection(self, eval_param, save_n_best):
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
            protein = copy(population[i])
            bests.append(protein)

        q = distribution(eval_param, pop_size)
        for _ in range(pop_size):
            n, r = 0, random.uniform(0, q[-1])
            while r > q[n]:
                n += 1
            protein = copy(population[n])
            npopulation.append(protein)

        npopulation = sorted(npopulation, key=lambda x: x.value)[0:-save_n_best]
        npopulation.extend(bests)

        random.shuffle(npopulation)

        self._population = npopulation

    def compute(self):
        for_computing = []

        # already_computed = 0
        for protein in self._population:
            sequence = protein.sequence
            if sequence not in self._computed:
                for_computing.append(sequence)
                self._computed[sequence] = None
            # else:
            #     already_computed += 1
        # print(f"Already computed: {already_computed}")

        with open(".tempfile", "w") as ouf:
            for sequence in for_computing:
                ouf.write(sequence + "\n")
        os.rename(".tempfile", self._output_file)

        # print(f"For computing: {len(for_computing)}")
        while not os.path.exists(self._input_file):
            time.sleep(5)

        # computed = 0
        with open(self._input_file) as inf:
            for sequence in for_computing:
                value = float(inf.readline())
                self._computed[sequence] = value
                # computed += 1
        # print(f"Computed: {computed}")

        for protein in self._population:
            sequence = protein.sequence
            value = self._computed[sequence]
            protein.set_value(value)

        os.remove(self._output_file)
        os.remove(self._input_file)

    def save_to_file(self):
        with open(self._save_file, "a") as ouf:
            for sequence, value in self._computed.items():
                if sequence not in self._saved:
                    line = f'{sequence}\t{value}\n'
                    ouf.write(line)
                    self._saved.add(sequence)

    def load_from_file(self, path):
        with open(path, "r") as inf:
            for line in inf.readlines():
                sequence, value = line.split()
                self._computed[sequence] = value
                self._saved.add(sequence)

    def is_stable_protein(self, protein):
        if self._checker is not None:
            return self._checker.check(protein)
        return True

    def get_best_protein(self) -> Protein:
        best_protein = max(self.population, key=lambda x: x.value)
        return best_protein

    def generate_population(self, default_sequence, pop_size):
        population = []

        if os.path.exists(self._save_file):
            with open(self._save_file, "r") as inf:
                line = inf.readline()
                while line:
                    sequence, value = line.split()
                    value = float(value)
                    protein = Protein(sequence)
                    protein.set_value(value)

                    population.append(protein)

                    self._saved.add(sequence)
                    self._computed[sequence] = value

                    line = inf.readline()
                population = sorted(population, key=lambda x: x.value, reverse=True)[:pop_size]

        while len(population) < pop_size:
            protein = Protein(default_sequence)
            population.append(protein)

        self._population = population

    def print_current_population(self):
        for protein in self._population:
            print(protein.sequence, protein.value)


def get_best_protein(population: List[Protein]) -> Protein:
    best_protein = max(population, key=lambda x: x.value)
    return best_protein
