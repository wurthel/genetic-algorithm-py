import random
import copy
import time
import os
import shutil
from typing import List, Dict
from data import Gene, Protein, CHARGED, NON_CHARGED
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


def generate_population(sequence: str, pop_size: int, compute_lmb_ouf: str, compute_lmb_inf: str) -> List[Protein]:
    proteins = list()
    for _ in range(pop_size):
        protein = Protein(sequence)
        proteins.append(protein)

    if os.path.exists(compute_lmb_ouf) and os.path.exists(compute_lmb_inf):
        lmb_ouf_file = open(compute_lmb_ouf, "r")
        lmb_inf_file = open(compute_lmb_inf, "r")

        for (protein, seq, val) in zip(proteins, lmb_ouf_file, lmb_inf_file):
            protein.update_genes_from_sequence(seq)
            protein.value = float(val)

        lmb_ouf_file.close()
        lmb_inf_file.close()

        # os.remove(compute_lmb_ouf)
        # os.remove(compute_lmb_inf)

    return proteins


def crossover(population: List[Protein], cros_prob: float) -> None:
    proteins_for_cross = list()
    for protein in population:
        r = random.random()
        if r < cros_prob:
            proteins_for_cross.append(protein)

    if len(proteins_for_cross) % 2 == 1:
        proteins_for_cross = proteins_for_cross[0:-1]

    random.shuffle(proteins_for_cross)

    for p1, p2 in zip(proteins_for_cross[0:-1:2], proteins_for_cross[1::2]):
        for g1, g2 in zip(p1.genes, p2.genes):
            assert g1.cros_prob == g2.cros_prob and g1.position == g2.position
            r = random.random()
            if r < g1.cros_prob and g1.type != g2.type:
                g2.value, g1.value = g1.value, g2.value
                p1.value, p2.value = None, None


def mutation(population: List[Protein], genes: List[Gene], mut_prob: float) -> None:
    for protein in population:
        r = random.random()
        if r < mut_prob:
            for gene in protein.genes:
                r = random.random()
                if r < gene.mut_prob:
                    gene.value = random.choice(NON_CHARGED) if gene.type == "CHARGED" else random.choice(CHARGED)
                    protein.value = None


def evaluate(eval_param: float, n: int) -> float:
    return eval_param * pow(1 - eval_param, n - 1)


def selection(population: List[Protein], eval_param: float) -> (List[Protein], List[Protein]):
    for protein in population:
        if protein.value is None:
            print('WARNING: selection is missing. The population contains uncomputed lambda.')
            return [], population

    best_proteins, new_population = list(), list()
    population = sorted(population, key=lambda x: x.value, reverse=True)
    for _ in range(3):
        protein = population.pop(0)
        best_proteins.append(copy.deepcopy(protein))
        new_population.append(copy.deepcopy(protein))

    pop_size = len(population)
    q = list(map(lambda n: sum(map(lambda m: evaluate(eval_param, m), range(1, n + 1))), range(1, pop_size + 1)))
    for _ in range(pop_size):
        n, r = 0, random.uniform(0, q[-1])
        while r > q[n]:
            n += 1
        protein = population[n]
        new_population.append(copy.deepcopy(protein))
    random.shuffle(new_population)

    return best_proteins, new_population


def compute_lambda(population: List[Protein], pattern_seq: str,
                   computed_proteins: Dict[str, float],
                   compute_lmb_ouf: str, compute_lmb_inf: str) -> None:
    proteins_for_computing = list()
    for protein in population:
        variance = ''.join(protein.variance.values())
        if variance in computed_proteins:
            value = computed_proteins[variance]
            protein.update_value(value)
        else:
            proteins_for_computing.append(protein)

    if proteins_for_computing:
        with open('.tempfile', 'w') as ouf:
            for protein in proteins_for_computing:
                seq = protein.sequence
                ouf.write(seq + '\n')
        os.rename('.tempfile', compute_lmb_ouf)

        while not os.path.exists(compute_lmb_inf):
            time.sleep(5)

        with open(compute_lmb_inf, 'r') as inf:
            for protein in proteins_for_computing:
                protein.mlambda = float(inf.readline())

        os.remove(compute_lmb_ouf)
        os.remove(compute_lmb_inf)


def save_computing(population: List[Protein], computed_proteins: Dict[str, float], path: str) -> None:
    with open(path, 'a') as ouf:
        for protein in population:
            sequence = protein.sequence
            value = protein.value
            if sequence in computed_proteins or value is None:
                continue
            computed_proteins[sequence] = value
            line = f'{sequence}\t{value}\n'
            ouf.write(line)


def read_computed_proteins(path: str) -> Dict[str, float]:
    computed_proteins = dict()
    with open(path, "r") as inf:
        for line in inf.readlines():
            sequence, value = line.split()
            computed_proteins[sequence] = float(value)
    return computed_proteins


def get_best_protein(population: List[Protein]) -> Protein:
    best_protein = max(population, key=lambda x: x.value)
    return best_protein
