import random
import copy
import time
import os
import shutil
from typing import List, Dict
from data import Gene, Protein

def read_genes(config) -> List[Gene]:
    genes = list()
    for (position, params) in config.items():
        gene = Gene()
        params = params.split(' ')
        gene.position = int(position)
        gene.variants = params[0]
        gene.cros_prob = float(params[1])
        gene.mut_prob  = float(params[2])
        genes.append(gene)
    return genes

def generate_population(pop_size: int, genes: List[Gene],
                       compute_lmb_ouf: str, compute_lmb_inf: str) -> List[Protein]:
    proteins = list()
    for _ in range(pop_size):
        protein = Protein(genes)
        proteins.append(protein)

    if os.path.exists(compute_lmb_ouf) and os.path.exists(compute_lmb_inf):
        lmb_ouf_file = open(compute_lmb_ouf, 'r')
        lmb_inf_file = open(compute_lmb_inf, 'r')
        lmb_ouf_lines = list(filter(None, [line.strip() for line in lmb_ouf_file]))
        lmb_inf_lines = list(filter(None, [line.strip() for line in lmb_inf_file]))
        
        for (protein, l1, l2) in zip(proteins, lmb_ouf_file, lmb_inf_file):
            protein.update_variance_from_sequence(genes, l1)
            protein.mlambda = float(l2)

        lmb_ouf_file.close()
        lmb_inf_file.close()

        os.remove(compute_lmb_ouf)
        os.remove(compute_lmb_inf)

    return proteins

def crossover(population: List[Protein], genes: List[Gene], cros_prob: float) -> None:
    proteins_for_cross = list()
    for protein in population:
        r = random.random()
        if r < cros_prob:
            proteins_for_cross.append(protein)

    if len(proteins_for_cross) % 2 == 1:
        proteins_for_cross = proteins_for_cross[0:-1]

    random.shuffle(proteins_for_cross)

    for a, b in zip(proteins_for_cross[0:-1:2], proteins_for_cross[1::2]):
        for gene in genes:
            r = random.random()
            if r < gene.cros_prob:
                pos = gene.position
                x, y = a.variance[pos], b.variance[pos]
                a.update_variance(pos, y)
                b.update_variance(pos, x)

def mutation(population: List[Protein], genes: List[Gene], mut_prob: float) -> None:
    for protein in population:
        r = random.random()
        if r < mut_prob:
            for gene in genes:
                r = random.random()
                if r < gene.mut_prob:
                    protein.update_variance(gene.position, gene.select_aminoacid())

def evaluate(eval_param: float, n: int) -> float:
    return eval_param * pow(1 - eval_param, n - 1)

def selection(population: List[Protein], eval_param: float) -> (List[Protein], List[Protein]):
    for protein in population:
        if protein.mlambda == None:
            print('WARNING: selection is missing. The population contains uncomputed lambda.')
            return [], population

    best_proteins, other_proteins = list(), list()
    population = sorted(population, key=lambda protein: protein.mlambda, reverse=True)
    for _ in range(3):
        protein = population.pop(0)
        protein = copy.deepcopy(protein)
        best_proteins.append(protein)

    pop_size = len(population)
    q = list(map(lambda n: sum(map(lambda m: evaluate(eval_param, m), range(1, n + 1))), range(1, pop_size + 1)))
    for _ in range(pop_size):
        n, r = 0, random.uniform(0, q[-1])
        while r > q[n]:
            n += 1
        protein = copy.deepcopy(population[n])
        other_proteins.append(protein)
    random.shuffle(other_proteins)

    return best_proteins, other_proteins

def compute_lambda(population: List[Protein], pattern_seq: str,
                  computed_proteins: Dict[str, float],
                  compute_lmb_ouf: str, compute_lmb_inf: str) -> None:
    proteins_for_computing = list()
    for protein in population:
        variance = protein.get_variance()
        if variance in computed_proteins:
            protein.mlambda = computed_proteins[variance]
        else:
            proteins_for_computing.append(protein)

    if proteins_for_computing:
        with open('.tempfile', 'w') as ouf:
            for protein in proteins_for_computing:
                seq = protein.get_sequence(pattern_seq)
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
            variance = protein.get_variance()
            lmb = protein.mlambda
            if variance in computed_proteins or lmb == None:
                continue
            computed_proteins[variance] = lmb
            line = f'{variance}\t{lmb}\n'
            ouf.write(line)

def read_computed_proteins(path: str) -> Dict[str, float]:
    computed_proteins = dict()
    with open(path, 'r') as inf:
        for line in inf.readlines():
            variance, value = line.split()
            computed_proteins[variance] = float(value)
    return computed_proteins

def get_best_protein(population: List[Protein]) -> Protein:
    bestProtein = max(population, key=lambda p: p.mlambda)
    return bestProtein