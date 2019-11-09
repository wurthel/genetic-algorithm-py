import random
import copy
import time
import configparser
import os
import shutil
from typing import List, Dict
from data import Gene, Protein

def ReadGenes(config) -> List[Gene]:
    genes = list()
    for (position, params) in config.items():
        gene = Gene()
        params = params.split(' ')
        gene.Position = int(position)
        gene.Variants = params[0]
        gene.CrosProb = float(params[1])
        gene.MutProb  = float(params[2])
        genes.append(gene)
    return genes

def GeneratePopulation(popSize: int, genes: List[Gene]) -> List[Protein]:
    proteins = list()
    for _ in range(popSize):
        protein = Protein(genes)
        proteins.append(protein)
    return proteins

def Crossover(population: List[Protein], genes: List[Gene], crosProb: float) -> None:
    proteinsForCross = list()
    for protein in population:
        r = random.random()
        if r > crosProb:
            continue
        proteinsForCross.append(protein)

    if len(proteinsForCross) % 2 == 1:
        proteinsForCross = proteinsForCross[0:-1]

    random.shuffle(proteinsForCross)
    
    for a, b in zip(proteinsForCross[0:-1:2], proteinsForCross[1::2]):
        for gene in genes:
            r = random.random()
            if r > gene.CrosProb:
                continue
            pos = gene.Position
            x, y = a.Variance[pos], b.Variance[pos]
            a.UpdateVariance(pos, y)
            b.UpdateVariance(pos, x)

def Mutation(population: List[Protein], genes: List[Gene], mutProb: float) -> None:
    for protein in population:
        r = random.random()
        if r > mutProb:
            continue
        for gene in genes:
            r = random.random()
            if r > gene.MutProb:
                continue
            protein.UpdateVariance(gene.Position, gene.SelectAminoacid())

def Eval(evalParam: float, n: int) -> float:
    return evalParam * pow(1 - evalParam, n - 1)

def Selection(population: List[Protein], evalParam: float) -> List[Protein]:
    newPopulation = list()
    popSize = len(population)
    
    population = sorted(population, key=lambda protein: protein.Lambda, reverse=True)
    q = list(map(lambda n: sum(map(lambda m: Eval(evalParam, m), range(1, n + 1))), range(1, popSize + 1)))
    for _ in range(popSize):
        n = 0
        r = random.uniform(0, q[-1])
        while r > q[n]:
            n += 1
        newPopulation.append(copy.deepcopy(population[n]))
    
    return newPopulation

def ComputeLambda(population: List[Protein], patternSeq: str, 
                  computedProteins: Dict[str, float],
                  computeLambdaOuf: str, computeLambdaInf: str) -> None:
    needComputing = False
    for protein in population:
        variance = protein.GetVariance()
        if variance in computedProteins:
            protein.Lambda = computedProteins[variance]
            continue
        needComputing = True

    if needComputing:
        with open('.tempfile', 'w') as ouf:
            for protein in population:
                if protein.Lambda != None:
                    continue
                seq = protein.GetSequence(patternSeq)
                ouf.write(seq + '\n\n')
        os.rename('.tempfile', computeLambdaOuf)
        
        while not os.path.exists(computeLambdaInf):
            time.sleep(1)

        with open(computeLambdaInf, 'r') as inf:
            for protein in population:
                if protein.Lambda != None:
                    continue
                protein.Lambda = float(inf.readline())

        os.remove(computeLambdaOuf)
        os.remove(computeLambdaInf)
        
def SaveComputing(population: List[Protein], computedProteins: Dict[str, float], path: str) -> None:
    with open(path, 'a') as ouf:
        for protein in population:
            variance = protein.GetVariance()
            if variance in computedProteins:
                continue
            computedProteins[variance] = protein.Lambda
            line = f'{variance}\t{protein.Lambda}\n'
            ouf.write(line)

def ReadComputedProteins(path: str) -> Dict[str, float]:
    computedProteins = dict()
    with open(path, 'r') as inf:
        for line in inf.readlines():
            variance, value = line.split()
            computedProteins[variance] = float(value)
    return computedProteins

def GetBestProtein(population: List[Protein]) -> Protein:
    bestProtein = max(population, key=lambda p: p.Lambda)
    return bestProtein

config = configparser.ConfigParser()
config.read('config.ini')

genes = ReadGenes(config['GENES'])
patternSequence = config['TEMPLATE']['Sequence']

crosProb = float(config['PARAMS']['CrosProb'])
mutProb = float(config['PARAMS']['MutProb'])
evalParam = float(config['PARAMS']['EvalParam'])
popSize = int(config['PARAMS']['PopSize'])

computeLambdaInf = config['COMPUTING']['ComputeLambdaInf']
computeLambdaOuf = config['COMPUTING']['ComputeLambdaOuf']
computedProteinsPath = config['COMPUTING']['ComputedProteinsFileName']
resultFileName = config['COMPUTING']['ResultFileName']

computedProteins = dict()
if os.path.exists(computedProteinsPath):
    computedProteins = ReadComputedProteins(computedProteinsPath)

bestProtein, iteration, stopStep = None, 1, 1
population = GeneratePopulation(popSize, genes)
ComputeLambda(population, patternSequence, computedProteins, computeLambdaOuf, computeLambdaInf)
SaveComputing(population, computedProteins, computedProteinsPath)
bestProtein = copy.deepcopy(GetBestProtein(population))
print(f'Iter: {iteration}. The best lambda: {bestProtein.Lambda}')
while stopStep < 5:
    iteration += 1
    population = Selection(population, evalParam)
    Crossover(population, genes, crosProb)
    Mutation(population, genes, mutProb)
    ComputeLambda(population, patternSequence, computedProteins, computeLambdaOuf, computeLambdaInf)
    SaveComputing(population, computedProteins, computedProteinsPath)

    curBestProtein = GetBestProtein(population)
    if curBestProtein.Lambda <= bestProtein.Lambda:
        stopStep += 1
    else:
        bestProtein = copy.deepcopy(curBestProtein)
        stopStep = 0
    print(f'Iter: {iteration}. The best lambda: {bestProtein.Lambda}')

population = sorted(population, key=lambda protein: protein.Lambda, reverse=True)
with open(resultFileName, 'w') as ouf:
    for protein in population:
        line = f'{protein.Lambda}\n{protein.GetSequence(patternSequence)}\n\n'
        ouf.write(line)