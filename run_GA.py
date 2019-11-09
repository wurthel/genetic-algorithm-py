import configparser
from evolution import *

# PARSING CONFIG
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

# COMPUTING
computedProteins = dict()
bestVariance = None
if os.path.exists(computedProteinsPath):
    computedProteins = ReadComputedProteins(computedProteinsPath)
    bestVariance = max(computedProteins, key=computedProteins.get)

population = GeneratePopulation(popSize, genes)
iteration, stopStep = 0, 0
while stopStep < 5:
    iteration += 1
    if iteration > 1:
        population = Selection(population, evalParam)
    Crossover(population, genes, crosProb)
    Mutation(population, genes, mutProb)
    ComputeLambda(population, patternSequence, computedProteins, computeLambdaOuf, computeLambdaInf)
    curBestProtein = GetBestProtein(population)
    if bestVariance == None:
        bestVariance = curBestProtein.GetVariance()
    elif curBestProtein.Lambda <= computedProteins[bestVariance]:
        stopStep += 1
    else:
        bestVariance = curBestProtein.GetVariance()
        stopStep = 0
    SaveComputing(population, computedProteins, computedProteinsPath)
    print(f'Iter: {iteration}. The best result: {bestVariance}, {computedProteins[bestVariance]}')

# WRITING RESULTS
results = sorted(computedProteins.items(), key=lambda x: x[1], reverse=True)
protein = Protein(genes)
with open(resultFileName, 'w') as ouf:
    positions = protein.Variance.keys()
    for variance, lmb in results:
        protein.SetVariance(positions, variance)
        protein.Lambda = lmb
        line = f'{protein.Lambda}, {protein.GetVariance()}\n{protein.GetSequence(patternSequence)}\n\n'
        ouf.write(line)