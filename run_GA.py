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

population = GeneratePopulation(popSize, genes, computeLambdaOuf, computeLambdaInf)
SaveComputing(population, computedProteins, computedProteinsPath)

bestVariance = None
if os.path.exists(computedProteinsPath):
    computedProteins = ReadComputedProteins(computedProteinsPath)
    if computedProteins:
        bestVariance = max(computedProteins, key=computedProteins.get)

iteration, step, stopStep = 1, 0, 5
while step < stopStep:
    population = Selection(population, evalParam)
    Crossover(population, genes, crosProb)
    Mutation(population, genes, mutProb)
    ComputeLambda(population, patternSequence, computedProteins, computeLambdaOuf, computeLambdaInf)
    curBestProtein = GetBestProtein(population)
    if bestVariance == None:
        bestVariance = curBestProtein.GetVariance()
    elif curBestProtein.Lambda <= computedProteins[bestVariance]:
        step += 1
    else:
        bestVariance = curBestProtein.GetVariance()
        step = 0
    SaveComputing(population, computedProteins, computedProteinsPath)
    line =  f'Iter: {iteration}. The best result: {bestVariance}, {computedProteins[bestVariance]}'
    line += f' | {step}/{stopStep}'
    print(line)
    iteration += 1

# WRITING RESULTS
results = sorted(computedProteins.items(), key=lambda x: x[1], reverse=True)
protein = Protein(genes)
with open(resultFileName, 'w') as ouf:
    for variance, lmb in results:
        protein.SetVariance(variance)
        line = f'{lmb}, {variance}\n{protein.GetSequence(patternSequence)}\n\n'
        ouf.write(line)