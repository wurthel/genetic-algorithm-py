import random
import numpy as np
import copy
import configparser
from data import Gene, Protein

def ReadBros(config) -> [Gene]:
    bros = list()
    for (position, params) in config.items():
        gene = Gene()
        params = params.split(' ')
        gene.Position = position
        gene.Variants = params[0]
        gene.CrosProb = float(params[1])
        gene.MutProb  = float(params[2])
        bros.append(gene)
    return bros

def GeneratePopulation(n: int, bros: [Gene]) -> [Protein]:
    proteins = list()
    for _ in range(n):
        protein = Protein(bros)
        proteins.append(protein)
    return proteins

def Crossover(population: [Protein], bros: [Gene], crosProb: float):
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
        for gene in bros:
            r = random.random()
            if r > gene.CrosProb:
                continue
            pos = gene.Position
            a.Variance[pos], b.Variance[pos] = b.Variance[pos], a.Variance[pos]

def Mutation(population: [Protein], bros: [Gene], mutProb: float):
    proteinsForMut = list()
    for protein in population:
        r = random.random()
        if r > mutProb:
            continue
        proteinsForMut.append(protein)

    for protein in proteinsForMut:
        for gene in bros:
            r = random.random()
            if r > gene.MutProb:
                continue
            pos = gene.Position
            protein.Variance[pos] = gene.SelectAminoacid()

def Eval(evalParam: float, n: int) -> float:
    return evalParam * pow(1 - evalParam, n - 1)

def Selection(population: [Protein], evalParam: float) -> [Protein]:
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

def ComputeLambda(population: [Protein]) -> Protein:
    for protein in population:
        protein.Lambda = random.randint(400,500)

config = configparser.ConfigParser()
config.read('config.ini')

bros = ReadBros(config['BROS'])
crosProb = float(config['PARAMS']['CrosProb'])
mutProb = float(config['PARAMS']['MutProb'])
evalParam = float(config['PARAMS']['EvalParam'])
popSize = int(config['PARAMS']['PopSize'])

population = GeneratePopulation(popSize, bros)
ComputeLambda(population)
for _ in range(10):
    Crossover(population, bros, crosProb)
    Mutation(population, bros, mutProb)
    ComputeLambda(population)
    population = Selection(population, evalParam)

population = sorted(population, key=lambda protein: protein.Lambda, reverse=True)
for p in population:
    # print(p.GetSequence(TMP_SEQUENCE))
    print(p.GetVariance())
    print(p.Lambda)
    print()



# -- | НАЧАЛО. ЧТЕНИЕ И ВЫВОД.
# computeLambda :: [Protein] -> [Protein] -> IO [Protein]
# computeLambda pop all = do
#     -- |Реализовать алгоритм выбора белков @p_a@ и @p_b@
#     -- можно иначе -- сначала найти @p_b@, а потом искать @p_a@
#     -- как дополнение @p_b@ до @p@. Такая реализация
#     -- окажется быстрее. СДЕЛАТЬ ПОТОМ
#     let pa = [x | x <- pop , x `notElem` all]
#         pb = map (\x -> x {lambda = lambda (inPs x)}) (pop \\ pa)
#          where inPs x = fromJust $ find (== x) all 
    
#     print "Current population"
#     print $ zip (map variance pop) (map lambda pop)     
#     print "Will be computed"
#     print $ zip (map variance pa) (map lambda pa)
#     print "Won't be computed"
#     print $ zip (map variance pb) (map lambda pb)
#     print "All proteins"   
#     print $ zip (map variance all) (map lambda all)
#     print "----------------------------------------------"

#     (tmpNameOuf, tmpHandleOuf) <- openTempFile "." "temp"
#     mapM_ (hPutStrLn tmpHandleOuf . protein) pa
#     hClose tmpHandleOuf
#     renameFile tmpNameOuf computeLambdaOuf
#     wait True
#     handleInf <- openFile computeLambdaInf ReadMode
#     pa' <- mapM (\p' -> hGetLine handleInf >>= return . (\x -> p' { lambda = Just x}) . read) pa
#     hClose handleInf
#     removeFile computeLambdaOuf
#     removeFile computeLambdaInf

#     return (pa' <> pb) 
#     where 
#         wait False = return ()
#         wait True  = do 
#             e <- doesFileExist computeLambdaInf 
#             if e then wait False
#             else threadDelay timeWait >> wait True

# writeInProteinFile :: String -> [Protein] -> IO [Protein]
# writeInProteinFile msg ps = withFile resultFile AppendMode write >> return ps
#             where write hdl = hPutStrLn hdl msg >> mapM_ (hPrint hdl) ps
# -- | КОНЕЦ. ЧТЕНИЕ И ВЫВОД.