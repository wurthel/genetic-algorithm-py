import random
import numpy as np
import copy
from data import Gene, Protein

TMP_SEQUENCE =  "MLMTVFSSAPELALLGSTFAQVDPSNLSVSDSLTYGQFNL" \
                "VYNAFSFAIAAMFASALFFFSAQALVGQRYRLALLVSAIV" \
                "VSIAGYHYFRIFNSWDAAYVLENGVYSLTSEKFNDAYRYV" \
                "DWLLTVPLLLVETVAVLTLPAKEARPLLIKLTVASVLMIA" \
                "GYPGEISDDITTRIIWGTVSTIPFAYILYVLWVELSRSLV" \
                "RQPAAVQTLVRNMRWLLLLSWGVYPIAYLLPMLGVSGTSA" \
                "AVGVQVGYTIADVLAKPVFGLLVFAIALVKTKADQESSEP" \
                "HAAIGAAANKSGGSLIS"
TMP_VARIANCE =  [(1,"D"), (2,"T"), (3,"L"), (4,"W")]
i = 0

def ReadBros(path: str):
    bros = list()
    with open(path) as file:
        for line in file.readlines():
            x = line.split(' ')
            gene = Gene()
            gene.Position = int(x[0])
            gene.Variants = x[1]
            gene.CrosProb = float(x[2])
            gene.MutProb  = float(x[3])
            bros.append(gene)
    return bros

def GeneratePopulation(n: int, bros: [Gene]) -> [Protein]:
    proteins = list()
    for _ in range(n):
        proteins.append(GenerateProtein(bros))
    return proteins

def GenerateProtein(bros: [Gene]) -> Protein:
    global i
    protein = Protein(TMP_SEQUENCE, TMP_VARIANCE, i)
    i+=1
    protein.GenerateRandomVariance(bros)
    return protein

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
        a.UpdateSequence()
        b.UpdateSequence()

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
        protein.UpdateSequence()

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


bros = ReadBros('brosVariance')
population = GeneratePopulation(2, bros)
Crossover(population, bros, 1)
Mutation(population, bros, 1)
population = Selection(population, 0.05)
population = sorted(population, key=lambda protein: protein.Lambda, reverse=True)
for p in population:
    print(p.GetSequence())
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