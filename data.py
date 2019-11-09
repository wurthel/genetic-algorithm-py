import random

class Gene():
    def __init__(self, position = None, 
                 variants = None, 
                 crosProb = None, mutProb = None):
        self.Position = position
        self.Variants = variants
        self.CrosProb = crosProb
        self.MutProb = mutProb
    
    def SelectAminoacid(self) -> str:
        x = random.choice(self.Variants)
        return x

class Protein():
    def __init__(self, bros):
        self.Variance = dict()
        for gene in bros:
            self.Variance[gene.Position] = gene.SelectAminoacid()
        self.Lambda = None

    def GetSequence(self, patternSequence: str) -> str:
        patternSequence = list(patternSequence)
        for (position, aminoacid) in self.Variance.items():
            patternSequence[position - 1] = aminoacid
        return ''.join(patternSequence)
    
    def GetVariance(self) -> str:
        variance = ''
        for c in self.Variance.values():
            variance += c
        return variance