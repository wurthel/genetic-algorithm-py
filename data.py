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
    def __init__(self, sequence: str, variance: str, lmd: float):
        self.Sequence = list(sequence)
        self.Variance = dict(variance)
        self.Lambda = lmd

    def GenerateRandomVariance(self, bros: [Gene]):
        self.Variance = dict()
        for b in bros:
            x = b.SelectAminoacid()
            self.Variance[b.Position] = x
        self.UpdateSequence()

    def UpdateSequence(self):
        for (position, aminoacid) in self.Variance.items():
            self.Sequence[position - 1] = aminoacid

    def GetSequence(self) -> str:
        sequence = ''
        for c in self.Sequence:
            sequence += c
        return sequence
    
    def GetVariance(self) -> str:
        variance = ''
        for c in self.Variance.values():
            variance += c
        return variance