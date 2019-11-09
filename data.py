import random

class Gene():
    def __init__(self, position: int = None, variants: str = None, 
                 crosProb: float = None, mutProb: float = None) -> None:
        self.Position = position
        self.Variants = variants
        self.CrosProb = crosProb
        self.MutProb = mutProb
    
    def SelectAminoacid(self) -> str:
        x = random.choice(self.Variants)
        return x

class Protein():
    def __init__(self, genes) -> None:
        self.Variance = dict()
        for gene in genes:
            self.Variance[gene.Position] = gene.SelectAminoacid()
        self.Lambda = None

    def GetSequence(self, patternSequence: str) -> str:
        patternSequence = list(patternSequence)
        for (position, aminoacid) in self.Variance.items():
            patternSequence[position - 1] = aminoacid
        return ''.join(patternSequence)
    
    def UpdateVariance(self, position: int, aminoacid: str) -> None:
        self.Variance[position] = aminoacid
        self.Lambda = None

    def GetVariance(self) -> str:
        return ''.join(self.Variance.values())

    def SetVariance(self, positions: [int], aminoacids: str) -> None:
        for position, aminoacid in zip(positions, aminoacids):
            self.UpdateVariance(position, aminoacid)