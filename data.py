import random
from typing import List

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
        self.Lambda = None
        for gene in genes:
            self.Variance[gene.Position] = gene.SelectAminoacid()

    def UpdateVarianceFromSequence(self, genes: List[Gene], sequence: str) -> None:
        self.Variance = dict()
        self.Lambda = None
        for gene in genes:
            pos = gene.Position
            self.Variance[pos] = sequence[pos - 1]

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

    def SetVariance(self, aminoacids: str) -> None:
        positions = self.Variance.keys()
        for position, aminoacid in zip(positions, aminoacids):
            self.UpdateVariance(position, aminoacid)