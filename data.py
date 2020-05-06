import numpy as np

AMINOACIDS = "RHKDESTNQCGPAVILMFYW"
CHARGED = "RHKDE"
NON_CHARGED = "STNQCGPAVILMFYW"
POSITIVE_CHARGE = "RHK"
NEGATIVE_CHARGE = "DE"
POLAR = "STCNQKRYDEH"
NONPOLAR = "GAVILPMFW"


class Gene:
    def __init__(self, value) -> None:
        self.__value: str = value

    @property
    def charged(self):
        if self.value in CHARGED:
            return True
        return False

    @property
    def polared(self):
        if self.value in POLAR:
            return True
        return False

    @property
    def value(self):
        return self.__value

    def __copy__(self):
        copy = Gene(self.__value)
        return copy


class Protein:
    def __init__(self, sequence, value=None) -> None:
        self.__value = value
        self.__genes = [Gene(value=x) for x in sequence]
        self.__origin_sequence = sequence
        self.__num_changes = 0

    def __getitem__(self, item):
        return self.__genes[item]

    @property
    def charge(self):
        charge = 0
        for v in self.sequence:
            if v in POSITIVE_CHARGE:
                charge += 1
            if v in NEGATIVE_CHARGE:
                charge -= 1
        return charge

    @property
    def sequence(self):
        return ''.join([str(x.value) for x in self.__genes])

    @property
    def genes(self):
        return self.__genes

    @property
    def value(self):
        return self.__value

    @property
    def num_changes(self):
        return self.__num_changes

    def set_value(self, new_value):
        self.__value = new_value

    def update_gene(self, idx, gene):
        """

        :param idx: позиция гена в sequence. 0 <= idx < len(sequqence)
        :param gene: новый ген
        :return:
        """
        if gene.value != self.__genes[idx].value:
            if self.__genes[idx].value != self.__origin_sequence[idx]:
                self.__num_changes -= 1
            self.__genes[idx] = gene
            self.__value = None
            if self.__genes[idx].value != self.__origin_sequence[idx]:
                self.__num_changes += 1

    def __copy__(self):
        copy = Protein(self.sequence, self.value)
        copy._Protein__num_changes = self.__num_changes
        return copy
