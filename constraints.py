from data import Protein
from functools import partial
from numpy.linalg import norm


class Constraints:
    """
    Класс, соедржащий констрэинты
    """

    def __init__(self):
        self._constraints = []

    def add(self, f) -> None:
        """
        Добавляет новый коснтрэинт.
        Все функции должна обладать сигнатурой T -> bool.

        :param f: констрэинт.
        :return:
        """
        self._constraints.append(f)

    def check(self, x) -> bool:
        """
        Проверяет удовлетворение объекта всем констрэинтам

        :param x:
        :return:
        """
        for f in self._constraints:
            if not f(x):
                return False
        return True


def constaint_charge(p: Protein, max_charge: float) -> bool:
    if abs(p.charge) < max_charge:
        return True
    return False


def constraint_n_charged(p: Protein, max_n_charged: float) -> bool:
    n = len(p.genes)
    count = 0
    for i in range(0, n):
        if p[i].charged:
            count += 1
    if count < max_n_charged:
        return True
    return False


def constraint_distances(p: Protein, min_distance: float) -> bool:
    n = len(p.genes)
    for i in range(0, n):
        for j in range(i + 1, n):
            pi, pj = p[i], p[j]
            if pi.charged and pj.charged:
                dist = norm(pi.coordinates - pj.coordinates)
                if dist < min_distance:
                    return False
    return True


def constraint_included(p: Protein, aminoacids_set, positions_set) -> bool:
    for n, gene in enumerate(p.genes):
        if n + 1 in positions_set:
            if gene.value in aminoacids_set:
                return True
    return False


if __name__ == "__main__":
    from evolution import PositionsSet1
    import random
    from data import AMINOACIDS

    constraint = Constraints()

    f1 = partial(constraint_included, aminoacids_set="DE", positions_set=PositionsSet1)
    f2 = partial(constaint_charge, max_charge=7)
    f3 = partial(constraint_n_charged, max_n_charged=60)
    f4 = partial(constraint_distances, min_distance=5.0)

    constraint.add(f1)
    constraint.add(f2)
    constraint.add(f3)
    constraint.add(f4)

    import random
    from data import AMINOACIDS

    sequence = random.choices(AMINOACIDS, k=200)
    protein = Protein(sequence)

    results_f = [f(protein) for f in [f1, f2, f3, f4]]
    print(f"Constraint 1: {results_f[0]}")
    print(f"Constraint 2: {results_f[1]}")
    print(f"Constraint 3: {results_f[2]}")
    print(f"Constraint 4: {results_f[3]}")

    results_checker = constraint.check(protein)
    print(f"Constraints checker: {results_checker}")

    assert all(results_f) == results_checker

    print("OK!")
