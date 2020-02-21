import openbabel
from sys import argv
import numpy as np

ifname = argv[1]  # "1U19.pdb"
ofname = argv[2]  # "1U19-zmatrix.txt"

izresseq = int(argv[3])  # 296
zresnm = argv[4]  # "REL"
ozresseq = int(argv[5])  # 297

startnum = 8000
zenum = {
    -3: ("N", izresseq - 1, None),
    -2: ("CA", izresseq - 1, None),
    -1: ("C", izresseq - 1, None),
    0: ("N", izresseq, [-1, -2, -3]),
    1: ("CA", izresseq, [0, -1, -2]),
    2: ("C", izresseq, [1, 0, -1]),
    3: ("O", izresseq, [2, 1, 0]),
    4: ("CB", izresseq, [1, 2, 3]),
    5: ("CG", izresseq, [4, 1, 2]),
    6: ("CD", izresseq, [5, 4, 1]),
    7: ("CE", izresseq, [6, 5, 4]),
    8: ("NZ", izresseq, [7, 6, 5]),
    9: ("C15", izresseq, [8, 7, 6]),
    10: ("C14", izresseq, [9, 8, 7]),
    11: ("C13", izresseq, [10, 9, 8]),
    12: ("C20", izresseq, [11, 10, 9]),
    13: ("C12", izresseq, [11, 10, 9]),
    14: ("C11", izresseq, [13, 11, 10]),
    15: ("C10", izresseq, [14, 13, 11]),
    16: ("C9", izresseq, [15, 14, 13]),
    17: ("C19", izresseq, [16, 15, 14]),
    18: ("C8", izresseq, [16, 15, 14]),
    19: ("C7", izresseq, [18, 16, 15]),
    20: ("C6", izresseq, [19, 18, 16]),
    21: ("C5", izresseq, [20, 19, 18]),
    22: ("C18", izresseq, [21, 20, 19]),
    23: ("C4", izresseq, [21, 20, 19]),
    24: ("C3", izresseq, [23, 21, 20]),
    25: ("C2", izresseq, [24, 23, 21]),
    26: ("C1", izresseq, [25, 24, 23]),
    27: ("C17", izresseq, [26, 25, 24]),
    28: ("C16", izresseq, [26, 25, 24]),
}


def GetAtomByName(idx):
    aname = zenum[idx][0]
    residue = mol.GetResidue(zenum[idx][1] - 1)
    for obatom in openbabel.OBResidueAtomIter(residue):
        name = residue.GetAtomID(obatom).strip()
        if name == aname:
            return obatom
    assert False, "Atom not found"


obConversion = openbabel.OBConversion()
obConversion.SetInFormat('pdb')
mol = openbabel.OBMol()
obConversion.ReadFile(mol, ifname)

with open(ofname, "w") as outfile:
    outfile.write(f'RESSEQ = {str(ozresseq)}\n')
    outfile.write(f'RESNAME = {zresnm}\n')

    for i in range(len(zenum) - 3):
        a0 = GetAtomByName(i)
        a1 = GetAtomByName(zenum[i][2][0])
        a2 = GetAtomByName(zenum[i][2][1])
        a3 = GetAtomByName(zenum[i][2][2])

        # Distance
        v0 = np.array([a0.x(), a0.y(), a0.z()])
        v1 = np.array([a1.x(), a1.y(), a1.z()])
        dist = np.sqrt(np.dot(v0 - v1, v0 - v1))

        # ValAngl
        v2 = np.array([a2.x(), a2.y(), a2.z()])
        dist01 = dist
        dist12 = np.sqrt(np.dot(v1 - v2, v1 - v2))
        dist02 = np.sqrt(np.dot(v0 - v2, v0 - v2))
        valangl = (dist01 ** 2 + dist12 ** 2 - dist02 ** 2) / (2 * dist01 * dist12)
        valangl = np.arccos(valangl)

        # DihAngl
        v3 = np.array([a3.x(), a3.y(), a3.z()])
        b0 = -1.0 * (v1 - v0)
        b1 = v2 - v1
        b2 = v3 - v2
        b1 /= np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        dihangl = np.arctan2(y, x)

        # Numeration
        serial = startnum + i
        con = a1.GetIdx()
        valcon = a2.GetIdx()
        dihcon = a3.GetIdx()
        if i > 0:
            con = startnum + zenum[i][2][0]
        if i > 1:
            valcon = startnum + zenum[i][2][1]
        if i > 2:
            dihcon = startnum + zenum[i][2][2]

        # Write to
        line = f"{serial}\t{con}\t{dist}\t{valcon}\t{valangl}\t" \
               f"{dihcon}\t{dihangl}\t{zenum[i][0]}\n"
        outfile.write(line)

print("OKEY-DOKEY!")
