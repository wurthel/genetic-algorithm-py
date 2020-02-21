import numpy as np
import sys
import os

from copy import copy
from pathlib import Path
from scipy.linalg import expm, norm
from math import sqrt

ifnmol = sys.argv[1]  # "Poxr.pdb"
ifnzmx = sys.argv[2]  # "zmatrix.txt"
ofnmol = sys.argv[3]  # "Poxr-result.pdb"

opt_path = Path("/Users/wurthel/Desktop/GitHub/chemistry-scripts/Optimizer/")
opt_cmd = "bash relax.sh"
opt_mol_bef = "1M0L.pdb"
opt_mol_aft = "1M0L_0.pdb"
opt_resid = 'N'
opt_skip = 4

vdwr = {
    'H': 1.10,
    'O': 1.52,
    'N': 1.55,
    'C': 1.70,
    'S': 1.80
}


class Atom:
    def __init__(self):
        self.DataType = 'ATOM'  # "ATOM"/"HETATM"
        self.Name = ''  # Atom name
        self.AltLoc = ''  # Alternate location indicator.
        self.ResName = ''  # Residue name
        self.ChainId = 'A'  # Chain identifier
        self.ResSeq = 0  # Residue sequence number
        self.ResCode = ''  # Code for insertions of residues
        self.Coordin = np.array([0, 0, 0])  # (X,Y,Z) orthogonal Å coordinate
        self.Occup = 0.0  # Occupancy
        self.TempFac = 0.0  # Temperature factor
        self.Element = 'Xx'  # Element symbol
        self.Charge = 0.0  # Atom charge


class ZAtom:
    def __init__(self):
        self.Serial = 0
        self.Con = 0
        self.Dist = 0.0
        self.ValCon = 0
        self.ValAngle = 0.0
        self.DihCon = 0
        self.DihAngle = 0.0
        self.Name = 'Xx'


def M(axis, theta):
    return expm(np.cross(np.eye(3), axis / norm(axis) * theta))


def val_angle(vec1, vec2):
    dist01 = np.linalg.norm(vec1)
    dist02 = np.linalg.norm(vec2)
    dist12 = np.linalg.norm(vec1 - vec2)
    angle = (dist01 ** 2 + dist02 ** 2 - dist12 ** 2) / (2 * dist01 * dist02)
    return np.arccos(angle)


def read_pdb(filename):
    with open(filename, 'r') as file:
        molecule = dict()
        for line in file.readlines():
            data_type = line[0:6].strip()
            if data_type not in ['ATOM', 'HETATM']:
                continue
            atom = Atom()
            atom.DataType = line[0:6].strip()
            num = int(line[6:11])
            atom.Name = line[12:16].strip()
            atom.AltLoc = line[16].strip()
            atom.ResName = line[17:20].strip()
            atom.ChainId = line[21].strip()
            atom.ResSeq = int(line[22:26])
            atom.ResCode = line[26].strip()
            atom.Coordin = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
            atom.Occup = 0.0  # float(line[54:60])
            atom.Tempfac = 0.0  # float(line[60:66])
            atom.Element = atom.Name[0]  # line[76:78].strip()

            molecule[num] = atom

    return molecule


def write_pdb(molecule, filename):
    f = open(filename, "w")
    for (idx, atom) in molecule.items():
        line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
               '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}' \
               '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
            a0=atom.DataType, a1=idx, a2=atom.Name, a3=atom.AltLoc,
            a4=atom.ResName, a5=atom.ChainId, a6=atom.ResSeq,
            a7=atom.ResCode, a8=atom.Coordin, a9=atom.Occup,
            a10=atom.TempFac, a11=atom.Element, s=' '
        )
        f.write(line)
    f.close()


def read_zmatrix(filename):
    z_matrix = list()
    with open(ifnzmx, 'r') as file:
        resseq = int(file.readline().split()[2])
        resname = file.readline().split()[2]
        for line in file.readlines():
            words = line.split()
            atom = ZAtom()
            atom.Serial = int(words[0])
            atom.Con = int(words[1])
            atom.Dist = float(words[2])
            atom.ValCon = int(words[3])
            atom.ValAngle = float(words[4])
            atom.DihCon = int(words[5])
            atom.DihAngle = float(words[6])
            atom.Name = words[7].strip()

            z_matrix.append(atom)

    return resseq, resname, z_matrix


# Remove optFiles
if os.path.isfile(opt_path / opt_mol_bef):
    os.remove(opt_path / opt_mol_bef)
if os.path.isfile(opt_path / opt_mol_aft):
    os.remove(opt_path / opt_mol_aft)
if os.path.isfile(opt_path / opt_resid):
    os.remove(opt_path / opt_resid)

molecule = read_pdb(ifnmol)
z_resseq, z_resname, z_matrix = read_zmatrix(ifnzmx)

for N in range(len(z_matrix)):
    z_atom = z_matrix[N]

    atom = Atom()
    atom.ResSeq = z_resseq
    atom.ResName = z_resname
    atom.Name = z_atom.Name
    atom.Element = z_atom.Name[0]

    atom1 = molecule[z_atom.Con]
    atom2 = molecule[z_atom.ValCon]
    atom3 = molecule[z_atom.DihCon]

    oX = np.array([1, 0, 0])
    oY = np.array([0, 1, 0])
    oZ = np.array([0, 0, 1])
    v1 = copy(atom1.Coordin)
    v2 = copy(atom2.Coordin)
    v3 = copy(atom3.Coordin)

    v3 -= v1
    v2 -= v1

    angl1 = val_angle(np.array([v2[0], v2[1], 0]), oX)
    if v2[1] > 0:
        angl1 *= -1
    v2 = np.dot(M(oZ, angl1), v2)
    v3 = np.dot(M(oZ, angl1), v3)

    angl2 = val_angle(np.array([v2[0], v2[1], v2[2]]), oX)
    if v2[2] < 0:
        angl2 *= -1
    v2 = np.dot(M(oY, angl2), v2)
    v3 = np.dot(M(oY, angl2), v3)

    angl3 = val_angle(np.array([0, v3[1], v3[2]]), oZ)
    if v3[1] < 0:
        angl3 *= -1
    v2 = np.dot(M(oX, angl3), v2)
    v3 = np.dot(M(oX, angl3), v3)

    v = np.array([0, 0, 0], dtype=float)
    r01 = z_atom.Dist
    r12 = np.sqrt(np.dot(v2, v2))
    v[0] = np.cos(z_atom.ValAngle) * r01 * r12 / v2[0]
    v[2] = np.sqrt(r01 ** 2 - v[0] ** 2)
    v = np.dot(M(oX, -z_atom.DihAngle), v)

    v = np.dot(M(oX, -angl3), v)
    v = np.dot(M(oY, -angl2), v)
    v = np.dot(M(oZ, -angl1), v)
    v += v1

    atom.Coordin = copy(v)

    # Optimization
    resSeqList = set()
    for x in molecule.values():
        if x.ResSeq == z_resseq:
            continue
        if x.ResSeq == z_resseq - 1:
            continue
        if x.ResSeq == z_resseq + 1:
            continue
        r = atom.Coordin - x.Coordin
        dist = sqrt(np.dot(r, r))
        if dist < (vdwr[atom.Element] + vdwr[x.Element]) * 0.8:
            resSeqList.add(x.ResSeq)

    molecule[z_atom.Serial] = atom
    if resSeqList and (opt_skip <= 0):
        old_dir = os.getcwd()
        os.chdir(opt_path)
        write_pdb(molecule, opt_mol_bef)
        with open(opt_resid, 'a') as nFile:
            for n in resSeqList:
                nFile.write(str(n) + '\n')
        os.system(opt_cmd)
        Molecule = read_pdb(opt_mol_aft)
        os.chdir(old_dir)

    opt_skip -= 1

write_pdb(molecule, ofnmol)
print("OKEY-DOKEY!")
