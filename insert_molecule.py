import numpy as np
import sys
import os
import copy
from scipy.linalg import expm, norm
from math import sqrt

ifnmol = sys.argv[1] # "Poxr.pdb"
ifnzmx = sys.argv[2] # "zmatrix.txt"
ofnmol = sys.argv[3] # "Poxr-result.pdb"

optPath = '/Users/wurthel/Desktop/GitHub/chemistry-scripts/Optimizer/'
optScript = 'bash relax.sh'
optMolBef = '1M0L.pdb'
optMolAft = '1M0L_0.pdb'
optResid = 'N'
optSkip = 4

vdwr = {
    'H': 1.10,
    'O': 1.52,
    'N': 1.55,
    'C': 1.70,
    'S': 1.80
}

class Atom():
    def __init__(self):
        self.DataType = 'ATOM'  		  # "ATOM"/"HETATM"
        self.Name = ''          		  # Atom name
        self.AltLoc = ''        		  # Alternate location indicator. 
        self.ResName = ''       		  # Residue name
        self.ChainId = 'A'       		  # Chain identifier
        self.ResSeq = 0                   # Residue sequence number
        self.ResCode = ''                 # Code for insertions of residues
        self.Coordin = np.array([0,0,0])  # (X,Y,Z) orthogonal Å coordinate
        self.Occup = 0.0        		  # Occupancy
        self.TempFac = 0.0      		  # Temperature factor
        self.Element = 'Xx'      		  # Element symbol
        self.Charge = 0.0         		  # Atom charge

class ZAtom():
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
    return expm(np.cross(np.eye(3), axis/norm(axis)*theta))

def valangle(v1, v2):
    dist01 = np.sqrt(np.dot(v1, v1))
    dist02 = np.sqrt(np.dot(v2, v2))
    dist12 = np.sqrt(np.dot(v2 - v1, v2 - v1))
    angl = (dist01 ** 2 + dist02 ** 2 - dist12 ** 2) / (2 * dist01 * dist02)
    return np.arccos(angl)

def readPdb(filename):
    file = open(filename, 'r')
    Molecule = dict()
    for line in file.readlines():
        DataType = line[0:6].strip()
        if not DataType in ['ATOM', 'HETATM']:
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
        x, y, z = list(map(float, [line[30:38], line[38:46], line[46:54]]))
        atom.Coordin = np.array([x,y,z])
        atom.Occup = 0.0 # float(line[54:60])
        atom.Tempfac = 0.0 # float(line[60:66])
        atom.Element = atom.Name[0] # line[76:78].strip()
        Molecule[num] = atom
    file.close()
    return Molecule

def writePdb(molecule, filename):
     f = open(filename, "w")
     for (idx, atom) in molecule.items():
         line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
                '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}'\
                '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
             a0 = atom.DataType, a1 = idx, a2 = atom.Name, a3 = atom.AltLoc,
             a4 = atom.ResName, a5 = atom.ChainId, a6 = atom.ResSeq, 
             a7 = atom.ResCode, a8 = atom.Coordin, a9 = atom.Occup, 
             a10 = atom.TempFac, a11 = atom.Element, s = ' '
             )
         f.write(line)
     f.close()

def readZMatrix(filename):
    zMatrix = list()
    zFile = open(ifnzmx, 'r')
    zResseq = int(zFile.readline().split()[2])
    zResname = zFile.readline().split()[2]
    for line in zFile.readlines():
        words = line.split()
        zAtom = ZAtom()
        zAtom.Serial = int(words[0])
        zAtom.Con = int(words[1])
        zAtom.Dist = float(words[2])
        zAtom.ValCon = int(words[3])
        zAtom.ValAngle = float(words[4])
        zAtom.DihCon = int(words[5])
        zAtom.DihAngle = float(words[6])
        zAtom.Name = words[7].strip()
        zMatrix.append(zAtom)
    zFile.close()

    return zResseq, zResname, zMatrix

# Remove optFiles
if os.path.exists(optPath + '/' + optMolBef):
    os.remove(optPath  + '/' + optMolBef)
if os.path.exists(optPath  + '/' + optMolAft):
    os.remove(optPath  + '/' + optMolAft)
if os.path.exists(optPath + '/' + optResid):
    os.remove(optPath + optResid)

Molecule = readPdb(ifnmol)
zResseq, zResname, zMatrix = readZMatrix(ifnzmx)

for N in range(len(zMatrix)):
    zAtom = zMatrix[N]

    atom = Atom()
    atom.ResSeq = zResseq
    atom.ResName = zResname
    atom.Name = zAtom.Name
    atom.Element = zAtom.Name[0]

    atom1 = Molecule[zAtom.Con]
    atom2 = Molecule[zAtom.ValCon]
    atom3 = Molecule[zAtom.DihCon]

    oX = np.array([1, 0, 0])
    oY = np.array([0, 1, 0])
    oZ = np.array([0, 0, 1])
    v1 = copy.copy(atom1.Coordin)
    v2 = copy.copy(atom2.Coordin)
    v3 = copy.copy(atom3.Coordin)

    v3 -= v1
    v2 -= v1

    angl1 = valangle(np.array([v2[0], v2[1], 0]), oX)
    if v2[1] > 0: 
        angl1 *= -1
    v2 = np.dot(M(oZ, angl1), v2)
    v3 = np.dot(M(oZ, angl1), v3)

    angl2 = valangle(np.array([v2[0], v2[1], v2[2]]), oX)
    if v2[2] < 0:
        angl2 *= -1
    v2 = np.dot(M(oY, angl2), v2)
    v3 = np.dot(M(oY, angl2), v3)

    angl3 = valangle(np.array([0, v3[1], v3[2]]), oZ)
    if v3[1] < 0:
        angl3 *= -1
    v2 = np.dot(M(oX, angl3), v2)
    v3 = np.dot(M(oX, angl3), v3)

    v = np.array([0,0,0], dtype=float)
    r01 = zAtom.Dist
    r12 = np.sqrt(np.dot(v2, v2))
    v[0] = np.cos(zAtom.ValAngle) * r01 * r12 / v2[0]
    v[2] = np.sqrt(r01 ** 2 - v[0] ** 2)
    v = np.dot(M(oX, -zAtom.DihAngle), v)

    v = np.dot(M(oX, -angl3), v)
    v = np.dot(M(oY, -angl2), v)
    v = np.dot(M(oZ, -angl1), v)
    v += v1

    atom.Coordin = copy.copy(v)

    # Optimization
    resSeqList = set()
    for x in Molecule.values():
        if x.ResSeq == zResseq:
            continue
        if x.ResSeq == zResseq - 1:
            continue
        if x.ResSeq == zResseq + 1:
            continue
        r = atom.Coordin - x.Coordin
        dist = sqrt(np.dot(r, r))
        if dist < (vdwr[atom.Element] + vdwr[x.Element]) * 0.8:
            resSeqList.add(x.ResSeq)

    Molecule[zAtom.Serial] = atom
    if resSeqList and (optSkip <= 0):
        oldDir = os.getcwd()
        os.chdir(optPath)
        writePdb(Molecule, optMolBef)
        with open(optResid,'a') as nFile:
            for n in resSeqList: 
                nFile.write(str(n) + '\n')
        os.system(optScript)
        Molecule = readPdb(optMolAft)
        os.chdir(oldDir)

    optSkip -= 1

writePdb(Molecule, ofnmol)
print ("OKEY-DOKEY!")