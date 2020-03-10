import numpy as np
import sys
import os

from copy import copy
from pathlib import Path
from utils import *

ifnmol = sys.argv[1]  # "tests/D95E.pdb"
ifnzmx = sys.argv[2]  # "tests/zmatrix.txt"
ofnmol = sys.argv[3]  # "tests/D95E_result.pdb"

opt_path = Path("/home/nikolaev_d/vusal-dir/genetic-algorithm-py/Optimizer")
opt_cmd = "bash relax.sh"
opt_skip = 4


# Remove opt files
opt_mol_bef = "1M0L.pdb"
opt_mol_aft = "1M0L_0.pdb"
opt_resids = 'N'
if os.path.isfile(opt_path / opt_mol_bef):
    os.remove(opt_path / opt_mol_bef)
if os.path.isfile(opt_path / opt_mol_aft):
    os.remove(opt_path / opt_mol_aft)
if os.path.isfile(opt_path / opt_resids):
    os.remove(opt_path / opt_resids)

molecule = read_pdb(ifnmol)
zresseq, zresname, zmatriz = read_zmatrix(ifnzmx)
resSeqList = set()

for N in range(len(zmatriz)):
    zatom = zmatriz[N]

    atom = Atom()
    atom.ResSeq = zresseq
    atom.ResName = zresname
    atom.Name = zatom.Name
    atom.Element = zatom.Name[0]

    atom1 = molecule[zatom.Con]
    atom2 = molecule[zatom.ValCon]
    atom3 = molecule[zatom.DihCon]

    oX = np.array([1.0, 0.0, 0.0])
    oY = np.array([0.0, 1.0, 0.0])
    oZ = np.array([0.0, 0.0, 1.0])
    v1 = copy(atom1.Coordin)
    v2 = copy(atom2.Coordin)
    v3 = copy(atom3.Coordin)

    v3 -= v1
    v2 -= v1

    angle1 = val_angle(np.array([v2[0], v2[1], 0]), oX)
    if v2[1] > 0:
        angle1 *= -1
    v2 = np.dot(M(oZ, angle1), v2)
    v3 = np.dot(M(oZ, angle1), v3)

    angle2 = val_angle(np.array([v2[0], v2[1], v2[2]]), oX)
    if v2[2] < 0:
        angle2 *= -1
    v2 = np.dot(M(oY, angle2), v2)
    v3 = np.dot(M(oY, angle2), v3)

    angle3 = val_angle(np.array([0, v3[1], v3[2]]), oZ)
    if v3[1] < 0:
        angle3 *= -1
    v2 = np.dot(M(oX, angle3), v2)
    v3 = np.dot(M(oX, angle3), v3)

    v = np.array([0.0, 0.0, 0.0])
    r01 = zatom.Dist
    r12 = np.linalg.norm(v2)
    v[0] = np.cos(zatom.ValAngle) * r01 * r12 / v2[0]
    v[2] = np.sqrt(r01 ** 2 - v[0] ** 2)
    v = np.dot(M(oX, -zatom.DihAngle), v)

    v = np.dot(M(oX, -angle3), v)
    v = np.dot(M(oY, -angle2), v)
    v = np.dot(M(oZ, -angle1), v)
    v += v1

    atom.Coordin = copy(v)

    # Optimization
    for x in molecule.values():
        if x.ResSeq == zresseq:
            continue
        if x.ResSeq == zresseq - 1:
            continue
        if x.ResSeq == zresseq + 1:
            continue
        r = atom.Coordin - x.Coordin
        dist = np.linalg.norm(r)
        if dist < (vdwr[atom.Element] + vdwr[x.Element]) * 0.8:
            resSeqList.add(x.ResSeq)

    molecule[zatom.Serial] = atom
    if resSeqList and not (opt_skip > 0):
        old_dir = os.getcwd()
        os.chdir(opt_path)
        write_pdb(molecule, opt_mol_bef)
        with open(opt_resids, "w") as f:
            for n in resSeqList:
                f.write(f"{n}\n")
        os.system(opt_cmd)
        molecule = read_pdb(opt_mol_aft)
        os.chdir(old_dir)
    elif opt_skip > 0:
        opt_skip -= 1

write_pdb(molecule, ofnmol)
