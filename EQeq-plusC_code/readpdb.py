#!/usr/bin/env python
"""
This is the code for reading and writting pdb files. This file is adpated from the
pyMSMT module in the AmberTools software package.
"""
from __future__ import print_function

class Atom:
    def __init__(self, gtype, atid, atname, element, atomtype, crd, charge, resid, resname):
        self.gtype = gtype
        self.atid = atid
        self.atname = atname
        self.element = element
        self.atomtype = atomtype
        self.crd = crd
        self.charge = charge
        self.resid = resid
        self.resname = resname

class Residue:
    def __init__(self, resid, resname, resconter):
        self.resid = resid
        self.resname = resname
        self.resconter = resconter

class Molecule:
    def __init__(self, atoms, residues):
        self.atoms = atoms
        self.residues = residues

Metal1pdb = {('F', 'F'): 'F', ('BR', 'BR'): 'Br', ('CL', 'CL'): 'Cl', 
             ('IOD', 'I'): 'I', ('LI', 'LI'): 'Li', ('NA', 'NA'): 'Na',
             ('RB', 'RB'): 'Rb', ('TL', 'TL'): 'Tl', ('CS', 'CS'): 'Cs',
             ('K', 'K'): 'K', ('CU1', 'CU'): 'Cu', ('AG', 'AG'): 'Ag',
             ('AU', 'AU'): 'Au'}

Metal2pdb = {('CU', 'CU'): 'Cu', ('NI', 'NI'): 'Ni', ('PT', 'PT'): 'Pt',
             ('ZN', 'ZN'): 'Zn', ('CO', 'CO'): 'Co', ('PD', 'PD'): 'Pd',
             ('FE2', 'FE2'): 'Fe', ('MG', 'MG'): 'Mg', ('MN', 'MN'): 'Mn',
             ('HG', 'HG'): 'Hg', ('CD', 'CD'): 'Cd', ('YB2', 'YB2'): 'Yb',
             ('CA', 'CA'): 'Ca', ('PB', 'PB'): 'Pb', ('EU', 'EU'): 'Eu',
             ('SR', 'SR'): 'Sr', ('BA', 'BA'): 'Ba'}

Metal3pdb = {('AL', 'AL'): 'Al', ('FE', 'FE'): 'Fe', ('CR', 'CR'): 'Cr',
             ('IN', 'IN'): 'In', ('Y', 'Y'): 'Y', ('LA', 'LA'): 'La',
             ('CE', 'CE'): 'Ce', ('PR', 'PR'): 'Pr', ('SM', 'SM'): 'Sm',
             ('EU3', 'EU'): 'Eu', ('GD3', 'GD'): 'Gd', ('TB', 'TB'): 'Tb',
             ('LU', 'LU'): 'Lu', ('V', 'V'): 'V', ('ARS', 'AS'): 'As',
             ('RU', 'RU'): 'Ru'}

Metal4pdb = {('IR', 'IR'): 'Ir', ('MO', 'MO'): 'Mo'}

Metalpdb = Metal1pdb
Metalpdb.update(Metal2pdb)
Metalpdb.update(Metal3pdb)
Metalpdb.update(Metal4pdb)

METAL_PDB = Metalpdb

class FileError(Exception):
    def __init__(self, info='File Error'):
        self.info = info
    def __str__(self):
        return self.info

def get_atominfo_fpdb(fname):
    Atoms = {}
    Residues = {}

    atids = []
    resids = []
    resnamedict = {}
    conterdict = {}

    fp = open(fname, 'r')

    for line in fp:
        if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
            gtype = line[0:6].strip(" ")
            atid = int(line[6:11])
            atids.append(atid)
            atname = line[12:16].strip(" ")
            allocind = line[16:17]
            resname = line[17:20].strip(" ")
            chainid = line[21:22]
            resid = int(line[22:26])
            codeinsert = line[26:27]
            crdx = float(line[30:38])
            crdy = float(line[38:46])
            crdz = float(line[46:54])
            crd = (crdx, crdy, crdz)
            occp = line[54:60]
            tempfac = line[60:66]
            atomtype = line[76:78].strip(" ")
            charge = line[78:80]

            if (resname, atname) in list(METAL_PDB.keys()):
                element = METAL_PDB[(resname, atname)][0]
            elif atname[0:2].upper() in ['CL', 'BR']:
                element = atname[0].upper() + atname[1].lower()
            else:
                element = atname[0]

            if atid not in list(Atoms.keys()):
                Atoms[atid] = Atom(gtype, atid, atname, element, atomtype, crd, charge, resid, resname)
            else:
                raise FileError('There are more than one atom with atom id '
                                  '%d in the PDB file : %s .' %(atid, fname))

            if resid not in resids:
                resids.append(resid)
            if resid not in list(resnamedict.keys()):
                resnamedict[resid] = resname

    fp.close()

    resids.sort()

    for i in resids:
        preconter = []
        for j in atids:
            if (Atoms[j].resid == i) and (j not in preconter):
                preconter.append(j)
        preconter.sort()
        conterdict[i] = preconter

    for i in resids:
        resname = resnamedict[i]
        resconter = conterdict[i]
        Residues[i] = Residue(i, resname, resconter)

    del resnamedict
    del conterdict

    mol = Molecule(Atoms, Residues)

    return mol, atids, resids

