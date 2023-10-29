#!/usr/bin/env python
from EQeq_charge_plusC import get_EQeq_pC_chg
from data_para_dict import *

# Contains residue name and atom name
pdbf = 'ZnWAT6.pdb'

# Atom type dict: residue name, atname <-> atom type
at_dict = ATOM_TYPE_DICT

# Atom types <-> electronegativity, chemical hardness
fq_para_dict = cm5eqeq_para

# Tkk parameters
Tkk_para = cm5tkk_para

# Total charges is +2.0, which is important
qtot = 2.0

# Calculate the charges using the optimized parameters
chg_l = get_EQeq_pC_chg(pdbf, qtot, fq_para_dict, Tkk_para, hardness_equ='charmm', aa_dict=at_dict)
print(chg_l)

quit()

