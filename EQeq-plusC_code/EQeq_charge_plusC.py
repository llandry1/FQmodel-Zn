#!/usr/bin/env python
from readpdb import get_atominfo_fpdb
from data_para_dict import *
from numpy import zeros, ones, sqrt, exp, vstack
from numpy.linalg import solve as matrix_solver

# Unit conversion factors
kcal2kj = 4.184
ev2kcal = 23.0609
ev2kj = ev2kcal * kcal2kj
K = 1389.0 # Unit kJ/mol*A
K = K/ev2kj # Unit: eV*A

# Define an atom
class fq_Atom:
    def __init__(self, resid_atname, atm_typ, crd, ele_neg, J):
        self.resid_atname = resid_atname
        # resid_atname was used instead of atname is to differenciate different
        # atoms with the same atom name but within diferent residues
        self.atm_typ = atm_typ
        self.crd = crd
        self.ele_neg = ele_neg
        self.J = J

# Calculate the distance
def cal_dis(crd1, crd2):
    if len(crd1) != len(crd2):
        raise ValueError('The two coordinates should have the same number of elements')

    result = 0.0
    for i in range (0, len(crd1)):
        result += (crd1[i]-crd2[i])**2
    return sqrt(result)

# This is the EQeq hardness equation
# Eq 64 in the SI of the JPCL, 2012, 3, 2506-2511
def cal_Jij(J_i, J_j, r, die=1.0):
    J_ij = sqrt(J_i * J_j)
    Eo = exp(-(J_ij*r/K)**2) * (J_ij / K - (J_ij**2) * r / (K**2) - 1.0 / r)
    return (K / 2.0 / die) * (1.0 / r + Eo)

def cal2_Jij(r):
    Jij = 1.0 / 7.0 * 0.5 * K / r
    return Jij

###############################################################################
# This is the main program to calculate the Qeq charges
# The default method for calculating the hardness is the simple method
###############################################################################
def get_EQeq_pC_chg(pdbf, qtot, fq_para_dict, Tkk_para, hardness_equ='simple', aa_dict=ATOM_TYPE_DICT):

    mol, atids, resids = get_atominfo_fpdb(pdbf)

    # The total charge of the molecule
    q_total= float(qtot)

    # Create a molecule with all the atoms and their parameters in it
    Atoms = {}

    for i in atids:
        resname_i =  mol.atoms[i].resname
        atname_i = mol.atoms[i].atname
        resid_i = mol.atoms[i].resid
        resid_atname = str(resid_i) + str("-") + str(atname_i)
        crd_i = mol.atoms[i].crd

        if resname_i not in aa_dict.keys():
            raise ValueError("Could not find the resname %s in the "
                             "CHARMM amino acid dictionary" %resname_i)
        else:
            res_dict = aa_dict[resname_i]
            if atname_i not in res_dict:
                raise ValueError("Could not find the atom name %s in the "
                                 "residue %s in the CHARMM amino acid dictionary" %atname_i)
            else:
                atyp_i = res_dict[atname_i]
                if atyp_i not in fq_para_dict.keys():
                    raise ValueError("Could not find the atom type %s in the "
                                     "CHARMM FQ parameter dictionary" %atyp_i)
                else:
                    ele_neg_i = fq_para_dict[atyp_i][0]/ev2kcal
                    J_i = fq_para_dict[atyp_i][1]/ev2kcal


        Atoms[i] = fq_Atom(resid_atname, atyp_i, crd_i, ele_neg_i, J_i)
    
    # Perform the calculations in the matrix
    num_atoms = len(Atoms)
    J_matrix = zeros([num_atoms, num_atoms])
    a1 = zeros([num_atoms-1, num_atoms])
    a2 = zeros([num_atoms-1, num_atoms])
    b = zeros(num_atoms)
    b[0] = q_total

    for i in range(1, num_atoms+1):
        for j in range(i, num_atoms+1):
            if (i == j):
                Jij = Atoms[i].J
            else:
                crdi = Atoms[i].crd
                crdj = Atoms[j].crd
                J_i = Atoms[i].J
                J_j = Atoms[j].J
                disij = cal_dis(crdi, crdj)
                if hardness_equ == 'eqeq':
                    Jij = cal_Jij(J_i, J_j, disij)
                elif hardness_equ == 'simple':
                    Jij = cal2_Jij(disij)
            J_matrix[i-1, j-1] = Jij
            J_matrix[j-1, i-1] = Jij

    a1 = J_matrix[0:num_atoms-1]
    a2 = J_matrix[1:num_atoms]
    a = a1 - a2

    one_array = ones(num_atoms)
    a = vstack((one_array, a))

    for i in range(1, num_atoms):
        j = i + 1
        ele_negi = Atoms[i].ele_neg
        ele_negj = Atoms[j].ele_neg
        b[i] = ele_negj - ele_negi

    # Solve the main matrix, herein the FQchg_l is the list which the Eqeq
    # without the +C terms
    FQchg_l = matrix_solver(a, b)
    FQchg_l = [round(i, 4) for i in FQchg_l]

    # Get the charges due to the +C term
    plusC_l = get_plusC(Atoms, Tkk_para)

    # Get the final Eqeq+C charges
    FQchg_l_plusC_l = []

    for (x, y) in zip(FQchg_l, plusC_l):
        FQchg_l_plusC_l.append(x+y)

    FQchg_l_plusC_l = [round(i, 4) for i in FQchg_l_plusC_l]

    # Store the final Eqeq+C charges in an dictionary
    results_dict = {}

    #print("The charge of the atoms are:")
    for i in range(len(atids)):
        j = atids[i]
        #print(FQchg_l_plusC_l[i])
        results_dict[Atoms[j].resid_atname] = [Atoms[j].atm_typ, FQchg_l_plusC_l[i]]

        # IF YOU ONLY WANT SPECIFIC ATOM TYPES PRINTED
        #if [Atoms[j].atm_typ] == ['SCD']:
        #    print(FQchg_l_plusC_l[i])
        #if [Atoms[j].atm_typ] == ['NR2']:
        #    print(FQchg_l_plusC_l[i])
        #if [Atoms[j].atm_typ] == ['ZN']:
        #   print(FQchg_l_plusC_l[i])

    return results_dict

# This function is to get the +C term, the Tkk_para should follow a certain
# sequence based on the system one is working on
def get_plusC(Atoms, Tkk_para):

    num_atoms = len(Atoms)

    # A list to store the charges
    plusC_l = []

    # To loop over all the atoms
    for i in range(1, num_atoms+1):
        plusC_value = 0
        atyp_i = Atoms[i].atm_typ
        Ri = co_radii_dict[atyp_i]

        # This list is just for debugging purposes
        Tkk_l = []

        # To loop over all the atoms' influence on atom i
        for j in range(1, num_atoms+1):
            atyp_j = Atoms[j].atm_typ
            Rj = co_radii_dict[atyp_j]
            sumR = Ri + Rj

            Tkk_value = 0

            if (atyp_i, atyp_j) in Tkk_para.keys():
                Tkk_value = Tkk_para[(atyp_i, atyp_j)]
            elif (atyp_j, atyp_i) in Tkk_para.keys():
                Tkk_value = -Tkk_para[(atyp_j, atyp_i)]

            # The alpha value is from JCTC, 2012, 8, 527-541
            alpha = 2.474
                
            crdi = Atoms[i].crd
            crdj = Atoms[j].crd
            disij=cal_dis(crdi, crdj)
          
            Bkk_value = exp(-alpha * (disij - sumR))
            plusC_value += Bkk_value * Tkk_value

        plusC_l.append(plusC_value)

    return(plusC_l)

