#!/usr/bin/env python
#
#Dictionary for atom types
#
ATOM_TYPE_DICT = {
        'ASP': {"N":"NH1",
                "H":"H",
                "CA":"CT3", #CT1->CT3 because of cap  
                "HA":"HA3", #HA1->HA3 because of cap
                "HA2":"HA3", #HA1->HA3 because of cap
                "HA3":"HA3", #HA1->HA3 because of cap
                "CB":"CT2",                                        
                "HB2":"HA2", #HA->HA2 because of connection to CT2
                "HB3":"HA2", #HA->HA2 because of connection to CT2
                "CG":"CC",
                "OD1":"OC",
                "OD2":"OC",
                "C":"C",
                "O":"O"},
        'CYM': {"N":"NH1",
                "H":"H",
                "CA":"CT3", 
                "HA":"HA3", 
                "HA2":"HA3",
                "HA3":"HA3",
                "CB":"CT2",
                "HB2":"HA2",
                "HB3":"HA2",
                "SG":"SCD",
                "C":"C",
                "O":"O"},
        'GLU': {"N":"NH1",
                "H":"H",
                "CA":"CT3", #CT1->CT3 because of cap  
                "HA":"HA3", #HA1->HA3 because of cap
                "HA2":"HA3", #HA1->HA3 because of cap
                "HA3":"HA3", #HA1->HA3 because of cap
                "CB":"CT2",                                        
                "HB2":"HA2", #HA->HA2 because of connection to CT2
                "HB3":"HA2", #HA->HA2 because of connection to CT2
                "CG":"CT2",
                "HG2":"HA2", #HA->HA2 because of connection to CT2
                "HG3":"HA2", #HA->HA2 because of connection to CT2
                "CD":"CC",
                "OE1":"OC",
                "OE2":"OC",
                "C":"C",
                "O":"O"},
        'HID': {"N":"NH1",
                "H":"H",
                "CA":"CT3", #CT1->CT3 because of cap  
                "HA":"HA3", #HA1->HA3 because of cap
                "HA2":"HA3", #HA1->HA3 because of cap
                "HA3":"HA3", #HA1->HA3 because of cap
                "CB":"CT2",                                        
                "HB2":"HA2", #HA->HA2 because of connection to CT2
                "HB3":"HA2", #HA->HA2 because of connection to CT2
                "CG":"CPHG",
                "ND1":"NR1",
                "HD1":"H",
                "CE1":"CPH2",
                "HE1":"HR1",
                "NE2":"NR2",
                "CD2":"CPHD", 
                "HD2":"HR3",
                "C":"C",
                "O":"O"},
        'HIE': {"N":"NH1",
                "H":"H",
                "CA":"CT3", #CT1->CT3 because of cap  
                "HA":"HA3", #HA1->HA3 because of cap
                "HA2":"HA3", #HA1->HA3 because of cap
                "HA3":"HA3", #HA1->HA3 because of cap
                "CB":"CT2",                                       
                "HB2":"HA2", #HA->HA2 because of connection to CT2
                "HB3":"HA2", #HA->HA2 because of connection to CT2
                "CG":"CPHG",
                "ND1":"NR2",
                "CE1":"CPH2", 
                "HE1":"HR1",
                "NE2":"NR1",
                "HE2":"H",
                "CD2":"CPHD",
                "HD2":"HR3",
                "C":"C",
                "O":"O"},
        'SER': {"N":"NH1",
                "H":"H",
                "CA":"CT3", #CT1->CT3 because of cap  
                "HA":"HA3", #HA1->HA3 because of cap
                "HA2":"HA3", #HA1->HA3 because of cap
                "HA3":"HA3", #HA1->HA3 because of cap
                "CB":"CSD",
                "HB2":"HSD", #HA->HA2 because of connection to CT2
                "HB3":"HSD", #HA->HA2 because of connection to CT2
                "OG":"OH1",
                "HG":"H",
                "C":"C",
                "O":"O"},
        'SED': {"N":"NH1",
                "H":"H",
                "CA":"CT3", 
                "HA":"HA3", 
                "HA2":"HA3",
                "HA3":"HA3",
                "CB":"CSD", 
                "HB2":"HSD",
                "HB3":"HSD",
                "OG":"OH1",
                "C":"C",
                "O":"O"},
        'ZN': {"ZN":"ZN"
                },
        'WAT': {"OA":"OW",
                "HA":"HW",
                "HA2":"HW"
                }
    }

#
# Following are parameters
#

# Covalent radii
co_radii_dict = {
    "CT3":  0.68,
    "HA3":  0.23,
    "CT2":  0.68,
    "HA2":  0.23,
    "SCD":  1.02,
    "CPHG": 0.68,
    "NR1":  0.68,
    "H":    0.23,
    "CPH2": 0.68,
    "HR1":  0.23,
    "NR2":  0.68,
    "CPHD": 0.68,
    "HR3":  0.23,
    "ZN":   1.45,
    "CC":   0.68,
    "OC":   0.68,
    "OW":   0.68,
    "HW":   0.23,
    "CSD":  0.68,
    "HSD":  0.23,
    "OH1":  0.68,
    }

# CM5 EQeq parameters
cm5eqeq_para = {
    "CT3":  [319.65,240.34],
    "HA3":  [232.61,490.24],
    "CT2":  [303.04,240.52],
    "HA2":  [228.16,491.67],
    "SCD":  [369.20,257.31],
    "CPHG": [252.63,238.97],
    "NR1":  [328.11,175.89],
    "H":    [112.17,467.59],
    "CPH2": [237.60,232.62],
    "HR1":  [197.32,513.96],
    "NR2":  [364.87,268.18],
    "CPHD": [272.58,249.65],
    "HR3":  [209.09,514.29],
    "ZN":   [158.55,192.15],
    "CC":   [238.30,184.66],
    "OC":   [369.53,273.77],
    "OW":   [430.30,311.98],
    "HW":   [106.65,464.26],
    "CSD":  [282.57,213.84],
    "HSD":  [235.83,502.95],
    "OH1":  [409.06,314.70],
    }

#CM5 Tkk parameters
cm5tkk_para = {
        ('ZN', 'SCD'): -0.0345,
        ('ZN', 'NR2'): 0.0159,
        ('ZN', 'OW'): 0.0216,
        ('ZN', 'OC'): 0.0083,
        ('ZN', 'OH1'): 0.0077
        }

