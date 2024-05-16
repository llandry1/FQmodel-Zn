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
        'CR': {"CR":"CR"
                },
        'MN': {"MN":"MN"
                },
        'FE': {"FE":"FE"
                },
        'CO': {"CO":"CO"
                },
        'NI': {"NI":"NI"
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
    "CR":   1.35,
    "MN":   1.35,
    "FE":   1.34,
    "CO":   1.33,
    "NI":   1.50
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
   "HA3":  [233.88,497.29],
   "CT2":  [301.52,215.37],
   "HA2":  [230.73,494.99],
   "SCD":  [363.90,242.55],
   "CPHG": [256.24,239.08],
   "NR1":  [326.95,163.98],
   "H":    [112.55,479.41],
   "CPH2": [241.53,232.72],
   "HR1":  [202.81,519.32],
   "NR2":  [365.48,259.89],
   "CPHD": [276.60,246.05],
   "HR3":  [213.72,520.57],
   "ZN":   [162.38,190.79],
   "CR":   [159.58,159.65],
   "MN":   [179.81,171.11],
   "FE":   [184.06,179.33],
   "CO":   [193.92,184.94],
   "NI":   [178.44,169.72],
   "CC":   [239.40,197.92],
   "OC":   [370.60,266.93],
   "OW":   [428.23,294.62],
   "HW":   [111.48,474.63],
   "CSD":  [282.80,214.69],
   "HSD":  [236.38,502.27],
   "OH1":  [400.48,297.08]
    }

#CM5 Tkk parameters
cm5tkk_para = {
        ('ZN', 'SCD'): -0.0331,
        ('ZN', 'NR2'): 0.0144,
        ('ZN', 'OW'): 0.0194,
        ('ZN', 'OC'): 0.0074,
        ('ZN', 'OH1'): 0.0092,
        ('CR', 'SCD'): -0.0573,
        ('CR', 'NR2'): 0.0295,
        ('CR', 'OW'): 0.0312,
        ('CR', 'OC'): 0.0019,
        ('CR', 'OH1'): 0.0004,
        ('MN', 'SCD'): -0.0523,
        ('MN', 'NR2'): 0.0300,
        ('MN', 'OW'): 0.0344,
        ('MN', 'OC'): -0.0016,
        ('MN', 'OH1'): -0.0141,
        ('FE', 'SCD'): -0.0456,
        ('FE', 'NR2'): 0.0252,
        ('FE', 'OW'): 0.0304,
        ('FE', 'OC'): -0.0009,
        ('FE', 'OH1'): -0.0140,
        ('CO', 'SCD'): -0.0517,
        ('CO', 'NR2'): 0.0161,
        ('CO', 'OW'): 0.0235,
        ('CO', 'OC'): -0.0056,
        ('CO', 'OH1'): 0.0044,
        ('NI', 'SCD'): -0.0292,
        ('NI', 'NR2'): 0.0112,
        ('NI', 'OW'): 0.0170,
        ('NI', 'OC'): -0.0014,
        ('NI', 'OH1'): 0.0033
        }

