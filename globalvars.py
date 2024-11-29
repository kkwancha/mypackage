import os

MYPACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATABASE_DIR = os.path.join(os.path.dirname(__file__), "DataBase")
DATABASE_DIR = DEFAULT_DATABASE_DIR  # Current database path


ATOMIC_NUM_TO_ELEMENT = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
    35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce",
    59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
    67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta", 74: "W",
    75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb",
    83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
    99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db",
    106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn",
    113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
}

ELEMENT_TO_ATOMIC_NUM = {v: k for k, v in ATOMIC_NUM_TO_ELEMENT.items()}

# Covalent radii from Alvarez (2008)
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 
    'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 
    'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16, 
    'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 
    'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 
    'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 
    'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 
    'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69
}

TRANSITION_METALS = {'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 
                     'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
                     'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'}


# van der Waals radii for elements
# Data from DOI: 10.1039/C3DT50599E, Dalton Trans., 2013, 42, 8617-8636
VDW_RADII = {'H': 1.2, 'He': 1.43, 'Li': 2.12, 'Be': 1.98, 'B': 1.91,
          'C': 1.77, 'N': 1.66, 'O': 1.50, 'F': 1.46, 'Ne': 1.58, 'Na': 2.50,
          'Mg': 2.51, 'Al': 2.25, 'Si': 2.19, 'P': 1.90, 'S': 1.89,
          'Cl': 1.82, 'Ar': 1.83, 'K': 2.73, 'Ca': 2.62, 'Sc': 2.58,
          'Ti': 2.46, 'V': 2.42, 'Cr': 2.45, 'Mn': 2.45, 'Fe': 2.44,
          'Co': 2.40, 'Ni': 2.40, 'Cu': 2.38, 'Zn': 2.39, 'Ga': 2.32,
          'Ge': 2.29, 'As': 1.88, 'Se': 1.82, 'Br': 1.86, 'Kr': 2.25,
          'Rb': 3.21, 'Sr': 2.84, 'Y': 2.75, 'Zr': 2.52, 'Nb': 2.56,
          'Mo': 2.45, 'Tc': 2.44, 'Ru': 2.46, 'Rh': 2.44, 'Pd': 2.15,
          'Ag': 2.53, 'Cd': 2.49, 'In': 2.43, 'Sn': 2.42, 'Sb': 2.47,
          'Te': 1.99, 'I': 2.04, 'Xe': 2.06, 'Cs': 3.48, 'Ba': 3.03,
          'La': 2.98, 'Ce': 2.88, 'Pr': 2.92, 'Nd': 2.95, 'Sm': 2.90,
          'Eu': 2.87, 'Gd': 2.83, 'Tb': 2.79, 'Dy': 2.87, 'Ho': 2.81,
          'Er': 2.83, 'Tm': 2.79, 'Yb': 2.80, 'Lu': 2.74, 'Hf': 2.63,
          'Ta': 2.53, 'W': 2.57, 'Re': 2.49, 'Os': 2.48, 'Ir': 2.41,
          'Pt': 2.29, 'Au': 2.32, 'Hg': 2.45, 'Tl': 2.47, 'Pb': 2.60,
          'Bi': 2.54, 'Ac': 2.8, 'Th': 2.93, 'Pa': 2.88, 'U': 2.71,
          'Np': 2.82, 'Pu': 2.81, 'Am': 2.83, 'Cm': 3.05, 'Bk': 3.4,
          'Cf': 3.05, 'Es': 2.7}

# Bondi van der Waals radii for elements
# From: doi:10.1021/j100881a503
# Accessed: http://www.knowledgedoor.com/2/elements_handbook/bondi_van_der_waals_radius.html 4/8/2021
BONDIVDW = {'Ar': 1.88, 'As': 1.85, 'Br': 1.85, 'Cd': 1.62, 'C': 1.70,
            'Cl': 1.75, 'Cu': 1.4, 'F': 1.47, 'Ga': 1.87, 'Au': 1.66,
            'He': 1.40, 'H': 1.20, 'In': 1.93, 'I': 1.98, 'Kr': 2.02,
            'Pb': 2.02, 'Li': 1.82, 'Mg': 1.73, 'Hg': 1.70, 'Ne': 1.54,
            'Ni': 1.63, 'N': 1.55, 'O': 1.52, 'Pd': 1.63, 'P': 1.80, 'Pt': 1.7,
            'K': 2.75, 'Se': 1.90, 'Si': 2.10, 'Ag': 1.72, 'Na': 2.27, 'S': 1.80,
            'Te': 2.06, 'Tl': 1.96, 'Sn': 1.96, 'U': 1.86, 'Xe': 2.16, 'Zn': 1.39}

# Default database path
DEFAULT_DATABASE_DIR = os.path.join(os.path.dirname(__file__), "DataBase")
DATABASE_DIR = DEFAULT_DATABASE_DIR  # Current database path

def set_database_path(path):
    """
    Set a custom database path.
    """
    global DATABASE_DIR
    DATABASE_DIR = path
    print(f"Database path set to: {DATABASE_DIR}")