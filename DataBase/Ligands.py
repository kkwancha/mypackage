import numpy as np
import pandas as pd
import networkx as nx
from openbabel import openbabel
import os
from ..globalvars import MYPACKAGE_DIR, ATOMIC_NUM_TO_ELEMENT, COVALENT_RADII, TRANSITION_METALS, ELEMENT_TO_ATOMIC_NUM
from ..AnalyzeGeom import geomobj

# input geom: SMILE, mol2, xyz, mol, sdf
# input catoms

# def add_ligand(MolLigand, catoms):


