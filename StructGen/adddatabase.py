import numpy as np
import networkx as nx
from openbabel import openbabel
import os
# import yaml
from ..globalvars import MYPACKAGE_DIR, ATOMIC_NUM_TO_ELEMENT, COVALENT_RADII, TRANSITION_METALS, ELEMENT_TO_ATOMIC_NUM
from ..config import DEFAULT_DIR
from ..AnalyzeGeom import geomobj

DEFAULT_DATABASE_DIR = os.path.join(MYPACKAGE_DIR, 'DataBase')
DATABASE_DIR = DEFAULT_DATABASE_DIR

# input geom: SMILE, mol2, xyz
# for smile, generate 3d structure using openbabel

# input catoms: list of connecting atoms, 0-based index

# # Directory for storing ligands and metadata
# LIGAND_DIR = os.path.join("ligands")
# METADATA_FILE = os.path.join(LIGAND_DIR, "Ligands.yaml")


# def set_database_path(path):
#     """
#     Set a custom database path.
#     """
    
#     fullpath = os.path.abspath(path)
#     global DATABASE_DIR
#     DATABASE_DIR = fullpath
#     print(f"Database path set to: {DATABASE_DIR}")
    
# def get_current_database_path():
#     print(DATABASE_DIR)
    
# def copy_to_custom_db():
#     print(DEFAULT_DATABASE_DIR)
    
# # get ligand

# # def add_ligand(xyzfile, )


# def initialize_ligand_directory():
#     """Ensure the ligand directory and metadata file exist."""
#     if not os.path.exists(LIGAND_DIR):
#         os.makedirs(LIGAND_DIR)
#     if not os.path.exists(METADATA_FILE):
#         with open(METADATA_FILE, 'w') as f:
#             yaml.dump({}, f)

# # def save_structure(mol, file_name, file_format):
# #     """Save the molecular structure to a file."""
# #     if not mol.write(file_format, file_name):
# #         raise ValueError(f"Failed to write structure to {file_name}.")

# def add_ligand(ligand_name, catoms, structure=None, file_format=None, description=None):
#     """
#     Add a new ligand to the ligand database.

#     Parameters:
#     - ligand_name: str, name of the ligand
#     - catoms: list of int, indices of connecting atoms (0-based)
#     - structure: str, SMILES string, MOL2 content, or XYZ content
#     - file_format: str, format of the structure (smiles, mol2, xyz)
#     - description: str, optional description of the ligand
#     """
#     initialize_ligand_directory()

#     # Load existing metadata
#     with open(METADATA_FILE, 'r') as f:
#         metadata = yaml.safe_load(f) or {}

#     if ligand_name in metadata:
#         raise ValueError(f"Ligand {ligand_name} already exists in metadata.")

#     # Generate 3D structure if input is SMILES
#     if file_format == "smiles":
#         ob_conversion = openbabel.OBConversion()
#         ob_conversion.SetInAndOutFormats("smi", "mol2")
#         mol = openbabel.OBMol()
#         if not ob_conversion.ReadString(mol, structure):
#             raise ValueError("Failed to parse SMILES string.")
#         pybel_mol = pybel.Molecule(mol)
#         pybel_mol.make3D()
#         structure = pybel_mol.write("xyz")
#         file_format = "xyz"

#     # Save structure file
#     file_path = os.path.join(LIGAND_DIR, f"{ligand_name}.{file_format}")
#     with open(file_path, 'w') as f:
#         f.write(structure)

#     # Update metadata
#     metadata[ligand_name] = {
#         "format": file_format,
#         "catoms": catoms,
#         "description": description or f"Ligand {ligand_name}."
#     }

#     # Save updated metadata
#     with open(METADATA_FILE, 'w') as f:
#         yaml.dump(metadata, f, default_flow_style=False)

#     print(f"Ligand {ligand_name} added successfully.")

# # # Example usage
# # add_ligand(
# #     ligand_name="L002PH",
# #     catoms=[0, 1],
# #     structure="C1=CC=CC=C1",  # Example SMILES for benzene
# #     file_format="smiles",
# #     description="Aromatic ligand example."
# # )