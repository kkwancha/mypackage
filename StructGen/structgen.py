import glob, os
import numpy as np
import networkx as nx
from ..globalvars import MYPACKAGE_DIR
from ..AnalyzeGeom import geommath, geomobj

TEMPLATE_GEOM = {'li'  : 'li.xyz',
                 'lig' : 'lig.xyz',
                 'sqp' : 'sqp.xyz',
                 'tsoa': 'tsoa.xyz',
                 'tst' : 'tst.xyz',
                 'tsre': 'tsre.xyz',
                 'i4':   'i4.xyz',
                 'I4_new':'I4_new.xyz',
                 'tsoa1coord': 'tsoa1coord.xyz'}

def create_template_Mol(geom):
    geom_template_path = TEMPLATE_GEOM[geom]
    file_path = os.path.join(MYPACKAGE_DIR, 'DataBase/Geometries/', geom_template_path)
    template_Mol = geomobj.Mol()
    template_Mol.readfile('xyz', file_path)
    print('hi')
    return template_Mol

def add_metal(metal_center):
    """
    Modify metal centers in template_Mol based on the specified metal_center list.

    Parameters
    ----------
    metal_center : list of str
        List of metal element symbols to replace atoms with symbol 'M'.

    Returns
    -------
    template_Mol : Mol
        The modified molecule with updated metal centers.
    """
    template_Mol = create_template_Mol()
    atoms_template_M = [atom for atom in template_Mol.atoms if atom.sym == 'M']
    if len(metal_center) < len(atoms_template_M):
        raise ValueError("Not enough elements in metal_center to match the atoms with symbol 'M'.")
    for idx, atom in enumerate(atoms_template_M):
        atom.sym = metal_center[idx % len(metal_center)]
    return template_Mol  # Return the modified template_Mol

def getbondlength_fromsym(sym1, sym2):
    BONDLENGTH = {
        frozenset(['Pd', 'P']): 2.310,
        frozenset(['Ni', 'P']): 2.213,
        frozenset(['Pd', 'C']): 1.994,
        frozenset(['Ni', 'C']): 1.890
    }
    
    bond_pair = frozenset([sym1, sym2])  # Use frozenset to ignore order of elements
    dist = BONDLENGTH.get(bond_pair)
    
    if dist is None:
        raise ValueError(f"No bond distance found for pair ({sym1}, {sym2})")
    
    return dist

# def genstruct(geom, metal_center, ligand_path, catoms):
#     # 1. Get mol template
#     template_Mol = create_template_Mol(geom)
#     template_atoms = template_Mol.atoms
    
#     # 2. Get atoms to replace
#     template_atoms_M = [atom for atom in template_atoms if atom.sym == 'M']
#     template_atoms_X = [atom for atom in template_atoms if atom.sym == 'X']
#     print(template_atoms_M)
    
#     # 3. Define metal center
#     metal_center = ['Pd'] # input: list of metal (either mononuclear, dinuclear)
#     if len(metal_center) < len(template_atoms_M):
#         raise ValueError("Not enough elements in metal_center to match the atoms with symbol 'M'.")
#     for idx, atom in enumerate(template_atoms_M):
#         atom.sym = metal_center[idx % len(metal_center)]
    
#     # 4. Get ligand Mol
#     lig_Mol = cpx.MolLigand()
#     lig_Mol.readfile('mol2', ligand_path)
#     lig_Mol.set_catoms(catoms)
    
#     # 5. Get M-L bond length
#     metal_sym = template_atoms[0].sym
#     ligand_sym = list(lig_Mol.getAtom_catoms())[0].sym # monodentate ligand
#     ML_bondlength = getbondlength_fromsym(metal_sym, ligand_sym)
    
#     # 6. Scale the M-L bond length
#     pointM = template_atoms[0].coord
#     pointL = template_atoms[1].coord # iterate for other points
#     pointL = geommath.scaled_distance(pointM, pointL, ML_bondlength)
    
#     # 7. Translate ligand to the connecting point
#     catom = list(lig_Mol.getAtom_catoms())[0] # monodentate ligand
#     translated_vector = pointL - catom.coord
#     translated_coord = lig_Mol.coords + translated_vector
    
#     # 8. Find lone pair vector
#     adjcatoms = lig_Mol.getAtom_adjcatom()[catoms[0]] # monodentate ligand
#     adjcatom_coords = lig_Mol.getcoord_fromatomlist(adjcatoms)
#     adjcatom_vecs = adjcatom_coords - catom.coord
#     lonepair_vec = geommath.normalize(-np.sum(adjcatom_vecs, axis=0))
    
#     # 9. Find bonding vector
#     template_bondvec = geommath.normalize(pointM - pointL)
    
#     # 10. Rotate ligands to bond the template
#     rotation_axis = np.cross(lonepair_vec, template_bondvec)
#     rotation_angle = geommath.angle_between_vectors(lonepair_vec, template_bondvec)
#     R = geommath.rotation_matrix(rotation_axis, rotation_angle)
#     rotated_coords = geommath.rotate_coordinates(translated_coord, pointL, R)
    
#     # 11. Add ligands to the template complex
#     for sym, coord in zip(lig_Mol.getlistofsym(), rotated_coords):
#         num_atoms = template_Mol.natoms
#         template_Mol.addAtom(cpx.Atom(num_atoms,sym,coord))
    
#     # 12. Delete template atoms
#     template_Mol.deleteAtoms_bysym(['X'])
    
#     return template_Mol


def genstruct_M(geom, list_metal_center, list_ligand_path, list_ligand_catoms):
    """
    Generate a molecular structure by defining the metal center and adding multiple ligands to a template.

    Parameters
    ----------
    geom : str
        The name of the geometry template (e.g., 'tsoa').
    metal_center : list of str
        List of metal symbols (e.g., ['Pd'] for mononuclear, ['Pd', 'Pd'] for dinuclear).
    ligand_paths : list of str
        Paths to the ligand MOL2 files.
    catoms : list of list of int
        List of atom indices for each ligand to set as coordination atoms. Each entry is a list of indices.

    Returns
    -------
    cpx.Mol
        The generated complex molecule with the specified metal and ligands.
    """
    # 1. Get mol template
    template_Mol = create_template_Mol(geom)
    template_atoms = template_Mol.atoms

    # 2. Get atoms to replace
    template_atoms_M = [atom for atom in template_atoms if atom.sym == 'M']
    template_atoms_X = [atom for atom in template_atoms if atom.sym == 'X']
    
    # 3. Define metal center
    if len(list_metal_center) < len(template_atoms_M):
        raise ValueError("Not enough elements in list_metal_center to match the atoms with symbol 'M'.")
    for idx, atom in enumerate(template_atoms_M):
        atom.sym = list_metal_center[idx % len(list_metal_center)]
    
    # 4. add ligand
    if len(list_ligand_path) != len(template_atoms_X):
        raise ValueError("Not enough elements in list_ligand_path to match the atoms with symbol 'X'.")
    
    for idx, template_atom_X in enumerate(template_atoms_X):
        ligand_path = list_ligand_path[idx]
        ligand_catoms = list_ligand_catoms[idx]
        # Step 1: Load the ligand molecule
        lig_Mol = geomobj.MolLigand()
        lig_Mol.readfile('mol2', ligand_path)
        lig_Mol.set_catoms(list_ligand_catoms[idx])
        
        # Step 2: Get the bond length for the metal-ligand interaction
        metal_sym = template_atoms_M[0].sym # mononuclear
        ligand_sym = list(lig_Mol.getAtom_catoms())[0].sym # monodentate ligand
        try:
            ML_bondlength = getbondlength_fromsym(metal_sym, ligand_sym)
        except ValueError as e:
            raise ValueError(f"Failed to find bond length for ({metal_sym}, {ligand_sym}): {e}")
        
        # Step 3: Scale the M-L bond length
        pointM = template_atoms_M[0].coord
        pointX = template_atom_X.coord
        pointX = geommath.scaled_distance(pointM, pointX, ML_bondlength)
        
        # Step 4: Translate ligand to the connecting point
        catom = list(lig_Mol.getAtom_catoms())[0]  # Monodentate ligand
        translated_vector = pointX - catom.coord
        translated_coord = lig_Mol.coords + translated_vector
        
        # Step 5: Find lone pair vector
        adjcatoms = lig_Mol.getAtom_adjcatom()[ligand_catoms[0]]  # Monodentate ligand
        adjcatom_coords = lig_Mol.getcoord_fromatomlist(adjcatoms)
        adjcatom_vecs = adjcatom_coords - catom.coord
        lonepair_vec = geommath.normalize(-np.sum(adjcatom_vecs, axis=0))

        # Step 6: Find bonding vector
        template_bondvec = geommath.normalize(pointM - pointX)

        # Step 7: Rotate ligands to bond the template
        rotation_axis = np.cross(lonepair_vec, template_bondvec)
        rotation_angle = geommath.angle_between_vectors(lonepair_vec, template_bondvec)
        R = geommath.rotation_matrix(rotation_axis, rotation_angle)
        rotated_coords = geommath.rotate_coordinates(translated_coord, pointX, R)

        # Step 8: Add rotated ligand atoms to the template molecule
        for sym, coord in zip(lig_Mol.getlistofsym(), rotated_coords):
            num_atoms = template_Mol.natoms
            template_Mol.addAtom(geomobj.Atom(num_atoms, sym, coord))
    
    # 12. Delete template atoms
    template_Mol.deleteAtoms_bysym(['X'])
    
    return template_Mol

# lig = [name for name, count in zip(lig_name, lig_occ) for _ in range(count)]
def genstruct(geom, list_metal_center, list_ligand_path, list_ligand_catoms):
    """
    Generate a molecular structure by defining the metal center and adding multiple ligands to a template.

    Parameters
    ----------
    geom : str
        The name of the geometry template (e.g., 'tsoa').
    metal_center : list of str
        List of metal symbols (e.g., ['Pd'] for mononuclear, ['Pd', 'Pd'] for dinuclear).
    ligand_paths : list of str
        Paths to the ligand MOL2 files.
    catoms : list of list of int
        List of atom indices for each ligand to set as coordination atoms. Each entry is a list of indices.

    Returns
    -------
    cpx.Mol
        The generated complex molecule with the specified metal and ligands.
    """
    # 1. Get mol template
    template_Mol = create_template_Mol(geom)
    template_atoms = template_Mol.atoms

    # 2. Get atoms to replace
    template_atoms_M = [atom for atom in template_atoms if atom.sym == 'M']
    template_atoms_X = [atom for atom in template_atoms if atom.sym == 'X']
    
    # 3. Define metal center
    if len(list_metal_center) < len(template_atoms_M):
        raise ValueError("Not enough elements in list_metal_center to match the atoms with symbol 'M'.")
    for idx, atom in enumerate(template_atoms_M):
        atom.sym = list_metal_center[idx % len(list_metal_center)]
    
    # 4. add ligand
    if len(list_ligand_path) != len(template_atoms_X):
        raise ValueError("Not enough elements in list_ligand_path to match the atoms with symbol 'X'.")
    
    for idx, template_atom_X in enumerate(template_atoms_X):
        ligand_path = list_ligand_path[idx]
        ligand_catoms = list_ligand_catoms[idx]
        # Step 1: Load the ligand molecule
        lig_Mol = geomobj.MolLigand()
        lig_Mol.readfile('mol2', ligand_path)
        lig_Mol.set_catoms(list_ligand_catoms[idx])
        
        # Step 2: Get the bond length for the metal-ligand interaction
        metal_sym = template_atoms_M[0].sym # mononuclear
        ligand_sym = list(lig_Mol.getAtom_catoms())[0].sym # monodentate ligand
        try:
            ML_bondlength = getbondlength_fromsym(metal_sym, ligand_sym)
        except ValueError as e:
            raise ValueError(f"Failed to find bond length for ({metal_sym}, {ligand_sym}): {e}")
        
        # Step 3: Scale the M-L bond length
        pointM = template_atoms_M[0].coord
        pointX = template_atom_X.coord
        pointX = geommath.scaled_distance(pointM, pointX, ML_bondlength)
        
        # Step 4: Translate ligand to the connecting point
        catom = list(lig_Mol.getAtom_catoms())[0]  # Monodentate ligand
        translated_vector = pointX - catom.coord
        translated_coord = lig_Mol.coords + translated_vector
        
        # Step 5: Find lone pair vector
        adjcatoms = lig_Mol.getAtom_adjcatom()[ligand_catoms[0]]  # Monodentate ligand
        adjcatom_coords = lig_Mol.getcoord_fromatomlist(adjcatoms)
        adjcatom_vecs = adjcatom_coords - catom.coord
        lonepair_vec = geommath.normalize(-np.sum(adjcatom_vecs, axis=0))

        # Step 6: Find bonding vector
        template_bondvec = geommath.normalize(pointM - pointX)

        # Step 7: Rotate ligands to bond the template
        rotation_axis = np.cross(lonepair_vec, template_bondvec)
        rotation_angle = geommath.angle_between_vectors(lonepair_vec, template_bondvec)
        R = geommath.rotation_matrix(rotation_axis, rotation_angle)
        rotated_coords = geommath.rotate_coordinates(translated_coord, pointX, R)
        
        # Step 8: Update the ligand coordinates with rotated coordinates
        for i, coord in enumerate(rotated_coords):
            lig_Mol.atoms[i].coord = coord
            
        # Step 9: Add the ligand to the template molecule
        template_Mol.addMol(lig_Mol)

    # 12. Delete template atoms
    template_Mol.deleteAtoms_bysym(['X'])
    
    return template_Mol, lig_Mol

def genstruct2(geom, list_metal_center, list_ligand_path, list_ligand_catoms):
    """
    Generate a molecular structure by defining the metal center and adding multiple ligands to a template.

    Parameters
    ----------
    geom : str
        The name of the geometry template (e.g., 'tsoa').
    metal_center : list of str
        List of metal symbols (e.g., ['Pd'] for mononuclear, ['Pd', 'Pd'] for dinuclear).
    ligand_paths : list of str
        Paths to the ligand MOL2 files.
    catoms : list of list of int
        List of atom indices for each ligand to set as coordination atoms. Each entry is a list of indices.

    Returns
    -------
    cpx.Mol
        The generated complex molecule with the specified metal and ligands.
    """
    # 1. Get mol template
    template_Mol = create_template_Mol(geom)
    template_atoms = template_Mol.atoms

    # 2. Get atoms to replace
    template_atoms_M = [atom for atom in template_atoms if atom.sym == 'M']
    template_atoms_X = [atom for atom in template_atoms if atom.sym == 'X']
    
    # 3. Define metal center
    if len(list_metal_center) < len(template_atoms_M):
        raise ValueError("Not enough elements in list_metal_center to match the atoms with symbol 'M'.")
    for idx, atom in enumerate(template_atoms_M):
        atom.sym = list_metal_center[idx % len(list_metal_center)]
    
    # 4. add ligand
    if len(list_ligand_path) != len(template_atoms_X):
        print(list_ligand_path)
        print(template_atoms_X)
        raise ValueError("Not enough elements in list_ligand_path to match the atoms with symbol 'X'.")
    
    for idx, template_atom_X in enumerate(template_atoms_X):
        ligand_path = list_ligand_path[idx]
        ligand_catoms = list_ligand_catoms[idx]
        # Step 1: Load the ligand molecule
        lig_Mol = geomobj.MolLigand()
        lig_Mol.readfile('mol2', ligand_path)
        lig_Mol.set_catoms(list_ligand_catoms[idx])
        
        # Step 2: Get the bond length for the metal-ligand interaction
        metal_sym = template_atoms_M[0].sym # mononuclear
        ligand_sym = list(lig_Mol.getAtom_catoms())[0].sym # monodentate ligand
        try:
            ML_bondlength = getbondlength_fromsym(metal_sym, ligand_sym)
        except ValueError as e:
            raise ValueError(f"Failed to find bond length for ({metal_sym}, {ligand_sym}): {e}")
        
        # Step 3: Scale the M-L bond length
        pointM = template_atoms_M[0].coord
        pointX = template_atom_X.coord
        pointX = geommath.scaled_distance(pointM, pointX, ML_bondlength)
        
        # Step 4: Translate ligand to the connecting point
        catom = list(lig_Mol.getAtom_catoms())[0]  # Monodentate ligand
        translated_vector = pointX - catom.coord
        translated_coord = lig_Mol.coords + translated_vector
        
        # Step 5: Find lone pair vector
        adjcatoms = lig_Mol.getAtom_adjcatom()[ligand_catoms[0]]  # Monodentate ligand
        adjcatom_coords = lig_Mol.getcoord_fromatomlist(adjcatoms)
        adjcatom_vecs = adjcatom_coords - catom.coord
        lonepair_vec = geommath.normalize(-np.sum(adjcatom_vecs, axis=0))

        # Step 6: Find bonding vector
        template_bondvec = geommath.normalize(pointM - pointX)

        # Step 7: Rotate ligands to bond the template
        rotation_axis = np.cross(lonepair_vec, template_bondvec)
        rotation_angle = geommath.angle_between_vectors(lonepair_vec, template_bondvec)
        R = geommath.rotation_matrix(rotation_axis, rotation_angle)
        rotated_coords = geommath.rotate_coordinates(translated_coord, pointX, R)
        
        # Step 8: Update the ligand coordinates with rotated coordinates
        for i, coord in enumerate(rotated_coords):
            lig_Mol.atoms[i].coord = coord
            
        # Step 9: Add the ligand to the template molecule
        template_Mol.addMol(lig_Mol)

    # 12. Delete template atoms
    template_Mol.deleteAtoms_bysym(['X'])
    
    return template_Mol