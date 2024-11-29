import numpy as np
import pandas as pd
import os
from ..globalvars import ATOMIC_NUM_TO_ELEMENT, COVALENT_RADII, TRANSITION_METALS, ELEMENT_TO_ATOMIC_NUM

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

class xyzinfo:
    """
    Attributes:
    ------------
    logfile (str):
        The path to the log file containing atomic information to be processed.
    xyz_coord (list):
        A list of atomic coordinates, where each item is a list representing 
        [atom_idx, atomic_num, element, x, y, z]:
        - atom_idx (int): The index of the atom (0-based).
        - atomic_num (int): The atomic number of the atom.
        - element (str): The element symbol of the atom.
        - x, y, z (float): The Cartesian coordinates of the atom.
    xyz_df (pandas.DataFrame or None):
        A DataFrame containing atomic information with columns:
        - 'atom_idx': The index of the atom.
        - 'element': The element symbol of the atom.
        - 'x', 'y', 'z': The Cartesian coordinates of the atom.
        - 'cov_radius': The covalent radius of the atom.
        Initially set to None and populated by `_create_xyz_df()`.
    distance_matrix (pandas.DataFrame or None):
        A symmetric matrix representing the distances between all pairs of atoms.
        Initially set to None and populated by `_create_distance_matrix()`.
    lines_forward (list of str):
        A list containing all lines from the log file, read in a forward direction.
    lines_num (int):
        The total number of lines in the log file.
    lines_reverse (list of str):
        A list containing all lines from the log file, read in reverse order.

    Methods Called During Initialization:
    ------------
    get_xyzdata():
        Extracts atomic coordinates from the log file and stores them in `self.xyz_coord`.
    _create_xyz_df():
        Creates a DataFrame to store the atomic information from `self.xyz_coord` in `self.xyz_df`.
    _create_distance_matrix():
        Calculates and stores the distance matrix for the atoms in `self.distance_matrix`.
    """
    def __init__(self, readfrom, xyz_input):
        self.readfrom = readfrom
        self.xyz_input = xyz_input
        self.xyz_content = str()
        self._what_input_type()
        self.xyz_coord = []
        self.get_xyz_coord()
        self.xyz_df = None
        self.get_xyz_df()
        self.distance_matrix = None
        self.get_distance_matrix()
        self.cutoff = None
        
    def _what_input_type(self):
        if self.readfrom == 'file':
            self._read_from_file()
        elif self.readfrom == 'str':
            self._read_from_str()
        else:
            raise ValueError("readfrom must be 'file' or 'str'")
        
    def _read_from_file(self):
        if not os.path.exists(self.xyz_input):
            raise FileNotFoundError(f"No such file: '{self.xyz_input}'")
        with open(self.xyz_input, 'r') as f:
            self.xyz_content = f.read().strip()
        return self.xyz_content
    
    def _read_from_str(self):
        content = self.xyz_input
        content = content.strip()
        self.xyz_content = content
        return self.xyz_content
    
    def get_xyz_coord(self):
        lines = self.xyz_content.splitlines()
        for idx, line in enumerate(lines[2:]):  # Skip the first two lines (header information)
            coord_eachrow = line.split()
            atom_idx = idx
            atom_symbol = coord_eachrow[0]
            atomic_num = ELEMENT_TO_ATOMIC_NUM.get(atom_symbol, np.nan)
            atom_x = float(coord_eachrow[1])
            atom_y = float(coord_eachrow[2])
            atom_z = float(coord_eachrow[3])
            self.xyz_coord.append([atom_idx, atomic_num, atom_symbol, atom_x, atom_y, atom_z])
        return self.xyz_coord
    
    def get_xyz_df(self):
        self.xyz_df = pd.DataFrame(self.xyz_coord, columns=['atom_idx', 'atomic_num', 'element', 'x', 'y', 'z'])
        self.xyz_df['cov_radius'] = self.xyz_df['element'].apply(lambda x: COVALENT_RADII.get(x, 0.0))
        return self.xyz_df
    
    def get_coord_byidx(self, idx):
        coord = self.xyz_df.loc[idx, ['x', 'y', 'z']].astype(float).to_numpy()
        return coord
    
    def get_distance_matrix(self):
        atom_coordinates = self.xyz_df[['x', 'y', 'z']].to_numpy()
        diff = atom_coordinates[:, np.newaxis, :] - atom_coordinates[np.newaxis, :, :]
        distance_matrix = np.sqrt(np.sum(diff**2, axis=-1))
        self.distance_matrix = pd.DataFrame(distance_matrix)
        return self.distance_matrix
    
    def get_bond_existed(self, cutoff=0.1):
        if self.distance_matrix is None:
            raise ValueError("Distance matrix has not been created. Ensure create_distance_matrix() is called.")
        atom_indices = self.xyz_df['atom_idx'].values
        elements = self.xyz_df['element'].values
        cov_radii = self.xyz_df['cov_radius'].values
        bond_existed = []
        for i in range(len(atom_indices)):
            for j in range(i + 1, len(atom_indices)):
                bond_length = self.distance_matrix.iloc[i, j]
                covalent_sum = cov_radii[i] + cov_radii[j] + cutoff

                if bond_length < covalent_sum:
                    bond_existed.append((atom_indices[i], elements[i], atom_indices[j], elements[j], bond_length))

        return bond_existed
    
    def calculate_distance(self, atomA_idx, atomB_idx):
        try:
            atomA_coord = self.get_coord_byidx(atomA_idx)
            atomB_coord = self.get_coord_byidx(atomB_idx)
        except KeyError as e:
            raise ValueError(f"Invalid atom index: {e}")
        dist = np.linalg.norm(atomA_coord - atomB_coord)
        return dist
    
    def calculate_angle(self, atomA_idx, atomB_idx, atomC_idx):
        try:
            atomA_coord = self.get_coord_byidx(atomA_idx)
            atomB_coord = self.get_coord_byidx(atomB_idx)
            atomC_coord = self.get_coord_byidx(atomC_idx)
        except KeyError as e:
            raise ValueError(f"Invalid atom index: {e}")
        vector_BA = atomA_coord - atomB_coord
        vector_BC = atomC_coord - atomB_coord
        norm_BA = np.linalg.norm(vector_BA)
        norm_BC = np.linalg.norm(vector_BC)
        if norm_BA == 0 or norm_BC == 0:
            raise ValueError("One of the vectors has zero length, which makes the angle calculation undefined.")
        uBA = vector_BA / norm_BA
        uBC = vector_BC / norm_BC
        dot_product = np.dot(uBA, uBC)
        dot_product_clipped = np.clip(dot_product, -1.0, 1.0)
        angle_radians = np.arccos(dot_product_clipped)
        angle_degrees = np.degrees(angle_radians)
        return float(angle_degrees)
    
    import numpy as np

    def calculate_dihedral(self, atomA_idx, atomB_idx, atomC_idx, atomD_idx):
        """
        Calculate the dihedral angle between four points A, B, C, and D in 3D space.
        Returns the angle in degrees.
        """
        try:
            atomA_coord = self.get_coord_byidx(atomA_idx)
            atomB_coord = self.get_coord_byidx(atomB_idx)
            atomC_coord = self.get_coord_byidx(atomC_idx)
            atomD_coord = self.get_coord_byidx(atomD_idx)
        except KeyError as e:
            raise ValueError(f"Invalid atom index: {e}")
        vector_BA = atomA_coord - atomB_coord
        vector_BC = atomC_coord - atomB_coord
        vector_CD = atomD_coord - atomC_coord
        n1 = np.cross(vector_BA, vector_BC)
        n2 = np.cross(vector_BC, vector_CD)
        n1 = n1 / np.linalg.norm(n1)
        n2 = n2 / np.linalg.norm(n2)
        dot_product = np.dot(n1, n2)
        dot_product = np.clip(dot_product, -1.0, 1.0)
        angle_radians = np.arccos(dot_product)
        angle_degrees = np.degrees(angle_radians)
        return angle_degrees
    
    def atomidx_from_element(self, element):
        list_indices = self.xyz_df[self.xyz_df['element'] == element].index.tolist()
        # list_indices = [3, 4, 14, 15, 16, 19, 21, 24, 28, 32] of the element
        return list_indices
    
    def atomnearby_from_idx(self, atomA_idx, cutoff=0.1):
        bond_existed = self.get_bond_existed(cutoff=cutoff)
        nearby_atoms = set()
        
        for bond in bond_existed:
            if atomA_idx in bond:
                atom1 = (int(bond[0]),bond[1])
                atom2 = (int(bond[2]),bond[3])
                nearby_atoms.add(atom1)
                nearby_atoms.add(atom2)
        nearby_atoms = {atom for atom in nearby_atoms if atom[0] != atomA_idx}
        # nearby_atoms = {(0, 'P'), (1, 'Ni'), (2, 'O'), (3, 'C'), (4, 'C'),}
        return nearby_atoms
    
    def get_unique_elements(self):
        unique_elements = set(self.xyz_df['element'].unique())
        return unique_elements

    def find_transition_metal(self):
        unique_elements = self.get_unique_elements()
        metals = unique_elements.intersection(TRANSITION_METALS)
        if len(metals) == 0:
            return None
        return metals # set of metal {'Ni', 'Pd'}
    
    def isbond(self, atomA_idx, atomB_idx, cutoff=0.1):
        cov_radii_A = self.xyz_df.loc[atomA_idx, 'cov_radius']
        cov_radii_B = self.xyz_df.loc[atomB_idx, 'cov_radius']
        distance = self.calculate_distance(atomA_idx, atomB_idx)
        if distance < (cov_radii_A + cov_radii_B + cutoff):
            return True
        else:
            return False
    
    def sum_cov_radii(self, elementA, elementB, cutoff=0.0):
        cov_radii_A = COVALENT_RADII[elementA]
        cov_radii_B = COVALENT_RADII[elementB]
        distance = cov_radii_A + cov_radii_B + cutoff
        return distance
    
    def find_ligand_mononuclear(self, getonlyidx=False):
        metal_element = list(self.find_transition_metal())[0] # only one metal center
        metal_idx = self.atomidx_from_element(metal_element)[0]
        bondtometal = self.atomnearby_from_idx(metal_idx) # cutoff = 0.1
        if getonlyidx:
            bondtometal = {pair[0] for pair in bondtometal}
        return bondtometal
        
    def find_nearest_ofelement_toatom(self, center_atom_index, target_element, cutoff=0.1):
        """
        Finds the nearest atom of a specified element to a given atom.

        Parameters:
        -----------
        center_atom_index : int
            The index of the central atom.
        target_element : str
            The element symbol of the target atom to find.

        Returns:
        --------
        nearest_atom : tuple
            A tuple (atom_index, element) representing the nearest atom of the specified element.
        """
        bond_to_center_atom = self.atomnearby_from_idx(center_atom_index, cutoff)
        atoms_with_target_element = [atom for atom in bond_to_center_atom if atom[1] == target_element]
        nearest_atom = min(atoms_with_target_element, key=lambda atom: self.calculate_distance(center_atom_index, atom[0]))
        return nearest_atom
    
    def is_planar(self, list_of_atoms, tolerance=5):
        """
        Check if all atoms in the list are coplanar.
        The method calculates the dihedral angles formed by the first three atoms in the list
        and each of the remaining atoms. If all angles are within the tolerance (near 0° or 180°),
        the atoms are coplanar.
        
        Parameters:
        - list_of_atoms: List of atom indices
        - tolerance: The allowed deviation from 0 or 180 degrees for planarity (in degrees)
        
        Returns True if all atoms are coplanar within the given tolerance, False otherwise.
        """
        if len(list_of_atoms) < 4:
            raise ValueError("At least 4 atom indices are required to check for planarity.")
        atomA_idx, atomB_idx, atomC_idx = list_of_atoms[0], list_of_atoms[1], list_of_atoms[2]
        for i in range(3, len(list_of_atoms)):
            atomD_idx = list_of_atoms[i]
            dihedral = self.calculate_dihedral(atomA_idx, atomB_idx, atomC_idx, atomD_idx)
            if abs(dihedral) > tolerance and abs(dihedral - 180) > tolerance:
                return False  # If any atom is not in planar, return False
        return True
    
    def is_sqp(self, tolerance=5):
        """
        Check if the four atoms A, B, C, and D are planar.
        Returns True if the atoms are coplanar within a given tolerance (in degree), False otherwise.
        """
        if len(self.find_ligand_mononuclear()) != 4:
            raise ValueError("Square planar geometry requires exactly 4 ligands.")
        
        ligatoms = list(self.find_ligand_mononuclear())
        
        print(ligatoms)
        
        atomA_idx = ligatoms[0][0]  # First ligand's coordinates
        atomB_idx = ligatoms[1][0]  # Second ligand's coordinates
        atomC_idx = ligatoms[2][0]  # Third ligand's coordinates
        atomD_idx = ligatoms[3][0]  # Fourth ligand's coordinates
    
        angle = self.calculate_dihedral(atomA_idx, atomB_idx, atomC_idx, atomD_idx)
        print(angle)
        if abs(angle) < tolerance or abs(angle - 180) < tolerance:
            return True
        return False
    
    def find_normv(self, list_of_atoms):
        """
        Find the normalized normal vector of the plane defined by three atoms.
        Returns a consistent normal vector regardless of the order of atoms.
        """
        if len(list_of_atoms) != 3:
            raise ValueError('You need exactly 3 atoms to find the normal vector of a plane.')
        atomA_idx, atomB_idx, atomC_idx = list_of_atoms
        try:
            atomA_coord = self.get_coord_byidx(atomA_idx)
            atomB_coord = self.get_coord_byidx(atomB_idx)
            atomC_coord = self.get_coord_byidx(atomC_idx)
        except KeyError as e:
            raise ValueError(f"Invalid atom index: {e}")
        vector_BA = atomA_coord - atomB_coord
        vector_BC = atomC_coord - atomB_coord
        normv = np.cross(vector_BA, vector_BC)
        normv = normv / np.linalg.norm(normv)
        if normv[0] < 0:  # Flip if the x-component is negative
            normv = -normv
        
        return normv
    
    def angle_twovec(self, vector_BA, vector_BC):
        norm_BA = np.linalg.norm(vector_BA)
        norm_BC = np.linalg.norm(vector_BC)
        if norm_BA == 0 or norm_BC == 0:
            raise ValueError("One of the vectors has zero length, which makes the angle calculation undefined.")
        uBA = vector_BA / norm_BA
        uBC = vector_BC / norm_BC
        dot_product = np.dot(uBA, uBC)
        dot_product_clipped = np.clip(dot_product, -1.0, 1.0)
        angle_radians = np.arccos(dot_product_clipped)
        angle_degrees = np.degrees(angle_radians)
        return float(angle_degrees)