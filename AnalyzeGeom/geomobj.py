import numpy as np
import pandas as pd
import networkx as nx
from openbabel import openbabel
import os
from ..globalvars import ATOMIC_NUM_TO_ELEMENT, COVALENT_RADII, TRANSITION_METALS, ELEMENT_TO_ATOMIC_NUM

def remove_digits(element_with_numbering):
    # Iterate through the string and keep only non-digit characters
    element = ''.join([char for char in element_with_numbering if not char.isdigit()])
    return element

def normalize(vector):
    """ Normalize a vector. """
    return vector / np.linalg.norm(vector)

class Atom:
    def __init__(self, idx, sym, coord):
        self.idx = idx
        self.sym = sym
        self.coord_x = coord[0]
        self.coord_y = coord[1]
        self.coord_z = coord[2]
        self.coord = np.array(coord)
        
    def __repr__(self):
        return f"Atom({self.sym}{self.idx})"
    
    @property
    def atomic_num(self):
        return ELEMENT_TO_ATOMIC_NUM.get(self.sym, np.nan)
    
    @property
    def cov_radii(self):
        return COVALENT_RADII.get(self.sym, 0.0)
    
    def distance_to(self, other):
        return np.linalg.norm(self.coord - other.coord)

class Mol:
    def __init__(self):
        self.atoms = []
        self.bonds = dict()
        self.xyz_df = None
        self.distance_matrix = None
        self.bond_matrix = None
        self.G_Graph = nx.Graph()
        self.G_MultiGraph = nx.MultiGraph()

    @property
    def natoms(self):
        """Dynamically calculate the number of atoms based on the length of self.atoms."""
        return len(self.atoms)
    
    @property
    def nbonds(self):
        """Dynamically calculate the number of bonds based on the length of self.atoms."""
        return len(self.bonds)

    @property
    def coords(self):
        """Dynamically retrieve coordinates from atoms."""
        return np.array([atom.coord for atom in self.atoms])

    def __repr__(self):
        return f"Mol(natoms={self.natoms})"
    
    def readfile(self, format, inputgeom):
        """
        Read geometry data from a file.
        
        Parameters:
        - inputgeom: Path to the geometry file.
        - file_format: Format of the file (e.g., 'xyz', 'mol2').
        """
        if not os.path.exists(inputgeom):
            raise FileNotFoundError(f"No such file: '{inputgeom}'")
        with open(inputgeom, 'r') as file:
            content = file.read().strip()
        self._parse_read(content, format)
        return self
    
    def readstr(self, inputgeom, file_format):
        """
        Read geometry data from a string.
        
        Parameters:
        - inputgeom: The geometry data as a string.
        - file_format: Format of the data (e.g., 'xyz', 'mol2').
        """
        content = inputgeom.strip()
        self._parse_read(content, file_format)
        return self
    
    def _parse_read(self, content, format):
        """
        Parse the geometry content based on the specified file format.
        
        Parameters:
        - content: The geometry content (string).
        - file_format: Format of the data (e.g., 'xyz', 'mol2').
        """
        if format == 'xyz':
            self.readxyz(content)
        elif format == 'mol2':
            self.readmol2(content)
        else:
            raise ValueError(f"Unsupported format: '{format}'")
    
    def readxyz(self, content):
        """
        Parse the XYZ content to create Atom objects.

        Parameters
        ----------
        content : str
            The XYZ format content as a string.
        """
        lines = content.splitlines()
        try:
            num_atoms = int(lines[0].strip())
            self.atoms = []  # Initialize or clear the atoms list

            for idx, line in enumerate(lines[2:num_atoms + 2], start=0):  # Start at index 0
                parts = line.split()
                if len(parts) != 4:
                    raise ValueError(f"Invalid line in XYZ data: {line}")
                element = parts[0]
                x, y, z = map(lambda x: x.replace('–','-'), parts[1:4])
                x, y, z = map(float, [x, y, z])
                self.atoms.append(Atom(idx, element, [x, y, z]))  # Create Atom instances directly in self.atoms
        except (ValueError, IndexError) as e:
            print(e)
            raise ValueError("Invalid XYZ format") from e
        # self.bond_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        return self
        
    
    def readmol2_1(self, content):
        """
        Parse the MOL2 content to create Atom objects and bond matrix.

        Parameters
        ----------
        content : str
            The MOL2 format content as a string.
        """
        lines = content.splitlines()
        section_molecule = None
        section_atom = None
        section_bond = None
        num_atoms = 0
        num_bonds = 0

        for idx, line in enumerate(lines):
            if line.startswith('@<TRIPOS>MOLECULE'):
                section_molecule = idx
            elif line.startswith('@<TRIPOS>ATOM'):
                section_atom = idx
            elif line.startswith('@<TRIPOS>BOND'):
                section_bond = idx

        if section_atom is None or section_bond is None:
            raise ValueError("ATOM or BOND section not found in MOL2 content.")

        # Parse the MOLECULE section
        try:
            num_atoms = int(lines[section_molecule + 2].split()[0])
            num_bonds = int(lines[section_molecule + 2].split()[1])
        except (IndexError, ValueError) as e:
            raise ValueError(f"Could not read number of atoms from the file: {e}")
        
        # Parse the ATOM section
        self.atoms = np.empty(num_atoms, dtype=object)
        if section_atom is not None:
            coords = lines[section_atom+1 : section_atom+1+num_atoms]
            for line in coords:
                coord_list = line.split()
                atom_idx = int(coord_list[0]) - 1  # Convert to 0-based index
                element = remove_digits(coord_list[1])
                x, y, z = map(float, coord_list[2:5])
                self.atoms[atom_idx] = Atom(atom_idx, element, [x, y, z])
        else:
            raise ValueError("ATOM section not found in the file")

        # Parse the BOND section
        self.bond_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        for line in lines[section_bond+1:section_bond+1+num_bonds]:
            parts = line.split()
            if len(parts) < 3:
                continue
            atom1_idx = int(parts[1]) - 1  # Convert to 0-based index
            atom2_idx = int(parts[2]) - 1  # Convert to 0-based index
            bond_type = parts[3]
            # Determine bond order based on bond type
            if bond_type == "1":
                bond_order = 1.0  # Single bond
            elif bond_type == "2":
                bond_order = 2.0  # Double bond
            elif bond_type == "3":
                bond_order = 3.0  # Triple bond
            elif bond_type == "ar":
                bond_order = 1.5  # Aromatic bond
            else:
                bond_order = 0.0  # Undefined or unsupported bond type; ignore
            self.bond_matrix[atom1_idx, atom2_idx] = bond_order
            self.bond_matrix[atom2_idx, atom1_idx] = bond_order  # Symmetric for undirected bonds
        return self
    
    def readmol2(self, content):
        """
        Parse the MOL2 content to create Atom objects and bond matrix.

        Parameters
        ----------
        content : str
            The MOL2 format content as a string.
        """
        lines = content.splitlines()
        section_MOLECULE = None
        section_ATOM = None
        section_BOND = None
        num_atoms = 0
        num_bonds = 0

        for idx, line in enumerate(lines):
            if line.startswith('@<TRIPOS>MOLECULE'):
                section_molecule = idx
            elif line.startswith('@<TRIPOS>ATOM'):
                section_ATOM = idx
            elif line.startswith('@<TRIPOS>BOND'):
                section_BOND = idx

        if section_ATOM is None or section_BOND is None:
            raise ValueError("ATOM or BOND section not found in MOL2 content.")

        # Parse the MOLECULE section
        try:
            num_atoms = int(lines[section_molecule + 2].split()[0])
            num_bonds = int(lines[section_molecule + 2].split()[1])
        except (IndexError, ValueError) as e:
            raise ValueError(f"Could not read number of atoms from the file: {e}")
        
        # Parse the ATOM section
        self.atoms = np.empty(num_atoms, dtype=object)
        if section_ATOM is not None:
            for line in lines[section_ATOM+1 : section_ATOM+1+num_atoms]:
                coord_list = line.split()
                atom_idx = int(coord_list[0]) - 1  # Convert to 0-based index
                element = remove_digits(coord_list[1])
                x, y, z = map(float, coord_list[2:5])
                self.atoms[atom_idx] = Atom(atom_idx, element, [x, y, z])
        else:
            raise ValueError("ATOM section not found in the file")

        # Parse the BOND section
        self.bond = np.empty(num_bonds, dtype=object)
        if section_BOND is not None:
            for line in lines[section_BOND+1:section_BOND+1+num_bonds]:
                parts = line.split()
                if len(parts) < 3:
                    continue
                atom1_idx = int(parts[1]) - 1  # Convert to 0-based index
                atom2_idx = int(parts[2]) - 1  # Convert to 0-based index
                bond_type = parts[3]
                # Determine bond order based on bond type
                if bond_type == "1":
                    bond_order = 1.0  # Single bond
                elif bond_type == "2":
                    bond_order = 2.0  # Double bond
                elif bond_type == "3":
                    bond_order = 3.0  # Triple bond
                elif bond_type == "ar":
                    bond_order = 1.5  # Aromatic bond
                else:
                    bond_order = 0.0  # Undefined or unsupported bond type; ignore
                self.add_bond(atom1_idx, atom2_idx, bond_type)
        return self

    
    def writefile(self, file_format, outputgeom):
        """
        Write geometry data to a file.

        Parameters:
        - outputgeom: Path to the output geometry file.
        - file_format: Format of the file (e.g., 'xyz', 'mol2').
        """
        with open(outputgeom, 'w') as file:
            content = self._parse_write(file_format)
            file.write(content)
        
    def writestr(self, file_format):
        """
        Generate geometry data as a string in the specified format.

        Parameters:
        - file_format: Format of the data (e.g., 'xyz', 'mol2').

        Returns:
        - A string containing the geometry data.
        """
        content = self._parse_write(file_format)
        return content

    def _parse_write(self, file_format):
        """
        Generate the geometry content based on the specified file format.

        Parameters:
        - file_format: Format of the data (e.g., 'xyz', 'mol2').

        Returns:
        - The content string in the specified format.
        """
        if file_format == 'xyz':
            return self.writexyz()
        elif file_format == 'mol2':
            return self.writemol2()
        else:
            raise ValueError(f"Unsupported format: '{file_format}'")

    def writexyz(self):
        """
        Generate XYZ format content for the geometry.

        Returns:
        - A string containing the XYZ format of the geometry.
        """
        num_atoms = len(self.atoms)
        lines = [f"{num_atoms}", ""]  # XYZ format starts with number of atoms and a comment line
        for atom in self.atoms:
            lines.append(f"{atom.sym}\t{atom.coord[0]:.6f}\t{atom.coord[1]:.6f}\t{atom.coord[2]:.6f}")
        return "\n".join(lines)

    def writemol2(self):
        """
        Generate MOL2 format content for the geometry.

        Returns:
        - A string containing the MOL2 format of the geometry.
        """
        lines = []
        
        # MOLECULE
        lines.append('@<TRIPOS>MOLECULE')
        mol_name = 'hi'
        lines.append(mol_name)
        num_atoms = self.natoms
        num_bonds = len(self.bonds)  # Use the length of self.bonds
        num_substr = 1
        lines.append(f'{num_atoms} {num_bonds} {num_substr}')
        mol_type = 'SMALL'
        lines.append(mol_type)
        charge_type = 'PartialCharges'
        lines.append(charge_type)
        lines.append('****')
        comment = 'catoms'
        lines.append(comment)
        lines.append('')
        
        # ATOM
        lines.append("@<TRIPOS>ATOM")
        for atom in self.atoms:
            lines.append(f"{atom.idx+1} {atom.sym} {atom.coord[0]:.6f} {atom.coord[1]:.6f} {atom.coord[2]:.6f} {atom.sym} 1")
        
        # BOND
        lines.append("@<TRIPOS>BOND")
        for bond_idx, (atom_pair, bond_type) in enumerate(self.bonds.items(), start=1):
            atom1, atom2 = atom_pair
            lines.append(f"{bond_idx} {atom1+1} {atom2+1} {bond_type}")
        
        return "\n".join(lines)
    
    def xyztoOBMol(self):
        OB_Conversion = openbabel.OBConversion()
        OB_Conversion.SetInFormat('xyz')
        OB_Mol = openbabel.OBMol()
        OB_Conversion.ReadString(OB_Mol, self.writexyz())
        return OB_Mol

    def get_bond_OB(self, savetoself=False):
        OB_Mol = self.xyztoOBMol()
        bonds = dict()
        for bond in openbabel.OBMolBondIter(OB_Mol):
            atomA = bond.GetBeginAtom()
            atomB = bond.GetEndAtom()
            bond_order = bond.GetBondOrder()
            is_aromatic = bond.IsAromatic()
            atomA_idx = atomA.GetIdx() - 1
            atomB_idx = atomB.GetIdx() - 1
            if is_aromatic:
                bond_type = 'ar'
            elif bond_order == 1:
                bond_type = '1'
            elif bond_order == 2:
                bond_type = '2'
            elif bond_order == 3:
                bond_type = '3'
            else:
                bond_type = 'unknown'
                
            bond_key = tuple(sorted((atomA_idx, atomB_idx)))
            bonds[bond_key] = bond_type
            
            if savetoself:
                self.add_bond(atomA_idx, atomB_idx, bond_type)

        return bonds
    
    def get_atom_hybridizations(self):
        OB_Mol = self.xyztoOBMol()
        hybridizations = {}
        for atom in openbabel.OBMolAtomIter(OB_Mol):
            atom_index = atom.GetIdx() - 1  # Zero-based indexing
            hyb = atom.GetHyb()  # Get hybridization
            if hyb == 1:
                hyb_str = 'sp'
            elif hyb == 2:
                hyb_str = 'sp2'
            elif hyb == 3:
                hyb_str = 'sp3'
            else:
                hyb_str = 'other'  # For any hybridization not covered
            hybridizations[atom_index] = hyb_str
        return hybridizations
            
    def get_coord_byidx(self, idx):
        return np.array(self.atoms[idx].coord)
    
    def atom_byidx(self, idx):
        return self.atoms[idx]
    
    def get_listofatomprop(self, prop):
        """
        Collects the specified property from each Atom object in self.atoms.

        Parameters
        ----------
        prop : str
            The property name to retrieve from each Atom object.

        Returns
        -------
        np.ndarray
            An array containing the values of the specified property from each Atom.

        Raises
        ------
        ValueError
            If an Atom does not have the specified property.
        """
        prop_list = []  # Collect the specified properties
        for atom in self.atoms:
            try:
                # Access the property or attribute directly
                prop_value = getattr(atom, prop)
                prop_list.append(prop_value)
            except AttributeError:
                raise ValueError(f"Atom object has no property '{prop}'")
        
        return np.array(prop_list)  # Return the list as a NumPy array
    
    def _build_xyz_df(self):
        """
        Build a DataFrame that contains atom indices and their coordinates.
        """
        data = {
            'atom_idx': [atom.idx for atom in self.atoms],
            'element': [atom.sym for atom in self.atoms],
            'x': [atom.coord[0] for atom in self.atoms],
            'y': [atom.coord[1] for atom in self.atoms],
            'z': [atom.coord[2] for atom in self.atoms],
        }
        self.xyz_df = pd.DataFrame(data)

    def get_distance_matrix(self):
        """
        Calculate and return the distance matrix for the molecule.

        Returns
        -------
        distance_matrix : np.ndarray
            A 2D numpy array where distance_matrix[i, j] is the distance between atom i and atom j.
        """
        # Number of atoms in the molecule
        num_atoms = len(self.atoms)

        # Initialize an empty distance matrix
        distance_matrix = np.zeros((num_atoms, num_atoms))

        # Calculate pairwise distances using calc_distance
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):  # Only calculate the upper triangle
                distance = self.calc_distance(i, j)
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance  # Symmetric matrix

        return distance_matrix
    
    def get_bond_existed(self, cutoff=0.1):
        """
        Determine existing bonds based on the covalent radii of atoms.

        Parameters
        ----------
        cutoff : float, optional
            Additional distance tolerance added to the sum of covalent radii to consider a bond.

        Returns
        -------
        bond_existed : list
            A list of bonds, where each bond is represented by (atom_i, atom_j, bond_length).
        """
        distance_matrix = self.get_distance_matrix()
        num_atoms = len(self.atoms)
        
        bond_existed = []

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                bond_length = distance_matrix[i, j]
                cov_radii_i = self.atoms[i].cov_radii
                cov_radii_j = self.atoms[j].cov_radii
                covalent_sum = cov_radii_i + cov_radii_j + cutoff
                if bond_length < covalent_sum:
                    bond_existed.append((self.atoms[i], self.atoms[j], bond_length))

        return bond_existed
    
    def calc_distance(self, atomA_idx, atomB_idx):
        try:
            atomA_coord = self.get_coord_byidx(atomA_idx)
            atomB_coord = self.get_coord_byidx(atomB_idx)
        except KeyError as e:
            raise ValueError(f"Invalid atom index: {e}")
        dist = np.linalg.norm(atomA_coord - atomB_coord)
        return dist
    
    def calc_angle(self, atomA_idx, atomB_idx, atomC_idx):
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
        return angle_degrees
    
    def calc_dihedral(self, atomA_idx, atomB_idx, atomC_idx, atomD_idx):
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
        list_indices = []
        for atom in self.atoms:
            if atom.sym == element:
                list_indices.append(atom.idx)
        # list_indices = [3, 4, 14, 15, 16, 19, 21, 24, 28, 32] of the element
        return list_indices
    
    def atomnearby_from_idx(self, atomA_idx, cutoff=0.1):
        """
        Find atoms nearby a specified atom based on covalent radii and a distance cutoff.

        Parameters
        ----------
        atomA_idx : int
            The index of the atom to find nearby atoms for.
        cutoff : float, optional
            Additional distance tolerance added to the sum of covalent radii to consider bonding.

        Returns
        -------
        nearby_atoms : set
            A set of Atom instances that are within bonding distance from the specified atom.
        """
        bond_existed = self.get_bond_existed(cutoff)
        nearby_atoms = set()
        for atom_i, atom_j, bond_length in bond_existed:
            if atomA_idx == atom_i.idx:
                nearby_atoms.add(atom_j)
            elif atomA_idx == atom_j.idx:
                nearby_atoms.add(atom_i)
        return nearby_atoms
    
    def get_unique_elements(self):
        unique_elements = {atom.sym for atom in self.atoms}
        return unique_elements
    
    def get_transition_metal(self):
        unique_elements = self.get_unique_elements()
        metals = unique_elements.intersection(TRANSITION_METALS)
        if len(metals) == 0:
            return None
        return metals  # Set of metals, e.g., {'Ni', 'Pd'}
    
    def getAtom_transition_metal(self):
        atom_list = []
        for atom in self.atoms:
            if atom.sym in TRANSITION_METALS:
                atom_list.append(atom)
        return atom_list
    
    def isbond(self, atomA_idx, atomB_idx, cutoff=0.1):
        cov_radii_A = self.atoms[atomA_idx].cov_radii
        cov_radii_B = self.atoms[atomB_idx].cov_radii
        distance = self.calc_distance(atomA_idx, atomB_idx)
        isbond = distance < (cov_radii_A + cov_radii_B + cutoff)
        return isbond

    def sum_cov_radii(self, elementA, elementB, cutoff=0.0):
        cov_radii_A = COVALENT_RADII[elementA]
        cov_radii_B = COVALENT_RADII[elementB]
        sum_radii = cov_radii_A + cov_radii_B + cutoff
        return sum_radii

    def find_ligand_mononuclear(self, getonlyidx=False):
        metal_element = list(self.get_transition_metal())[0]
        metal_idx = self.atomidx_from_element(metal_element)[0]
        bondtometal = self.atomnearby_from_idx(metal_idx)
        if getonlyidx:
            bondtometal = {atom.idx for atom in bondtometal}
        return bondtometal

    def find_nearest_ofelement_toatom(self, center_atom_index, target_element, cutoff=0.1):
        bond_to_center_atom = self.atomnearby_from_idx(center_atom_index, cutoff)
        atoms_with_target_element = [atom for atom in bond_to_center_atom if atom.sym == target_element]
        nearest_atom = min(atoms_with_target_element, key=lambda atom: self.calc_distance(center_atom_index, atom.idx))
        return nearest_atom

    def is_planar(self, list_of_atoms, tolerance=5):
        if len(list_of_atoms) < 4:
            raise ValueError("At least 4 atom indices are required to check for planarity.")
        atomA_idx, atomB_idx, atomC_idx = list_of_atoms[0], list_of_atoms[1], list_of_atoms[2]
        for i in range(3, len(list_of_atoms)):
            atomD_idx = list_of_atoms[i]
            dihedral = self.calc_dihedral(atomA_idx, atomB_idx, atomC_idx, atomD_idx)
            if abs(dihedral) > tolerance and abs(dihedral - 180) > tolerance:
                return False
        return True

    def is_sqp(self, tolerance=5):
        if len(self.find_ligand_mononuclear()) != 4:
            raise ValueError("Square planar geometry requires exactly 4 ligands.")
        ligatoms = list(self.find_ligand_mononuclear())
        atomA_idx, atomB_idx, atomC_idx, atomD_idx = ligatoms[0].idx, ligatoms[1].idx, ligatoms[2].idx, ligatoms[3].idx
        angle = self.calc_dihedral(atomA_idx, atomB_idx, atomC_idx, atomD_idx)
        return abs(angle) < tolerance or abs(angle - 180) < tolerance

    def find_normv(self, list_of_atoms):
        if len(list_of_atoms) != 3:
            raise ValueError('You need exactly 3 atoms to find the normal vector of a plane.')
        atomA_idx, atomB_idx, atomC_idx = list_of_atoms
        atomA_coord = self.get_coord_byidx(atomA_idx)
        atomB_coord = self.get_coord_byidx(atomB_idx)
        atomC_coord = self.get_coord_byidx(atomC_idx)
        vector_BA = atomA_coord - atomB_coord
        vector_BC = atomC_coord - atomB_coord
        normv = np.cross(vector_BA, vector_BC)
        normv = normv / np.linalg.norm(normv)
        if normv[0] < 0:
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
        return np.degrees(angle_radians)

    def reindex_atoms_and_bonds(self):
        """
        Reindex all atoms in self.atoms and update self.bonds to reflect new atom indices.
        """
        # Create a mapping from old indices to new indices
        old_to_new_indices = {}

        # Reindex atoms
        for new_idx, atom in enumerate(self.atoms):
            old_to_new_indices[atom.idx] = new_idx
            atom.idx = new_idx  # Update the atom's index

        # Reindex bonds based on updated atom indices
        new_bonds = {}
        for (old_atomA_idx, old_atomB_idx), bond_type in self.bonds.items():
            new_atomA_idx = old_to_new_indices[old_atomA_idx]
            new_atomB_idx = old_to_new_indices[old_atomB_idx]
            bond_key = tuple(sorted((new_atomA_idx, new_atomB_idx)))
            new_bonds[bond_key] = bond_type

        # Replace self.bonds with the reindexed bonds
        self.bonds = new_bonds
    
    def add_atom(self, atom, reindex=True):
        """
        Add a new atom to the molecule.

        Parameters
        ----------
        atom : Atom
            The atom to be added.
        reindex : bool, optional
            If True, reindex all atoms after adding. Defaults to True.
        """
        # Add atom with updated index
        atom.idx = len(self.atoms)
        self.atoms.append(atom)
        
        # Reindex if necessary
        if reindex:
            self.reindex_atoms_and_bonds()
            
    def addMol(self, submol, reindex=True):
        """
        Add all atoms and bond connections from submol to the current molecule.

        Parameters
        ----------
        submol : Mol
            The molecule to be added.
        reindex : bool, optional
            If True, reindex all atoms after adding. Defaults to True.
        """
        # Calculate the index offset for new atoms
        index_offset = len(self.atoms)

        # Step 1: Add atoms from submol to self.atoms, with adjusted indices
        new_atoms = []
        for atom in submol.atoms:
            # Create a new Atom with an updated index and add it to self.atoms
            new_index = atom.idx + index_offset
            new_atom = Atom(idx=new_index, sym=atom.sym, coord=atom.coord)
            new_atoms.append(new_atom)  # Append directly to self.atoms
        new_atoms_array = np.array(new_atoms, dtype=object)
        self.atoms = np.append(self.atoms, new_atoms_array)
        
        # Step 2: Add bonds from submol to self.bonds, with adjusted indices
        new_bonds = dict()
        for (atomA_idx, atomB_idx), bond_type in submol.bonds.items():
            # Adjust indices by adding the offset
            new_atomA_idx = atomA_idx + index_offset
            new_atomB_idx = atomB_idx + index_offset

            # Use the adjusted indices to add the bond to self.bonds
            bond_key = tuple(sorted((new_atomA_idx, new_atomB_idx)))
            new_bonds[bond_key] = bond_type
        self.bonds.update(new_bonds)
        
        # Step 3: Reindex atoms and bonds in self if requested
        if reindex:
            self.reindex_atoms_and_bonds()

    def reset_index_atoms_bonds(self):
        """
        Reindex all atoms in self.atoms and update self.bonds to reflect new atom indices.
        """
        # Create a mapping from old indices to new indices
        old_to_new_indices = {}

        # Reindex atoms
        for new_idx, atom in enumerate(self.atoms):
            old_to_new_indices[atom.idx] = new_idx
            atom.idx = new_idx  # Update the atom's index

        # Reindex bonds based on updated atom indices
        new_bonds = {}
        for (old_atomA_idx, old_atomB_idx), bond_type in self.bonds.items():
            new_atomA_idx = old_to_new_indices[old_atomA_idx]
            new_atomB_idx = old_to_new_indices[old_atomB_idx]
            bond_key = tuple(sorted((new_atomA_idx, new_atomB_idx)))
            new_bonds[bond_key] = bond_type

        # Replace self.bonds with the reindexed bonds
        self.bonds = new_bonds
        
    def deleteAtom(self, atom_idx, reindex=True):
        """
        Delete an atom and its associated bonds from the molecule.

        Parameters
        ----------
        atom_idx : int
            The index of the atom to delete.
        reindex : bool, optional
            If True, reindex all atoms and bonds after deletion. Defaults to True.
        """
        # Remove the atom from self.atoms
        self.atoms = [atom for atom in self.atoms if atom.idx != atom_idx]

        # Remove any bonds involving this atom
        self.bonds = {
            bond_key: bond_type
            for bond_key, bond_type in self.bonds.items()
            if atom_idx not in bond_key
        }

        # Reindex if necessary
        if reindex:
            self.reindex_atoms_and_bonds()
    
    def deleteAtoms(self, indices, reindex=True):
        """
        Deletes a list of Atoms from the molecule based on their indices.

        Parameters
        ----------
        indices : list of int
            The list of indices of Atoms to delete.

        Raises
        ------
        ValueError
            If any Atom with a specified index is not found.
        """
        indices_set = set(indices)
        self.atoms = [atom for atom in self.atoms if atom.idx not in indices_set]
        if reindex:
            self.reindex_atoms_and_bonds()
    
    def deleteAtoms_bysym(self, symbols):
        """
        Deletes all atoms with specified element symbols from the molecule.

        Parameters
        ----------
        symbols : list of str
            List of element symbols (e.g., ['Fe', 'X']) to delete from the molecule.

        Raises
        ------
        ValueError
            If no atoms with the specified symbols are found.
        """
        if not isinstance(symbols, list):
            raise ValueError("The symbols parameter must be a list of element symbols (strings).")
        symbols_set = set(symbols)
        atoms_to_keep = [atom for atom in self.atoms if atom.sym not in symbols_set]
        if len(atoms_to_keep) == len(self.atoms):
            raise ValueError(f"No atoms with symbols {symbols} found.")
        self.atoms = atoms_to_keep
        self._coords = np.array([atom.coord for atom in self.atoms]) if self.atoms else None
        self.reindex_atoms_and_bonds()
        
    def deleteHs(self):
        """
        Remove all hydrogen atoms and any associated bonds from the molecule.
        """
        hydrogen_indices = {atom.idx for atom in self.atoms if atom.sym == "H"}
        self.atoms = [atom for atom in self.atoms if atom.idx not in hydrogen_indices]
        self.bonds = {
            bond_key: bond_type
            for bond_key, bond_type in self.bonds.items()
            if bond_key[0] not in hydrogen_indices and bond_key[1] not in hydrogen_indices
        }
        self.reindex_atoms_and_bonds()
    
        
    def add_bond(self, atomA_idx, atomB_idx, bond_type):
        """
        Adds a bond between two atoms with a specified bond type.

        Parameters
        ----------
        atomA : Atom
            The first atom in the bond.
        atomB : Atom
            The second atom in the bond.
        bond_type : str
            The type of bond (e.g., '1', '2', 'ar').

        Returns
        -------
        self : Mol
            Returns the instance to allow method chaining.
        """
        bond_key = tuple(sorted((atomA_idx, atomB_idx)))
        self.bonds[bond_key] = bond_type
        sorted_bonds = dict(sorted(self.bonds.items(), key=lambda x: x[0][0]))
        self.bonds = sorted_bonds
        
        return self
    
    def delete_bond(self, atomA_idx, atomB_idx):
        """
        Deletes a bond between two atoms if it exists.

        Parameters
        ----------
        atomA_idx : int
            The index of the first atom in the bond.
        atomB_idx : int
            The index of the second atom in the bond.

        Returns
        -------
        bool
            Returns True if the bond was successfully deleted, False if the bond was not found.
        """
        bond_key = tuple(sorted((atomA_idx, atomB_idx)))
        if bond_key in self.bonds:
            del self.bonds[bond_key]
            return True
        else:
            return False
    

    def get_bond_fromidx(self, atom_idx):
        """
        Retrieve all bonds involving a specific atom index.

        Parameters
        ----------
        atom_idx : int
            The index of the atom for which to retrieve bonds.

        Returns
        -------
        list of tuples
            A list of tuples where each tuple contains the other atom's index
            and the bond type.
        """
        bonds_for_atom = dict()

        # Iterate through self.bonds and check if atom_idx is in the bond key
        for (atomA_idx, atomB_idx), bond_type in self.bonds.items():
            if atomA_idx == atom_idx:
                # If atom_idx is the first atom in the bond
                bond = {(atomA_idx, atomB_idx): bond_type}
                bonds_for_atom.update(bond)
            elif atomB_idx == atom_idx:
                # If atom_idx is the second atom in the bond
                bond = {(atomA_idx, atomB_idx): bond_type}
                bonds_for_atom.update(bond)

        return bonds_for_atom
    
    def get_bond_type(self, atomA_idx, atomB_idx):
        """
        Retrieves the bond type between two atoms, if it exists.

        Parameters
        ----------
        atom_a : Atom or int
            The first atom in the bond, or its index.
        atom_b : Atom or int
            The second atom in the bond, or its index.

        Returns
        -------
        str or None
            The bond type (e.g., '1', '2', 'ar') if the bond exists, or None if it doesn't.
        """
        bond_key = tuple(sorted((atomA_idx, atomB_idx)))
        return self.bonds.get(bond_key, None)
    
    def get_bond_length(self, atomA_idx, atomB_idx):
        """
        Get the bond length between two atoms specified by their indices.

        Parameters
        ----------
        atomA_idx : int
            The index of the first atom.
        atomB_idx : int
            The index of the second atom.

        Returns
        -------
        float
            The bond length (distance) between the two atoms.

        Raises
        ------
        ValueError
            If there is no bond between atomA_idx and atomB_idx in self.bonds.
        """
        # Ensure that the bond exists in self.bonds
        bond_key = tuple(sorted((atomA_idx, atomB_idx)))
        if bond_key not in self.bonds:
            raise ValueError(f"No bond found between atoms {atomA_idx} and {atomB_idx}")

        # Calculate and return the bond length using Euclidean distance
        return self.calc_distance(atomA_idx, atomB_idx)

    def getAtom_fromidx(self, idx):
        matching_atoms = [atom for atom in self.atoms if atom.idx == idx]
        if len(matching_atoms) > 1:
            raise ValueError(f"More than one atom found with index {idx}.")
        elif len(matching_atoms) == 0:
            raise ValueError(f"No atom found with index {idx}.")
        return matching_atoms[0]
    
    def getlistofsym(self):
        list_sym = []
        for atom in self.atoms:
            sym = atom.sym
            list_sym.append(sym)
        return list_sym
    
    def getcoord_fromatomlist(self, list_Atoms):
        """
        Get the coordinates of atoms in the provided list.

        Parameters
        ----------
        list_Atoms : list
            List of Atom objects from which to extract coordinates.

        Returns
        -------
        np.ndarray
            A 2D numpy array with each row representing the coordinates of an atom.
        
        Raises
        ------
        ValueError
            If the input is not a list.
        """
        if not isinstance(list_Atoms, list):
            raise ValueError("The input must be a list of Atom objects.")
        
        coords = np.vstack([atom.coord for atom in list_Atoms])
        return coords
    
    def generate_MultiGraph(self):
        """
        Generate and return a NetworkX MultiGraph from self.atoms and self.bonds.

        Returns
        -------
        G : networkx.MultiGraph
            The generated molecular graph where nodes represent atoms and edges represent bonds.
        """
        G = nx.MultiGraph()

        for atom in self.atoms:
            G.add_node(atom.idx, sym=atom.sym, coord=atom.coord)
        
        for (atomA_idx, atomB_idx), bond_type in self.bonds.items():
            G.add_edge(atomA_idx, atomB_idx, bond_type=bond_type)

        return G
    
    def generate_Graph(self):
        """
        Generate and return a NetworkX Graph from self.atoms and self.bonds.

        Returns
        -------
        G : networkx.Graph
            The generated molecular graph where nodes represent atoms and edges represent bonds.
            Only a single edge is allowed between any two nodes.
        """
        G = nx.Graph()

        # Add atoms as nodes
        for atom in self.atoms:
            G.add_node(atom.idx, sym=atom.sym, coord=atom.coord)

        for (atomA_idx, atomB_idx), bond_type in self.bonds.items():
            # If a bond already exists, it will not overwrite in an nx.Graph
            if not G.has_edge(atomA_idx, atomB_idx):
                G.add_edge(atomA_idx, atomB_idx, bond_type=bond_type)

        return G
    
    def get_graph_hash(self):
        """
        Generate the Weisfeiler-Lehman graph hash for the full molecule.

        Returns
        -------
        gh : str
            The Weisfeiler-Lehman graph hash of the molecule.
        """
        graph = self.generate_Graph()
        graph_hash = nx.weisfeiler_lehman_graph_hash(graph, node_attr='sym')
        return graph_hash
    
    def generate_bond_matrix(self):
        """
        Generate a bond (adjacency) matrix for the molecule.

        Returns
        -------
        bond_matrix : np.ndarray
            A 2D numpy array where bond_matrix[i, j] represents the bond type or bond order
            between atom i and atom j.
        """
        num_atoms = self.natoms
        bond_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        for (atomA_idx, atomB_idx), bond_type in self.bonds.items():
            # Set the bond order or type in both directions (i, j) and (j, i)
            bond_order = self._get_bond_order(bond_type)
            bond_matrix[atomA_idx, atomB_idx] = bond_order
            bond_matrix[atomB_idx, atomA_idx] = bond_order
        return bond_matrix

    def _get_bond_order(self, bond_type):
        """
        Convert bond type to bond order.

        Parameters
        ----------
        bond_type : str
            The type of bond (e.g., 'single', 'double', 'triple', 'ar').

        Returns
        -------
        int
            The bond order as an integer (e.g., 1 for single, 2 for double, etc.).
        """
        bond_order_map = {'1': 1.0,
                          '2': 2.0,
                          '3': 3.0,
                          'ar': 1.5,
                          'ts': 0.5
        }
        return bond_order_map.get(bond_type, 0)  # Default to 0 if bond type is unrecognized
    
    def calc_buriedV(self, radius=3.5, meshspace=0.1, bondiscale=1.17, hydrogens=False):
        # Delete hydrogens if specified
        if not hydrogens:
            self.deleteHs()
    
    def set_coordinate(self, index, new_coordinate):
        """
        Update the coordinates of a specific atom in the molecule.

        Args:
            index (int): The index of the atom to update.
            new_coordinate (tuple[float, float, float]): The new (x, y, z) coordinates for the atom.

        Raises:
            IndexError: If the provided index is out of range.
            ValueError: If `new_coordinate` is not a tuple of length 3.

        Example:
            >>> mol = Mol(['C', 'H'], [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)])
            >>> mol.set_coordinate(1, (1.5, 1.5, 1.5))
        """
        if index < 0 or index >= len(self.coords):
            raise IndexError("Index out of range for atom coordinates.")
        if len(new_coordinate) != 3:
            raise ValueError("Coordinates must be a tuple of 3 values (x, y, z).")

        self.coords[index] = tuple(new_coordinate)
        print(f"Updated atom at index {index} to new coordinates: {self.coords[index]}")
        
    def set_origin_by_index(self, index):
        """
        Translate all coordinates such that the atom at the given index becomes the origin.

        Args:
            index (int): Index of the atom to set as the origin.

        Raises:
            IndexError: If the index is out of range.
        """
        if index < 0 or index >= self.natoms:
            raise IndexError("Index out of range for atom coordinates.")

        # Get the coordinate of the atom at the given index
        origin = self.coords[index]

        # Translate all coordinates
        self.coords -= origin
        print(f"Translated all coordinates so atom at index {index} is now at the origin.")
        # Find the metal atom and center the coordinates around it
    
    # def adjust_bondlength(self, atom_fix_idx, atom_scaled_idx, scaled_distance):
    #     atom_fix_coord = self.coords[atom_fix_idx]
    #     atom_scaled_coord = self.coords[atom_scaled_idx]
    #     new_atom_scaled_coord = 
    
    
class MolLigand(Mol):
    def __init__(self):
        """
        Initialize the MolLigand object with optional attributes for name, path, and SMILES.
        
        Parameters
        ----------
        name : str, optional
            Name of the molecule.
        path : str, optional
            File path for the structure, if specified.
        smi : str, optional
            SMILES string for the molecule, only when needed.
        """
        super().__init__()
        self.catoms = []  # List of critical atoms, input from the user
        self.name = None
        self.smi = None  # SMILES string if specified
    
    def __repr__(self):
        return f"MolLigand(natoms={self.natoms}, catoms={self.catoms})"
    
    def set_catoms(self, catoms):
        """
        Set the critical atoms (catoms) for the molecule.
        
        Parameters
        ----------
        catoms : list
            List of atom indices or identifiers provided by the user.
        """
        if not isinstance(catoms, list):
            raise ValueError("catoms must be a list.")
        self.catoms = catoms
    
    def set_name(self, name):
        """
        Set the name of the molecule.
        
        Parameters
        ----------
        name : str
            The name to assign to the molecule.
        """
        if not isinstance(name, str):
            raise ValueError("The name must be a string.")
        self.name = name
        
    def getAtom_catoms(self):
        """
        Get the Atom instances corresponding to the indices in catoms.

        Returns
        -------
        List[Atom]
            List of Atom instances for the critical atoms (catoms).

        Raises
        ------
        ValueError
            If catoms is not set or if any index in catoms does not correspond to an atom in atoms.
        """
        if self.catoms is None or len(self.catoms) == 0:
            raise ValueError("catoms must be set before calling getcatoms_Atom.")

        catom_atoms = [atom for atom in self.atoms if atom.idx in self.catoms]

        # Ensure all catoms indices have corresponding Atom objects
        if len(catom_atoms) != len(self.catoms):
            raise ValueError("Not all indices in catoms correspond to atoms in the molecule.")

        return set(catom_atoms)
            
    def getAtom_adjcatom(self):
        """
        Get a dictionary of critical atoms (catoms) and their adjacent atoms.

        Returns
        -------
        dict
            A dictionary where each key is the index of a critical atom (catom) and the value
            is a list of adjacent atoms.

        Raises
        ------
        ValueError
            If catoms is not set or if any index in catoms does not correspond to an atom in atoms.
        """
        if not self.catoms:
            raise ValueError("catoms must be set before calling getAtom_adjcatom.")
        
        adjcatoms = {}
        for idx_catom in self.catoms:
            atomnearby = list(self.atomnearby_from_idx(idx_catom))
            adjcatoms[idx_catom] = atomnearby
        return adjcatoms
    
    # def lonepair_vec(self):
    
    
        
        
        
