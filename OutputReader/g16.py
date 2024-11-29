import numpy as np
import pandas as pd
from ..globalvars import ATOMIC_NUM_TO_ELEMENT, COVALENT_RADII, ELEMENT_TO_ATOMIC_NUM

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

class MainReader:
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
    def __init__(self, logfile):
        self.logfile = logfile
        self.job_finished = None
        self.runtime = None
        self.job = []
        self.natoms = 0
        self.xyz_comment = ''
        self.xyz_coord = []
        self.xyz_df = None
        self.distance_matrix = None
        self.xyz_string = str()
        self

        with open(logfile, 'r') as file:
            self.lines_forward = file.readlines()
        self.lines_num = len(self.lines_forward)
        self.lines_reverse = self.lines_forward[::-1]
        
        self.get_natoms()
        self.chk_jobfinish()
        # self.get_xyzdata()
        # self.create_xyz_df()
        # self.create_distance_matrix()
        
    def stoichiometry(self):
        lines = self.lines_forward
        for idx, line in enumerate(lines):
            if 'Stoichiometry' in line:
                formula = line.split()[1]
                break
        element_counts = {}
        i = 0
        while i < len(formula):
            # Check if current character is an uppercase letter (start of an element symbol)
            if formula[i].isupper():
                element = formula[i]  # Start of the element symbol
                # Check if the next character is a lowercase letter (second part of the element symbol)
                if i + 1 < len(formula) and formula[i + 1].islower():
                    element += formula[i + 1]
                    i += 1  # Move index to skip the lowercase letter
                # Move to the next character (which could be a number or another element)
                i += 1
                # Now we need to check if there is a number after the element symbol
                count = 0
                while i < len(formula) and formula[i].isdigit():
                    count = count * 10 + int(formula[i])
                    i += 1
                # If no number was found, default to 1
                if count == 0:
                    count = 1
                # Store the element and its count in the dictionary
                element_counts[element] = count
        return element_counts

    def get_natoms(self):
        element_counts = self.stoichiometry()
        self.natoms = sum(element_counts.values())
        return self.natoms
    
    def inputline(self):
        input_appearance = 0
        input_reading = False
        list_input = []
        input_line = ''  # Initialize input_line outside the loop

        for idx, line in enumerate(self.lines_forward):
            if line.startswith(' #'):
                if input_line:  # If there is something in input_line, append it to the list before resetting
                    list_input.append(input_line)
                input_line = ''  # Reset input_line for the new section
                input_appearance += 1
                input_reading = True
            elif '------' in line:
                input_reading = False
                if input_line:  # If we reached the end of a section, append the input line
                    list_input.append(input_line)
                input_line = ''  # Reset for the next section
            if input_reading:
                input_line += line.strip()  # Continue reading and appending the line

        # After the loop, ensure that any remaining input_line is added to list_input
        if input_line:
            list_input.append(input_line)

        return list_input  # Return the list of input sections
            
    def jobtype(self):
        list_input = self.inputline()
        self.job = []
        
        for section in list_input:
            section_lower = section.lower()
            if 'irc' in section_lower:
                self.job.append('irc')
            elif 'opt' in section_lower and 'freq' in section_lower:
                self.job.append('optfreq')
            elif 'opt' in section_lower:
                self.job.append('opt')
            elif 'freq' in section_lower:
                self.job.append('freq')
            else:
                self.job.append('energy')
        return self.job
    
    def chk_jobfinish(self):
        lines = self.lines_forward
        self.job_finished = [False] * len(self.inputline()) 
        occurence = 0
        for idx, line in enumerate(lines):
            if 'Normal termination' in line:
                if occurence < len(self.job_finished):
                    self.job_finished[occurence] = True
                occurence += 1 
        return self.job_finished
    
    def chk_runtime(self, showin='sec'):
        lines = self.lines_reverse
        runtimes = []
        for line in lines:
            if 'Elapsed time' in line:
                time_str = line
                time_components = time_str.split()
                days = int(time_components[2])
                hours = int(time_components[4])
                minutes = int(time_components[6])
                seconds = float(time_components[8])
                total_dhms = [days, hours, minutes, seconds]
                total_sec = days * 86400 + hours * 3600 + minutes * 60 + seconds
                if showin == 'sec':
                    runtimes.append(total_sec)
                elif showin == 'dhms':
                    runtimes.append(total_dhms)
        return runtimes  # Return the list of all runtimes
    
    def coord_xyz_endjob(self, xyzfiletowrite=None, writexyzfile=False):
        coordinates = []
        reading = False
        lines = self.lines_forward
        for line in lines:
            if 'Unable to Open any file for archive entry.' in line:
                print('hi')
                reading = True
                continue
            elif 'The archive entry for this job was punched.' in line:
                reading = False
                continue
            if reading:
                coordinates.append(line.strip())
        coordinates = ' '.join(coordinates)
        coordinates = coordinates.replace(' ','')
        coordinates = coordinates.split(r'\\')[3]
        coordinates = coordinates.replace("\\","\n")
        coordinates = coordinates[4:]
        coordinates = coordinates.replace(',','\t')
        xyz_coord_lines = coordinates.strip().split('\n')
        xyz_num = len(xyz_coord_lines)
        xyz_comment = ''
        xyz_coord = [line for line in xyz_coord_lines if line]
        
        if writexyzfile:
            if xyzfiletowrite is None:
                raise ValueError("name of xyzfile must be provided if writexyzfile=True")
            with open(xyzfiletowrite, 'w') as f:
                f.write(f'{xyz_num}\n')
                f.write(f'{xyz_comment}\n')
                for coord in xyz_coord:
                    f.write(f'{coord}\n')
        return xyz_num, xyz_comment, xyz_coord
          
    def keyword_occurence(self, keyword):
        """ Check if a keyword is found in a file

        Args:
            keyword: str
                keyword
        """
        occurence = 0
        occur_at = []
        for idx, line in enumerate(self.lines_forward):
            if keyword in line:
                occurence += 1
                occur_at.append(idx)
        return occurence, occur_at
    
    def read_orientation_table(self):
        natoms = self.natoms
        lines = self.lines_forward
        line_indices = self.keyword_occurence('Standard orientation')[1]
        
        if not line_indices:
            print("Error: 'Standard orientation' not found")
            return

        all_coords = []  # Store all coordinates

        for idx in line_indices:
            table = lines[idx + 5 : idx + 5 + natoms]
            coord_thisoccur = self._read_single_orientation_table(table)
            
            all_coords.append(coord_thisoccur)
            print(coord_thisoccur)

        return all_coords  # Return all coordinates
    
    def _read_single_orientation_table(self, table):
        # idx = starting index line of the orientation table, index from forward lines
        coord_thistable = []
        for row in table:
            coord_eachrow = row.split()

            if len(coord_eachrow) < 6:
                print(f"Error: Unexpected row format: {row}")
                continue

            atom_idx = int(coord_eachrow[0]) - 1
            atomic_num = int(coord_eachrow[1])
            atom_element = ATOMIC_NUM_TO_ELEMENT.get(atomic_num, 'Unknown')
            atom_x = float(coord_eachrow[-3])
            atom_y = float(coord_eachrow[-2])
            atom_z = float(coord_eachrow[-1])
            coord_thistable.append([atom_idx, atomic_num, atom_element, atom_x, atom_y, atom_z])
        return coord_thistable
    
    def get_xyzdata_inputorien(self, comment=''):
        """ Extract atomic coordinates from the log file.
            
        Attributes updated:
        ------------
            self.xyz_coord: list
                List of atomic coordinates [atom_idx, atomic_num, element, x, y, z].
                len(self.xyz_coord) = number of atoms.
        """
        for idx, line in enumerate(self.lines_reverse):
            if 'Input orientation' in line:
                idx_coord_begin = self.lines_num - 1 - idx
                break
        else:
            raise Exception("Input orientation section not found in the log file.")
        
        for idx, line in enumerate(self.lines_forward[idx_coord_begin + 5:]):
            if not '----' in line:
                coord_eachrow = line.split()
                atom_idx = int(coord_eachrow[0]) - 1
                atomic_num = int(coord_eachrow[1])
                atom_element = ATOMIC_NUM_TO_ELEMENT.get(atomic_num, 'Unknown')
                atom_x = float(coord_eachrow[-3])
                atom_y = float(coord_eachrow[-2])
                atom_z = float(coord_eachrow[-1])
                self.xyz_coord.append([atom_idx, atomic_num, atom_element, atom_x, atom_y, atom_z])
            else:
                break
            self.natom = len(self.xyz_coord)
            self.xyz_comment = comment
        return self.xyz_coord
    
    def get_xyzdata_stdorien(self, comment=''):
        """ Extract atomic coordinates from the log file.
            
        Attributes updated:
        ------------
            self.xyz_coord: list
                List of atomic coordinates [atom_idx, atomic_num, element, x, y, z].
                len(self.xyz_coord) = number of atoms.
        """
        for idx, line in enumerate(self.lines_reverse):
            if 'Standard orientation' in line:
                idx_coord_begin = self.lines_num - 1 - idx
                break
        else:
            raise Exception("Standard orientation section not found in the log file.")
        
        for idx, line in enumerate(self.lines_forward[idx_coord_begin + 5:]):
            if not '----' in line:
                coord_eachrow = line.split()
                atom_idx = int(coord_eachrow[0]) - 1
                atomic_num = int(coord_eachrow[1])
                atom_element = ATOMIC_NUM_TO_ELEMENT.get(atomic_num, 'Unknown')
                atom_x = float(coord_eachrow[3])
                atom_y = float(coord_eachrow[4])
                atom_z = float(coord_eachrow[5])
                self.xyz_coord.append([atom_idx, atomic_num, atom_element, atom_x, atom_y, atom_z])
            else:
                break
            self.natom = len(self.xyz_coord)
            self.xyz_comment = comment
        return self.xyz_coord
    
    def to_xyzstring(self, xyz_comment=''):
        """ Convert atomic coordinates to an XYZ string.
        
        Args:
        ------------
            xyz_comment (str, optional):
                A comment line for the XYZ file. Defaults to an empty string.
            
        Returns:
        ------------
            xyz_string (str):
                Atomic coordinates in XYZ format.
        """
        self.create_xyz_df()
        xyz_num = len(self.xyz_df)
        xyz_string = f"{xyz_num}\n"
        xyz_string += f"{xyz_comment}\n"
        for _, row in self.xyz_df.iterrows():
            xyz_string += f"{row['element']} {row['x']} {row['y']} {row['z']}\n"
        return xyz_string
    
    def to_xyzstring_fromanydf(self, xyz_df, xyz_comment=''):
        """ Convert atomic coordinates to an XYZ string.
        
        Args:
        ------------
            xyz_comment (str, optional):
                A comment line for the XYZ file. Defaults to an empty string.
            
        Returns:
        ------------
            xyz_string (str):
                Atomic coordinates in XYZ format.
        """
        xyz_num = len(xyz_df)
        xyz_string = f"{xyz_num}\n"
        xyz_string += f"{xyz_comment}\n"
        for _, row in xyz_df.iterrows():
            xyz_string += f"{row['element']} {row['x']} {row['y']} {row['z']}\n"
        return xyz_string

    def to_xyzfile(self, filename, xyz_comment=''):
        """ Write atomic coordinates to an XYZ file.
        
        Args:
        ------------
            filename (str):
                The name of the XYZ file to be written.
            xyz_comment (str, optional):
                A comment line for the XYZ file. Defaults to an empty string.
        """
        xyz_string = self.to_xyzstring(xyz_comment)
        with open(filename, 'w') as f:
            f.write(xyz_string)

    def create_xyz_df(self):
        """ Create DataFrame for atomic information
        
        Attributes updated:
        ------------
            self.xyz_df: DataFrame
                DataFrame of atomic information. Columns are ['atom_idx', 'atomic_num', 'element', 'x', 'y', 'z', 'cov_radius']
                atom_idx: index of atom
                atomic_num: atomic number of atom
                element: element symbol of atom
                x, y, z: atomic coordinates
                cov_radius: covalent radius of atom
        
        """
        self.get_xyzdata()
        self.xyz_df = pd.DataFrame(self.xyz_coord, columns=['atom_idx', 'atomic_num', 'element', 'x', 'y', 'z'])
        self.xyz_df['cov_radius'] = self.xyz_df['element'].apply(lambda x: COVALENT_RADII.get(x, 0.0))
        return self.xyz_df
    
    def create_distance_matrix(self):
        """ Create distance matrix from atomic coordinates.
        
        Attributes updated:
        ------------
            self.distance_matrix: DataFrame
                DataFrame of distance matrix. Columns and rows are atom indices.
        """
        atom_coordinates = self.xyz_df[['x', 'y', 'z']].to_numpy()
        diff = atom_coordinates[:, np.newaxis, :] - atom_coordinates[np.newaxis, :, :]
        distance_matrix = np.sqrt(np.sum(diff**2, axis=-1))
        self.distance_matrix = pd.DataFrame(distance_matrix)
        return self.distance_matrix
    
    def get_mulliken_charges(self):
        """ Extract Mulliken charges from the log file.

        Returns:
        ------------
            data_mulliken (list of lists):
                A list of Mulliken charges for each atom in the format:
                [atom_idx, element, charge]
                - atom_idx (int): The index of the atom (0-based).
                - element (str): The element symbol of the atom.
                - charge (float): The Mulliken charge of the atom.
        """
        data_mulliken = []

        # Search for 'Mulliken charges:' in lines_reverse to get the starting index
        for idx, line in enumerate(self.lines_reverse):
            if 'Mulliken charges:' in line:
                idx_mulliken = self.lines_num - 1 - idx
                # Extract the relevant lines from lines_forward
                mulliken_lines = self.lines_forward[idx_mulliken + 2:idx_mulliken + 2 + len(self.xyz_coord)]
                parsed_data_mulliken = []

                for line in mulliken_lines:
                    try:
                        atom_idx, element, charge = line.split()
                        atom_idx = int(atom_idx) - 1  # Convert to 0-based index
                        charge = float(charge)
                        parsed_data_mulliken.append([atom_idx, element, charge])
                    except ValueError:
                        print(f"Error parsing line: {line}")

                return parsed_data_mulliken

        return data_mulliken  # Return an empty list if no 'Mulliken charges:' section is found


    def get_orbital_energies(self):
        """ Extract orbital energies from the log file

        Returns:
        ------------
            energies_occupied: list
                List of occupied orbital energies
            energies_virtual: list
                List of virtual orbital energies
            homo: float
                Energy level of HOMO
            lumo: float
                Energy level of LUMO
            Egap: float
                Energy of LUMO - HOMO
        """
        for idx, line in enumerate(self.lines_reverse):
            if 'The electronic state' in line:
                idx_MO = self.lines_num - 1 - idx
                break
        energies_occupied = []
        energies_virtual = []
        for idx, line in enumerate(self.lines_forward[idx_MO:]):
            if 'Alpha  occ. eigenvalues' in line:
                energy_orbital = line.split()[4:]
                energies_occupied += [float(val) for val in energy_orbital]
                continue
            if 'Alpha virt. eigenvalues' in line:
                energy_orbital = line.split()[4:]
                energies_virtual += [float(val) for val in energy_orbital]
                energies_virtual = energies_virtual[:10]
                continue
        homo = max(energies_occupied) if energies_occupied else None
        lumo = min(energies_virtual) if energies_virtual else None
        Egap = lumo - homo if homo is not None and lumo is not None else None
        return energies_occupied, energies_virtual, homo, lumo, Egap

    def get_bond_existed(self, cutoff=0.1):
        """ Calculates bond lengths between atoms based on their coordinates and covalent radii.
        A bond is considered to exist if the distance between two atoms is less than the sum of their covalent
        radii plus the cutoff value.

        Args:
        ------------
            cutoff (float, optional):
                A threshold value added to the sum of covalent radii to determine if atoms are bonded.
                Default is 0.1 Å.

        Returns:
        ------------
            bond_exist (list of tuples):
                A list containing bond information, where each tuple represents:
                (atom_idx_1, element_1, atom_idx_2, element_2, bond_length)
        """
        self.create_distance_matrix()
        
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
        """ Calculate the distance between two atoms

        Args:
        ------------
            atomA_idx: int
                index of atom A
            atomB_idx: int
                index of atom B

        Returns:
        ------------
            dist: np.float64
                distance between atom A and atom B
        """
        atomA_coord = self.xyz_df.loc[atomA_idx, ['x', 'y', 'z']].to_numpy()
        atomB_coord = self.xyz_df.loc[atomB_idx, ['x', 'y', 'z']].to_numpy()
        dist = np.linalg.norm(atomA_coord - atomB_coord)
        return dist
    
    def calculate_angle(self, atomA_idx, atomB_idx, atomC_idx):
        """ Calculate the angle between three atoms (A-B-C).
        atomB is the central atom.

        Args:
        ------------
            atomA_idx (int): 
                Index of atom A.
            atomB_idx (int): 
                Index of atom B.
            atomC_idx (int): 
                Index of atom C.

        Returns:
        ------------
            angle_degrees (float):
                Angle between three atoms (A-B-C) in degrees.
        """
        atomA_coord = self.xyz_df.loc[atomA_idx, ['x', 'y', 'z']].to_numpy()
        atomB_coord = self.xyz_df.loc[atomB_idx, ['x', 'y', 'z']].to_numpy()
        atomC_coord = self.xyz_df.loc[atomC_idx, ['x', 'y', 'z']].to_numpy()
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

        
    def get_SCF(self):
        lines = self.lines_reverse
        job_finished = False
        found_energy = False
        try:
            for idx, line in enumerate(lines):
                if 'Normal termination' in line:
                    job_finished = True
                    continue
                if job_finished and 'SCF Done' in line:
                    found_energy = True
                    energySCF = float(line.split()[4])
                    break
            return energySCF
        except:
            return np.nan
    
    # def get_energy_and_geom(self):
    #     lines = self.lines_forward
    #     list_energy_geom = []
    #     for line in lines:
    #         if 'Standard orientation' in line:
                
    #         if 'SCF Done' in line:
    #             energySCF = float(line.split()[4])
    
class ReadSP(MainReader):
    def __init__(self, logfile):
        super().__init__(logfile)  # Inherit and initialize MainReader with logfile
        self.jobtype()
        if self.job[0] != 'energy':
            raise EnvironmentError("This is not the output of a single-point energy calculation")
        if self.job_finished[0] is False:
            raise RuntimeError("Job did not terminate normally (no 'Normal termination').")

        self.raw_coord = []
        self.get_final_coord()
        
    def get_final_coord(self, xyz_comment='', xyzfiletowrite=None, writexyzfile=False):
        # coordinates from the end of job
        coordinates = []
        reading = False
        lines = self.lines_forward
        for line in lines:
            if 'Unable to Open any file for archive entry.' in line:
                reading = True
                continue
            elif 'The archive entry for this job was punched.' in line:
                reading = False
                continue
            if reading:
                coordinates.append(line.strip())  # Collect lines between the two markers

        coordinates = ' '.join(coordinates)
        coordinates = coordinates.replace(' ', '')
        coordinates_parts = coordinates.split(r'\\')
        if len(coordinates_parts) > 3:
            coordinates = coordinates_parts[3].replace("\\", "\n")
        else:
            raise ValueError("Invalid coordinate data in the logfile.")
        coordinates = coordinates[4:]
        coordinates = coordinates.split('\n')
        xyz_coord = []
        for eachline in coordinates:
            coord = eachline.split(',')
            if len(coord) >= 4:  # Ensure the line contains sufficient data
                coord_element = coord[0]
                coord_x = coord[-3]
                coord_y = coord[-2]
                coord_z = coord[-1]
                self.raw_coord.append([coord_element, coord_x, coord_y, coord_z])
                xyz_coord.append(f'{coord_element}\t{coord_x}\t{coord_y}\t{coord_z}')
        
        # compare the number of atoms from the xyz_coord with the number of atoms expected from the stoichiometry (self.natoms)
        xyz_num = len(xyz_coord)
        if xyz_num == self.natoms:
            pass
        else:
            raise ValueError(f"Mismatch in number of atoms: Expected {self.natoms}, but got {xyz_num} from XYZ structure.")

        # Write to XYZ file if requested
        if writexyzfile:
            if xyzfiletowrite is None:
                raise ValueError("The name of the XYZ file must be provided if writexyzfile=True")
            with open(xyzfiletowrite, 'w') as f:
                f.write(f'{xyz_num}\n')
                f.write(f'{xyz_comment}\n')
                for coord in xyz_coord:
                    f.write(f'{coord}\n')
        # return xyz_num, xyz_comment, xyz_coord
        
    def create_xyz_df(self):
        df = pd.DataFrame(self.raw_coord, columns=['element', 'x', 'y', 'z'])
        df[['x', 'y', 'z']] = df[['x', 'y', 'z']].astype(float)
        df['atomic_num'] = df['element'].apply(lambda x: ELEMENT_TO_ATOMIC_NUM.get(x, 'Unknown'))
        df['cov_radius'] = df['element'].apply(lambda x: COVALENT_RADII.get(x, 0.0))
        df['atom_idx'] = df.index # 0-based index
        df = df[['atom_idx', 'atomic_num', 'element', 'x', 'y', 'z', 'cov_radius']]
        self.xyz_df = df
        return self.xyz_df
    
    def writexyzstr(self, xyz_comment=''):
        xyz_num = self.natoms
        self.create_xyz_df()
        xyz_df = self.xyz_df
        xyz_string = f"{xyz_num}\n"
        xyz_string += f"{xyz_comment}\n"
        for _, row in xyz_df.iterrows():
            xyz_string += f"{row['element']} {row['x']} {row['y']} {row['z']}\n"
        return xyz_string
    
    def writexyzfile(self, filename, xyz_comment=''):
        xyz_string = self.writexyzstr(xyz_comment)
        with open(filename, 'w') as f:
            f.write(xyz_string)
    
    def energySCF(self):
        lines = self.lines_reverse
        for line in lines:
            if 'SCF Done' in line:
                energySCF = float(line.split()[4])
                break
        return energySCF

# class ReadOpt:
    
# class ReadFreq:
    
class ReadIRC(MainReader):
    def __init__(self, logfile):
        super().__init__(logfile)
        self.jobtype()
        if self.job[0] != 'irc':
            raise EnvironmentError("This is not the output of an IRC calculation")
        if self.job_finished[0] is False:
            raise RuntimeError("Job did not terminate normally (no 'Normal termination').")

    def get_IRC_Energy(self):
        lines = self.lines_reverse
        job_finished = False
        found_IRC = False
        try:
            for idx, line in enumerate(lines):
                if 'Normal termination' in line:
                    job_finished = True
                    continue
                if job_finished and 'IRC' in line:
                    found_IRC = True
                    idx_end = idx + 6
                    continue
                if found_IRC and 'Summary of reaction path following' in line:
                    idx_begin = idx + 2
                    break
            energyTS = float(lines[idx_begin].split()[-1])
            if 'FORWARD' in lines[idx_begin + 3]:
                energyCritical = float(lines[idx_end].split()[-2])
            if 'REVERSE' in lines[idx_begin + 3]:
                energyCritical = float(lines[idx_begin - 5].split()[-2])
            else:
                pass
        except:
            return np.nan, np.nan, np.nan
        return {'E_TS': energyTS, 
                'dE_tominima': energyCritical, 
                'E_minima': energyTS + energyCritical}
    
    def get_IRC_struct_last(self):
        lines = self.lines_reverse
        for idx, line in enumerate(lines):
            if 'Input orientation' in line:
                break
        idx_start = self.lines_num - 1 - idx
        table = self.lines_forward[idx_start + 5 : idx_start + 5 + self.natoms]
        return self._read_single_orientation_table(table)
                # self._read_single_orientation_table()
#         if not self.job_finished:
#             raise Exception('Job not finished')
#         lines = self.lines_forward
#         idx_IRC_Structs = []
#         all_IRC_Structs = []
        
#         for idx, line in enumerate(lines):
#             if 'CURRENT STRUCTURE' in line:
#                 idx_start = idx + 6
#                 idx_IRC_Structs.append(idx_start)
#         print(idx_IRC_Structs)
#         print(self.natom)
        
#         for idx_IRC_start in idx_IRC_Structs:
#             print(idx_IRC_start)
#             print(lines[idx_IRC_start : idx_IRC_start + self.natom])
#             # for idx, line in enumerate(lines[idx_IRC_start : idx_IRC_start + self.natom]):
#             #     print(line)
#         #         xyz_coord = []
#         #         coord_eachrow = line.split()
#         #         print(coord_eachrow[0])
#         #         atom_idx = int(coord_eachrow[0]) - 1
#         #         atomic_num = int(coord_eachrow[1])
#         #         atom_symbol = ATOMIC_NUM_TO_ELEMENT.get(atomic_num, 'Unknown')
#         #         atom_x = float(coord_eachrow[-3])
#         #         atom_y = float(coord_eachrow[-2])
#         #         atom_z = float(coord_eachrow[-1])
#         #         xyz_coord.append([atom_idx, atomic_num, atom_symbol, atom_x, atom_y, atom_z])
                
#         #         xyz_df = pd.DataFrame(self.xyz_coord, columns=['atom_idx', 'atomic_num', 'element', 'x', 'y', 'z'])
#         #         xyz_df['cov_radius'] = self.xyz_df['element'].apply(lambda x: COVALENT_RADII.get(x, 0.0))
#         #         xyz_string = self.to_xyzstring_fromanydf(xyz_df)
#         #     all_IRC_Structs.append(xyz_string)
#         # return all_IRC_Structs