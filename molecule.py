#!/usr/bin/env python3
# molecule.py
#
# A base Molecule class for reading atoms from an XYZ file, building
# a covalent adjacency matrix using data from elements_table.py,
# and providing geometry utilities (distance, angle, dihedral).
#
# This example includes:
#  1) An inline Atom dataclass.
#  2) The Molecule class with the methods described.
#  3) A command-line "main" entry point for testing.
# 
# Wed Feb 19 2025
# - added Gaussian-style standard orientaion method
# - universal XYZ `read_xyz_data()` method for reading any type of XYZ 

import sys
import math
import numpy as np
from dataclasses import dataclass
import re
from typing import List, Tuple, Union, Optional

import pprint 

# Make sure elements_table.py is in the same directory or installable path
from elements_table import Elements, _DISTANCE_UNIT_ALIASES


@dataclass
class Atom:
    """
    A lightweight container for atomic data.
    """
    symbol: str
    x: float
    y: float
    z: float
    charge: float = 0.0


class Molecule:
    """
    A base class representing a molecular system with covalent bonding
    and fragment detection. Coordinates are in Å by convention.

    Key attributes:
        title (str): A descriptive name for the molecule.
        atoms (List[Atom]): List of Atom objects.
        bond_matrix (np.ndarray): 2D adjacency matrix of shape (n, n), storing bond information.
        fragments (Dict[int, List[int]]): Mapping of fragment_id -> list of atom indices.
    """

    def __init__(self, title: str = ""):
        self.title = title
        self.atoms = []               # type: List[Atom]
        self.bond_matrix = None       # type: np.ndarray
        self.fragments = {}           # type: Dict[int, List[int]]


    @classmethod
    def from_atoms(cls, atoms: List[Tuple[str, float, float, float]], title: str = "") -> "Molecule":
        """
        Creates a Molecule instance from a list of atom tuples.

        Each tuple should be of the form: (symbol, x, y, z)
        """
        mol = cls(title=title)
        for atom in atoms:
            symbol, x, y, z = atom
            mol.atoms.append(Atom(symbol, x, y, z))
        n = len(mol.atoms)
        mol.bond_matrix = np.zeros((n, n), dtype=float)
        return mol

    def read_xyz_data(self,
                      data: Optional[Union[str, List[str]]] = None,
                      file_name: Optional[str] = None,
                      units: str = "Ang") -> None:
        """
        Reads XYZ data from a file or directly from a string/list of lines, converts
        the coordinates to Ångströms (if necessary), and populates the molecule's atoms and title.
    
        Parameters:
            data: Either a multi-line string or a list of strings representing the XYZ data.
            file_name: Optional path to a file containing XYZ data. If provided, this takes precedence over 'data'.
            units: Units in which the coordinates are given. Defaults to "Ang" (Ångströms).
                   Other supported units include "pm", "nm", "bohr", etc.
        """
        # If a file name is provided, read from file.
        if file_name is not None:
            with open(file_name, 'r') as f:
                data = f.read()
    
        if data is None:
            raise ValueError("No XYZ data provided (neither data nor file_name).")
    
        # Normalize data to a list of lines.
        if isinstance(data, str):
            lines = data.splitlines()
        else:
            lines = data
    
        # Remove empty or whitespace-only lines.
        lines = [line.strip() for line in lines if line.strip()]
        if not lines:
            raise ValueError("No valid lines found in the provided data.")
    
        # Determine conversion factor to convert coordinates to Å.
        # The conversion factor is computed as 1.0 / factor from the unit alias dictionary from elements_table 
        unit_key = units.lower()
        if unit_key in _DISTANCE_UNIT_ALIASES:
            _, factor = _DISTANCE_UNIT_ALIASES[unit_key]
            conversion_factor = 1.0 / factor
        else:
            conversion_factor = 1.0
    
        title = ""
        atom_lines: List[str] = []
    
        # Check if the first line is numeric (possible header).
        try:
            expected_atom_count = int(lines[0])
            if len(lines) == expected_atom_count + 1:
                # No comment line.
                atom_lines = lines[1:]
            elif len(lines) >= expected_atom_count + 2:
                # If the second line is non-numeric, treat it as a comment/title.
                try:
                    int(lines[1])
                    # If it parses as an int, then no comment is present.
                    atom_lines = lines[1:expected_atom_count+1]
                except ValueError:
                    title = lines[1]
                    atom_lines = lines[2:expected_atom_count+2]
            else:
                # Fallback to header-less parsing.
                atom_lines = lines
        except ValueError:
            # First line isn't numeric; assume header-less format.
            atom_lines = lines
    
        # Validate that we have the expected number of atom lines, if applicable.
        if 'expected_atom_count' in locals() and len(atom_lines) != expected_atom_count:
            print(f"Warning: Expected {expected_atom_count} atoms, but found {len(atom_lines)} lines.")
    
        # Populate the molecule's attributes.
        self.title = title
        self.atoms = []
        for line in atom_lines:
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Line does not have enough entries for an atom: '{line}'")
            identifier = parts[0]
            if identifier.isdigit():
                symbol = Elements.symbol(int(identifier))
            else:
                symbol = identifier
            try:
                x, y, z = map(float, parts[1:4])
                # Convert the coordinates to Å.
                x *= conversion_factor
                y *= conversion_factor
                z *= conversion_factor
            except ValueError as e:
                raise ValueError(f"Error parsing coordinates in line: '{line}'. {e}")
            self.atoms.append(Atom(symbol, x, y, z))
    
        n = len(self.atoms)
        self.bond_matrix = np.zeros((n, n), dtype=float)

    @classmethod
    def from_xyz_data(cls,
                      data: Optional[Union[str, List[str]]] = None,
                      file_name: Optional[str] = None,
                      units: str = "Ang") -> "Molecule":
        """
        Convenience class method to construct a Molecule from XYZ data.

        This method creates a new Molecule, calls read_xyz_data() to populate it,
        and returns the Molecule.

        Parameters:
            data: Either a multi-line string or a list of strings representing the XYZ data.
            file_name: Optional path to a file containing XYZ data. If provided, this takes precedence over 'data'.
            units: Units in which the coordinates are given. Defaults to "Ang" (Ångströms).

        Returns:
            A new Molecule instance populated with the XYZ data.
        """
        mol = cls()
        mol.read_xyz_data(data=data, file_name=file_name, units=units)
        return mol


    def detect_bonds(self, tolerance: float = 0.3):
        """
        Fills self.bond_matrix by comparing interatomic distances
        to the sum of covalent radii (plus a tolerance).

        :param tolerance: Additional margin (in Å) on top of covalent radius sums.
                          Adjust as needed for borderline cases.
        """
        n = len(self.atoms)
        for i in range(n):
            for j in range(i + 1, n):
                dist_ij = self.distance(i, j)
                r1 = Elements.covalent_radius(self.atoms[i].symbol, order="single")
                r2 = Elements.covalent_radius(self.atoms[j].symbol, order="single")
                r_sum = r1 + r2 + tolerance

                if dist_ij <= r_sum:
                    # Simple approach: call it a single bond
                    self.bond_matrix[i, j] = 1
                    self.bond_matrix[j, i] = 1

                    # More advanced logic for double/triple bonds could go here

    def find_fragments(self):
        """
        Identify connected components ("fragments") in the bond_matrix.
        A DFS approach populates self.fragments = {frag_id: [atom_indices]}.
        """
        n = len(self.atoms)
        visited = set()
        frag_id = 0
        self.fragments = {}

        for start_atom in range(n):
            if start_atom not in visited:
                stack = [start_atom]
                connected = []
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        connected.append(current)

                        # look for neighbors
                        for neigh in range(n):
                            if self.bond_matrix[current, neigh] > 0.0 and neigh not in visited:
                                stack.append(neigh)

                # store the connected fragment
                self.fragments[frag_id] = connected
                frag_id += 1

    def distance(self, i: int, j: int) -> float:
        r""" 
        Returns the Euclidean distance between atoms i and j (in ang).

        $$
        d_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}
        $$
        """
        ax, ay, az = self.atoms[i].x, self.atoms[i].y, self.atoms[i].z
        bx, by, bz = self.atoms[j].x, self.atoms[j].y, self.atoms[j].z
        return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)

    def angle(self, i: int, j: int, k: int, degrees: bool = True) -> float:
        """
        Returns the angle at atom j formed by (i - j - k).
        If degrees=False, return radians.

        $$
        \\theta = \\cos^{-1} \\Bigl(
            \\frac{(\\mathbf{r}_i - \\mathbf{r}_j) \\cdot (\\mathbf{r}_k - \\mathbf{r}_j)}
                 {\\|\\mathbf{r}_i - \\mathbf{r}_j\\| \\, \\|\\mathbf{r}_k - \\mathbf{r}_j\\|}
        \\Bigr)
        $$
        """
        r_i = np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])
        r_j = np.array([self.atoms[j].x, self.atoms[j].y, self.atoms[j].z])
        r_k = np.array([self.atoms[k].x, self.atoms[k].y, self.atoms[k].z])

        v_ji = r_i - r_j
        v_jk = r_k - r_j

        dot_val = np.dot(v_ji, v_jk)
        mag_ji = np.linalg.norm(v_ji)
        mag_jk = np.linalg.norm(v_jk)

        cos_theta = dot_val / (mag_ji * mag_jk)

        # numerical safety
        cos_theta = max(min(cos_theta, 1.0), -1.0)

        theta_radians = np.arccos(cos_theta)
        return np.degrees(theta_radians) if degrees else theta_radians

    def dihedral(self, i: int, j: int, k: int, l: int, degrees: bool = True) -> float:
        """
        Returns the dihedral angle formed by atoms (i - j - k - l).
        If degrees=False, returns radians.

        A standard approach uses cross products to find normal vectors
        and the sign of the torsion from the 'atan2' of those vectors.
        """
        r_i = np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])
        r_j = np.array([self.atoms[j].x, self.atoms[j].y, self.atoms[j].z])
        r_k = np.array([self.atoms[k].x, self.atoms[k].y, self.atoms[k].z])
        r_l = np.array([self.atoms[l].x, self.atoms[l].y, self.atoms[l].z])

        b1 = r_i - r_j
        b2 = r_k - r_j
        b3 = r_l - r_k

        # normal vectors
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        # handle degenerate vectors (avoid /0)
        if np.linalg.norm(n1) == 0 or np.linalg.norm(n2) == 0:
            return 0.0  # or raise an exception

        # normalize
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        # cross product to get sign
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        angle_radians = np.arctan2(y, x)

        return np.degrees(angle_radians) if degrees else angle_radians


    def to_standard_orientation(self, use_atomic_numbers: bool = True) -> None:
        """
        Converts the molecule's atomic coordinates to a standard orientation.
        
        The procedure is as follows:
        
        1. **Centering:**  
           Compute the center of nuclear charge (or mass) and translate all coordinates so that
           the center is at the origin. That is,
           $$
           \\mathbf{T} = \\frac{\\sum_i w_i \\mathbf{r}_i}{\\sum_i w_i}, \\quad
           \\mathbf{r}_i' = \\mathbf{r}_i - \\mathbf{T},
           $$
           where $w_i$ is either the atomic number or atomic mass.
        
        2. **Inertia Tensor Construction:**  
           Build the inertia tensor using the weighted centered coordinates:
           $$
           I = \\sum_i w_i \\Bigl(\\|\\mathbf{r}_i'\\|^2\\, \\mathbf{I}_{3\\times3} - \\mathbf{r}_i' \\mathbf{r}_i'^{T}\\Bigr).
           $$
        
        3. **Diagonalization:**  
           Diagonalize $I$ to obtain eigenvalues and eigenvectors (the principal axes).
        
        4. **Rotation Matrix Construction:**  
           Adopt the following convention:
           - New $x$-axis is the eigenvector corresponding to the largest eigenvalue,
           - New $y$-axis is the eigenvector corresponding to the intermediate eigenvalue,
           - New $z$-axis is the eigenvector corresponding to the smallest eigenvalue.
           
           That is, if $\lambda_0 \le \\lambda_1 \\le \\lambda_2$, then set
           $$
           R = \\begin{pmatrix} \\mathbf{v}_{x} & \\mathbf{v}_{y} & \\mathbf{v}_{z} \\end{pmatrix}
           = \\begin{pmatrix} \\mathbf{v}_{\\lambda_2} & \\mathbf{v}_{\\lambda_1} & \\mathbf{v}_{\\lambda_0} \\end{pmatrix}.
           $$
           Finally, ensure that $\\det(R) = +1$ (right-handed coordinate system).
        
        5. **Coordinate Transformation:**  
           The final standard orientation coordinates are given by:
           $$
           \\mathbf{r}_{\\text{std}} = (\\mathbf{r} - \\mathbf{T}) \\cdot R.
           $$
           Update the molecule's atomic coordinates with these values.
        """

        # Step 1: Extract coordinates and compute the reference center.
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])
        if use_atomic_numbers:
            weights = np.array([Elements.atomic_number(atom.symbol) for atom in self.atoms])
        else:
            weights = np.array([Elements.mass(atom.symbol) for atom in self.atoms])
        total_weight = np.sum(weights)
        T = np.sum(weights[:, None] * coords, axis=0) / total_weight
        coords_centered = coords - T

        # Step 2: Construct the inertia tensor using centered coordinates.
        I = np.zeros((3, 3))
        for i, (x, y, z) in enumerate(coords_centered):
            w = weights[i]
            I[0, 0] += w * (y**2 + z**2)
            I[1, 1] += w * (x**2 + z**2)
            I[2, 2] += w * (x**2 + y**2)
            I[0, 1] -= w * x * y
            I[0, 2] -= w * x * z
            I[1, 2] -= w * y * z
        # Enforce symmetry:
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]

        # Step 3: Diagonalize the inertia tensor.
        # np.linalg.eigh returns eigenvalues in ascending order.
        eigenvalues, eigenvectors = np.linalg.eigh(I)

        # Step 4: Reorder eigenvectors to define the new axes.
        # Convention: new_x = eigenvector with largest eigenvalue,
        #             new_y = eigenvector with intermediate eigenvalue,
        #             new_z = eigenvector with smallest eigenvalue.
        idx = np.argsort(eigenvalues)  # idx[0]: smallest, idx[2]: largest.
        v_small = eigenvectors[:, idx[0]]
        v_mid   = eigenvectors[:, idx[1]]
        v_large = eigenvectors[:, idx[2]]
        # Build the rotation matrix (columns are new axes: x, y, z).
        R = np.column_stack((v_large, v_mid, v_small))

        # Step 5: Ensure the rotation matrix is right-handed.
        if np.linalg.det(R) < 0:
            R[:, 0] *= -1  # Flip the new x-axis if needed.

        # Step 6: Rotate the centered coordinates.
        coords_std = np.dot(coords_centered, R)

        # Update the atom coordinates,
        # swapping x and z axes to match Gaussian
        for i, atom in enumerate(self.atoms):
            atom.z, atom.y, atom.x = coords_std[i]


    def summary(self) -> str:
        """
        Returns a human-readable summary of the molecule:
        number of atoms, fragments, etc.
        """
        lines = []
        lines.append(f"Title: {self.title}")
        lines.append(f"Number of atoms: {len(self.atoms)}")
        if self.bond_matrix is not None:
            lines.append("Bond matrix is initialized.")
        if self.fragments:
            lines.append(f"Number of fragments: {len(self.fragments)}")
        else:
            lines.append("Fragments not yet determined.")
        return "\n".join(lines)


    def to_xyz(self) -> str:
        """
        Returns a well-formed XYZ format for the molecule as a string.

        The output format is as follows:

        <number of atoms>
        <title>
        <atom symbol> <x> <y> <z>

        Each atom line is formatted with the atomic symbol left-aligned in a field of 5 characters,
        and each coordinate is printed in a field of 17 characters with 12 decimal places.
        """
        lines = []
        # First line: the number of atoms.
        lines.append(f"{len(self.atoms)}")
        # Second line: use the molecule's title if provided
        comment = self.title 
        lines.append(comment)
        # Format each atom line.
        for atom in self.atoms:
            lines.append(f"{atom.symbol:<5}{atom.x:17.12f}{atom.y:17.12f}{atom.z:17.12f}")
        return "\n".join(lines)



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python molecule.py <xyz_file>")
        sys.exit(1)
    
    
    xyz_file = sys.argv[1]
    # Example of direct parsing using the instance method
    mol = Molecule().from_xyz_data(file_name=xyz_file)

    print("\nInput Orientation:")
    print(mol.to_xyz())

    mol.to_standard_orientation()
    print("\nStandard Orientation:")
    print(mol.to_xyz())

    mol.detect_bonds(tolerance=0.3)
    mol.find_fragments()
    print("\n" + mol.summary())

    # --- print bond matrix ---
    print("\nBond matrix:")
    n = len(mol.atoms)
    header = "    " + " ".join([f"{j:>2}" for j in range(n)])
    print(header)
    for i in range(n):
        row_str = " ".join(
            f"{int(mol.bond_matrix[i, j]):>2}" if mol.bond_matrix[i, j] != 0 else f"{'•':>2}"
            for j in range(n)
        )
        print(f"{i:>2}  {row_str}")
    
    # Example geometry checks
    if n >= 2:
        d_01 = mol.distance(0, 1)
        print(f"\nDistance between atom 0 and 1: {d_01:.3f} Å")

    if n >= 4:
        torsion_angle = mol.dihedral(0, 1, 2, 3)
        print(f"Dihedral angle (0-1-2-3): {torsion_angle:.2f} degrees")

