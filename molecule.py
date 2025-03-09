#!/usr/bin/env python3
# molecule.py
#
# A base Molecule class for reading atoms from an XYZ file, building
# a covalent adjacency matrix using data from elements_table.py,
# and providing geometry utilities (distance, angle, dihedral).
#
#

import sys
import math
import numpy as np
from dataclasses import dataclass
import re
from typing import List, Tuple, Union, Optional, Dict

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

    Additional metadata attributes:
        XYZ_Comment (str | None): The raw comment line from the XYZ file (if any).
        Energy (float | None): Parsed or user-provided energy (units unspecified).
        file_name (str | None): The file name (minus .xyz) from which data was read.
        frame_number (int | None): The index or ID for this molecule in a trajectory.

        metadata (Dict[str, float]): A dictionary for additional labeled numeric metadata,
                                     e.g. "coord", "time", etc.
    """

    def __init__(self, title: str = ""):
        self.title = title
        self.atoms = []               # type: List[Atom]
        self.bond_matrix = None       # type: np.ndarray
        self.fragments = {}           # type: Dict[int, List[int]]

        # New metadata fields
        self.XYZ_Comment = None       # Raw second line (or user-supplied comment)
        self.Energy = None            # Parsed or user-supplied energy
        self.file_name = None         # File name minus .xyz
        self.frame_number = None      # Frame index in a trajectory

        # Dictionary to hold arbitrary numeric metadata like coord, time, etc.
        self.metadata = {}            # type: Dict[str, float]


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
                      units: str = "Ang",
                      xyz_comment: Optional[str] = None,
                      energy: Optional[float] = None,
                      frame_number: Optional[int] = None,
                      metadata: Optional[Dict[str, float]] = None) -> None:
        """
        Reads XYZ data from a file or directly from a string/list of lines, converts
        the coordinates to Ångströms (if necessary), and populates the molecule's
        atoms, title, and metadata fields.

        Parameters:
            data:
                Either a multi-line string or a list of strings representing the XYZ data.
            file_name:
                Path to a file containing XYZ data. If provided, this takes precedence over 'data'.
                We also store the base file name (minus .xyz) in self.file_name.
            units:
                Units in which the coordinates are given. Defaults to "Ang" (Ångströms).
                Other supported units include "pm", "nm", "bohr", etc.
            xyz_comment:
                Optional explicit comment for the molecule, overriding anything read
                from the second line of the XYZ file.
            energy:
                Optional explicit energy, overriding anything parsed from the comment.
            frame_number:
                If provided, identifies which frame of a trajectory this data belongs to.
            metadata:
                Optional dictionary of additional numeric metadata to store or override,
                e.g. {"coord": 25.38, "time": 10.5}.
        """
        # If a file name is provided, read from file.
        if file_name is not None:
            with open(file_name, 'r') as f:
                data = f.read()
            # Strip off .xyz if present, store as self.file_name
            base = file_name
            if base.endswith(".xyz"):
                base = base[:-4]
            self.file_name = base

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
        unit_key = units.lower()
        if unit_key in _DISTANCE_UNIT_ALIASES:
            _, factor = _DISTANCE_UNIT_ALIASES[unit_key]
            conversion_factor = 1.0 / factor
        else:
            conversion_factor = 1.0

        # We'll store an intermediate variable that might become self.XYZ_Comment if none is passed.
        potential_comment_line = ""
        atom_lines: List[str] = []

        # Check if the first line is numeric (possible header).
        try:
            expected_atom_count = int(lines[0])
            # We have a numeric first line => standard XYZ format with atom count
            if len(lines) == expected_atom_count + 1:
                # There's no separate comment line, everything else is presumably atom lines
                atom_lines = lines[1:]
            elif len(lines) >= expected_atom_count + 2:
                # If the second line is non-numeric, treat it as a comment/title.
                try:
                    int(lines[1])
                    # If it parses as an int, then no comment line is present.
                    atom_lines = lines[1:expected_atom_count + 1]
                except ValueError:
                    # The second line is indeed some sort of comment
                    potential_comment_line = lines[1]
                    atom_lines = lines[2:expected_atom_count + 2]
            else:
                # Fallback to header-less parsing if something is mismatched
                atom_lines = lines

        except ValueError:
            # First line isn't numeric; assume header-less format => all lines are atom lines
            atom_lines = lines

        # If the user provided an explicit xyz_comment, use that, else use the potential_comment_line
        if xyz_comment is not None:
            self.XYZ_Comment = xyz_comment
        else:
            # If there's a second line from the file, store it here
            self.XYZ_Comment = potential_comment_line if potential_comment_line else None

        # Parse or store Energy
        if energy is not None:
            self.Energy = energy
        else:
            # Attempt to parse from the comment line if present
            if self.XYZ_Comment:
                self.Energy = self._parse_energy(self.XYZ_Comment)
            else:
                self.Energy = None

        # Set frame_number if provided
        self.frame_number = frame_number

        # We'll decide on the final self.title
        if not self.title:
            if not self.XYZ_Comment:
                # If no comment line at all, try file_name
                if self.file_name and self.frame_number is not None:
                    self.title = f"{self.file_name}_{self.frame_number:03d}"
                elif self.file_name:
                    self.title = self.file_name
                else:
                    self.title = ""
            else:
                # We do have some comment line, but it might or might not be meaningful.
                # We'll keep it simple: use the entire comment if not purely numeric.
                if self._looks_like_just_energy(self.XYZ_Comment) and self.file_name:
                    if self.frame_number is not None:
                        self.title = f"{self.file_name}_{self.frame_number:03d}"
                    else:
                        self.title = self.file_name
                else:
                    self.title = self.XYZ_Comment

        # Now parse the atom lines
        if 'expected_atom_count' in locals() and len(atom_lines) != expected_atom_count:
            print(f"Warning: Expected {expected_atom_count} atoms, but found {len(atom_lines)} lines.")

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

        # --- NEW: parse extra labeled metadata (excluding energy labels) ---
        discovered_metadata = {}
        if self.XYZ_Comment:
            discovered_metadata = self._parse_metadata(self.XYZ_Comment)

        # If user also provided metadata explicitly, merge/override
        if metadata:
            discovered_metadata.update(metadata)

        # Store everything in self.metadata
        self.metadata.update(discovered_metadata)


    @classmethod
    def from_xyz_data(cls,
                      data: Optional[Union[str, List[str]]] = None,
                      file_name: Optional[str] = None,
                      units: str = "Ang",
                      xyz_comment: Optional[str] = None,
                      energy: Optional[float] = None,
                      frame_number: Optional[int] = None,
                      metadata: Optional[Dict[str, float]] = None) -> "Molecule":
        """
        Convenience class method to construct a Molecule from XYZ data.

        This method creates a new Molecule, calls read_xyz_data() to populate it,
        and returns the Molecule.

        Parameters:
            data: Either a multi-line string or a list of strings representing the XYZ data.
            file_name: Optional path to a file containing XYZ data. If provided, this takes precedence over 'data'.
            units: Units in which the coordinates are given. Defaults to "Ang" (Ångströms).
            xyz_comment: Optional explicit comment line.
            energy: Optional explicit energy value.
            frame_number: Optional frame index.
            metadata: Optional dictionary of labeled numeric metadata (coord, time, etc.).
        """
        mol = cls()
        mol.read_xyz_data(data=data,
                          file_name=file_name,
                          units=units,
                          xyz_comment=xyz_comment,
                          energy=energy,
                          frame_number=frame_number,
                          metadata=metadata)
        return mol

    def compute_distance_matrix(self) -> np.ndarray:
        """
        Computes the pairwise Euclidean distance matrix between atoms using vectorized operations.

        If the coordinate array is defined as $C$ with shape (n, 3), the distance between any two atoms
        is given by:

        $$
        d_{ij} = \sqrt{\sum_{k=1}^{3} (C_{ik} - C_{jk})^2}
        $$

        Returns:
            A numpy ndarray of shape (n, n) with distances.
        """
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])
        # Broadcasting to get differences: (n,1,3) - (1,n,3)
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
        return dist_matrix

    def detect_bonds(self, tolerance: float = 0.3):
        """
        Fills self.bond_matrix using a vectorized approach by comparing all interatomic distances
        to the sum of covalent radii (plus a tolerance).

        Parameters:
            tolerance: Additional margin (in Å) on top of covalent radius sums.
        """
        n = len(self.atoms)
        # Compute the full pairwise distance matrix
        dist_matrix = self.compute_distance_matrix()

        # Create an array of covalent radii for each atom (assuming single bond order)
        radii = np.array([Elements.covalent_radius(atom.symbol, order="single") for atom in self.atoms])
        # Build a threshold matrix: each element (i,j) is radii[i] + radii[j] + tolerance
        threshold_matrix = radii[:, None] + radii[None, :] + tolerance

        # A bond is present if the distance is non-zero and less than or equal to the threshold
        self.bond_matrix = ((dist_matrix <= threshold_matrix) & (dist_matrix > 0)).astype(float)

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
        Returns the Euclidean distance between atoms i and j (in Å).

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
        \theta = \cos^{-1} \Biggl(
            \frac{(\mathbf{r}_i - \mathbf{r}_j) \cdot (\mathbf{r}_k - \mathbf{r}_j)}
                 {\|\mathbf{r}_i - \mathbf{r}_j\| \, \|\mathbf{r}_k - \mathbf{r}_j\|}
        \Biggr)
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

        1. **Centering**
        2. **Inertia Tensor Construction**
        3. **Diagonalization**
        4. **Rotation Matrix Construction**
        5. **Coordinate Transformation**
        """
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

        # Compute weights based on atomic numbers or masses.
        if use_atomic_numbers:
            weights = np.array([Elements.atomic_number(atom.symbol) for atom in self.atoms])
        else:
            weights = np.array([Elements.mass(atom.symbol) for atom in self.atoms])

        # Step 1: Compute center and centered coordinates.
        coords_centered, T = self._compute_center(coords, weights)

        # Step 2: Construct the inertia tensor.
        I = self._compute_inertia_tensor(coords_centered, weights)

        # Step 3: Diagonalize the inertia tensor.
        eigenvalues, eigenvectors = self._diagonalize_inertia_tensor(I)

        # Step 4: Build the rotation matrix from eigenvectors.
        R = self._build_rotation_matrix(eigenvalues, eigenvectors)

        # Step 5: Rotate the centered coordinates.
        coords_std = self._apply_transformation(coords_centered, R)

        # Update the atom coordinates.
        self._update_atom_coordinates(coords_std)

    def _compute_center(self, coords: np.ndarray, weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the center of mass (or nuclear charge) and returns the centered coordinates.

        Parameters:
            coords: Array of atomic coordinates.
            weights: Array of weights (atomic numbers or masses) for each atom.

        Returns:
            A tuple (coords_centered, T) where T is the computed center.
        """
        total_weight = np.sum(weights)
        T = np.sum(weights[:, None] * coords, axis=0) / total_weight
        coords_centered = coords - T
        return coords_centered, T

    def _compute_inertia_tensor(self, coords_centered: np.ndarray, weights: np.ndarray) -> np.ndarray:
        """
        Constructs the inertia tensor from centered coordinates and weights.

        Parameters:
            coords_centered: Centered atomic coordinates.
            weights: Array of weights for each atom.

        Returns:
            A 3x3 inertia tensor.
        """
        I = np.zeros((3, 3))
        for i, (x, y, z) in enumerate(coords_centered):
            w = weights[i]
            I[0, 0] += w * (y**2 + z**2)
            I[1, 1] += w * (x**2 + z**2)
            I[2, 2] += w * (x**2 + y**2)
            I[0, 1] -= w * x * y
            I[0, 2] -= w * x * z
            I[1, 2] -= w * y * z
        # Enforce symmetry
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        return I

    def _diagonalize_inertia_tensor(self, I: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Diagonalizes the inertia tensor.

        Returns:
            (eigenvalues, eigenvectors), sorted in ascending order of eigenvalues.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(I)
        return eigenvalues, eigenvectors

    def _build_rotation_matrix(self, eigenvalues: np.ndarray, eigenvectors: np.ndarray) -> np.ndarray:
        """
        Builds the rotation matrix from the eigenvectors, adopting the
        convention that the new x-axis = largest eigenvalue, y-axis = middle, z-axis = smallest.
        Ensures a right-handed system.
        """
        # Sort indices: idx[0] smallest, idx[2] largest
        idx = np.argsort(eigenvalues)
        v_small = eigenvectors[:, idx[0]]
        v_mid   = eigenvectors[:, idx[1]]
        v_large = eigenvectors[:, idx[2]]

        # Build rotation matrix (columns: new x, y, z)
        R = np.column_stack((v_large, v_mid, v_small))
        # Ensure right-handed coordinate system
        if np.linalg.det(R) < 0:
            R[:, 0] *= -1
        return R

    def _apply_transformation(self, coords: np.ndarray, R: np.ndarray) -> np.ndarray:
        """
        Applies the rotation (or transformation) matrix to the coordinates.
        """
        return np.dot(coords, R)

    def _update_atom_coordinates(self, coords_std: np.ndarray) -> None:
        """
        Updates the molecule's atoms with the new coordinates.
        Swaps axes to match Gaussian's standard orientation convention.
        """
        for i, atom in enumerate(self.atoms):
            # Swap axes: assign atom.z, atom.y, atom.x from coords_std[i]
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

    @staticmethod
    def _parse_energy(comment_line: str) -> Optional[float]:
        """
        Attempt to parse a floating-point energy from the given comment line using a
        tiered approach:

        1) Look for a labeled pattern like "Energy=...", "E=...", ignoring case.
        2) If not found, see if there's exactly one float in the entire line.
        3) Otherwise, return None.
        """
        if not comment_line.strip():
            return None

        # Step 1: labeled pattern (case-insensitive)
        labeled_pattern = r"(?i)(?:energy|e)\s*[:=]?\s*([-+]?\d+(\.\d+)?([eE][+-]?\d+)?)"
        match = re.search(labeled_pattern, comment_line)
        if match:
            return float(match.group(1))

        # Step 2: if no labeled pattern, see if there's exactly one float
        all_floats = re.findall(
            r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
            comment_line
        )
        if len(all_floats) == 1:
            return float(all_floats[0])

        return None

    @staticmethod
    def _looks_like_just_energy(comment_line: str) -> bool:
        """
        Heuristic check: if the line yields one parseable float or a labeled 'energy',
        and there's minimal extra text, treat it as "just an energy."
        """
        if not comment_line.strip():
            return True

        # If we found an energy using the parse method...
        possible_energy = Molecule._parse_energy(comment_line)
        if possible_energy is None:
            return False

        # Check how many floats are in the line
        all_floats = re.findall(
            r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
            comment_line
        )

        # If there's exactly 1 float, and the rest is basically "energy/E" text,
        # treat it as just energy
        if len(all_floats) == 1:
            without_float = re.sub(
                r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
                "",
                comment_line
            )
            if len(without_float.strip()) < 15:  # arbitrary threshold
                return True

        return False

    # ---------------------------------------------------------------------
    # NEW METHOD: parse labeled numeric metadata (excluding energy labels).
    # ---------------------------------------------------------------------
    @staticmethod
    def _parse_metadata(comment_line: str) -> Dict[str, float]:
        """
        Parses labeled numeric metadata fields from the comment line.
        For example, "coord= 25.38", "time=10.5", etc.

        Returns a dictionary {<field>: <float_value>}.

        Note: We do NOT parse 'energy' or 'e' here, to avoid conflicting
              with the special single-float energy logic.
        """
        results = {}
        if not comment_line.strip():
            return results

        # Known label variants for reaction coordinate or time, etc.
        known_labels = {
            "coord": ["coord", "coordinate"],
            "time": ["t", "time"],
        }
        # You can expand 'known_labels' as needed, e.g. add "step", "temp", etc.

        # Build one combined pattern: label + '=' + float
        # We'll do something like: (?i)\b(coord|coordinate|t|time)\s*[:=]\s*([float])
        label_group = "|".join(
            label for variants in known_labels.values() for label in variants
        )
        pattern = (
            rf"(?i)\b({label_group})\s*[:=]\s*"
            r"([-+]?\d+(\.\d+)?([eE][+-]?\d+)?)"
        )

        for match in re.finditer(pattern, comment_line):
            raw_label = match.group(1).lower()  # e.g. "coord", "coordinate", "t", etc.
            value_str = match.group(2)         # e.g. "25.38"
            # Map raw_label to our canonical key
            for canonical_field, variants in known_labels.items():
                if raw_label in variants:
                    results[canonical_field] = float(value_str)
                    break

        return results


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

    # Print bond matrix
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

    # Print new metadata fields
    print("\nMetadata:")
    print(f"  XYZ_Comment:  {mol.XYZ_Comment}")
    print(f"  Energy:       {mol.Energy}")
    print(f"  file_name:    {mol.file_name}")
    print(f"  frame_number: {mol.frame_number}")
    print(f"  extra metadata dictionary: {mol.metadata}")
