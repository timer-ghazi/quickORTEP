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

        If the coordinate array is defined as C with shape (n, 3), the distance between any two atoms
        is given by:

        d_{ij} = sqrt( sum_{k=1}^{3} (C_{ik} - C_{jk})^2 )

        Returns:
            A numpy ndarray of shape (n, n) with distances.
        """
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
        return dist_matrix

    # ------------------------------------------------------------------------
    # NEW: Return the distance matrix as "Total Connections" representation
    # ------------------------------------------------------------------------
    def to_total_connections(self) -> np.ndarray:
        """
        Returns the NxN matrix of all interatomic distances (the 'total connections'),
        which can be viewed as a fully redundant coordinate representation.
        """
        return self.compute_distance_matrix()

    # ------------------------------------------------------------------------
    # NEW: Class method to build a Molecule from a distance matrix
    # ------------------------------------------------------------------------
    @classmethod
    def from_total_connections(cls,
                               dist_matrix: np.ndarray,
                               symbols: List[str],
                               title: str = "From TC",
                               allow_inconsistent: bool = False) -> "Molecule":
        """
        Creates a Molecule by embedding the given NxN pairwise distance matrix dist_matrix
        into 3D using classical multidimensional scaling (a.k.a. distance geometry).

        Args:
            dist_matrix (np.ndarray):
                NxN array of pairwise distances.
            symbols (List[str]):
                List of length N with atomic symbols (e.g. ["C", "H", "H", "H"]).
            title (str):
                Title or label for the resulting molecule.
            allow_inconsistent (bool):
                If True, will proceed with an approximate 3D embedding even if the distances
                are not perfectly consistent with a Euclidean 3D geometry (negative eigenvalues, etc.).
                If False, will raise an error if the distance matrix is inconsistent with 3D.

        Returns:
            A new Molecule instance with the embedded coordinates.

        Raises:
            ValueError if dist_matrix is not NxN or if distances are inconsistent with 3D embedding
            (unless allow_inconsistent=True).
        """
        N = dist_matrix.shape[0]
        if dist_matrix.shape[1] != N:
            raise ValueError("dist_matrix must be square (NxN).")

        if len(symbols) != N:
            raise ValueError(f"symbols must have length {N}, but got {len(symbols)}.")

        # 1. Square the distances
        D2 = dist_matrix ** 2

        # 2. Double-center to form Gram matrix B
        J = np.eye(N) - np.ones((N, N)) / N
        B = -0.5 * (J @ D2 @ J)

        # 3. Eigen-decomposition
        eigvals, eigvecs = np.linalg.eigh(B)
        # Sort eigenvalues in descending order
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        # If we require consistent 3D embedding, check the 3rd eigenvalue
        if not allow_inconsistent:
            if eigvals[2] <= 0:
                raise ValueError("Distance matrix not embeddable in 3D (3rd eigenvalue <= 0).")

        # 4. Build 3D coords from top 3 eigenvalues/eigenvectors
        #    If allow_inconsistent, we clamp negative eigenvalues to zero
        top3 = np.maximum(eigvals[:3], 0.0)
        L = np.diag(np.sqrt(top3))
        V = eigvecs[:, :3]
        coords_3d = V @ L  # shape (N, 3)

        # 5. Build a new Molecule and populate coordinates.
        atom_data = []
        for i in range(N):
            symbol = symbols[i]
            x, y, z = coords_3d[i]
            atom_data.append((symbol, x, y, z))

        new_mol = cls.from_atoms(atom_data, title=title)
        return new_mol

    def detect_bonds(self, tolerance: float = 0.3):
        """
        Fills self.bond_matrix using a vectorized approach by comparing all interatomic distances
        to the sum of covalent radii (plus a tolerance).

        Parameters:
            tolerance: Additional margin (in Å) on top of covalent radius sums.
        """
        n = len(self.atoms)
        dist_matrix = self.compute_distance_matrix()

        # Create an array of covalent radii for each atom (assuming single bond order)
        radii = np.array([Elements.covalent_radius(atom.symbol, order="single") for atom in self.atoms])
        threshold_matrix = radii[:, None] + radii[None, :] + tolerance

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

                self.fragments[frag_id] = connected
                frag_id += 1

    def distance(self, i: int, j: int) -> float:
        """
        Returns the Euclidean distance between atoms i and j (in Å).
        """
        ax, ay, az = self.atoms[i].x, self.atoms[i].y, self.atoms[i].z
        bx, by, bz = self.atoms[j].x, self.atoms[j].y, self.atoms[j].z
        return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)

    def angle(self, i: int, j: int, k: int, degrees: bool = True) -> float:
        """
        Returns the angle at atom j formed by (i - j - k).
        If degrees=False, return radians.
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
        """
        r_i = np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])
        r_j = np.array([self.atoms[j].x, self.atoms[j].y, self.atoms[j].z])
        r_k = np.array([self.atoms[k].x, self.atoms[k].y, self.atoms[k].z])
        r_l = np.array([self.atoms[l].x, self.atoms[l].y, self.atoms[l].z])

        b1 = r_i - r_j
        b2 = r_k - r_j
        b3 = r_l - r_k

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        if np.linalg.norm(n1) == 0 or np.linalg.norm(n2) == 0:
            return 0.0

        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        angle_radians = np.arctan2(y, x)

        return np.degrees(angle_radians) if degrees else angle_radians

    def to_standard_orientation(self, use_atomic_numbers: bool = True) -> None:
        """
        Converts the molecule's atomic coordinates to a standard orientation via
        a moment of inertia alignment procedure.
        """
        coords = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

        if use_atomic_numbers:
            weights = np.array([Elements.atomic_number(atom.symbol) for atom in self.atoms])
        else:
            weights = np.array([Elements.mass(atom.symbol) for atom in self.atoms])

        coords_centered, T = self._compute_center(coords, weights)
        I = self._compute_inertia_tensor(coords_centered, weights)
        eigenvalues, eigenvectors = self._diagonalize_inertia_tensor(I)
        R = self._build_rotation_matrix(eigenvalues, eigenvectors)
        coords_std = self._apply_transformation(coords_centered, R)
        self._update_atom_coordinates(coords_std)

    def _compute_center(self, coords: np.ndarray, weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        total_weight = np.sum(weights)
        T = np.sum(weights[:, None] * coords, axis=0) / total_weight
        coords_centered = coords - T
        return coords_centered, T

    def _compute_inertia_tensor(self, coords_centered: np.ndarray, weights: np.ndarray) -> np.ndarray:
        I = np.zeros((3, 3))
        for i, (x, y, z) in enumerate(coords_centered):
            w = weights[i]
            I[0, 0] += w * (y**2 + z**2)
            I[1, 1] += w * (x**2 + z**2)
            I[2, 2] += w * (x**2 + y**2)
            I[0, 1] -= w * x * y
            I[0, 2] -= w * x * z
            I[1, 2] -= w * y * z
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        return I

    def _diagonalize_inertia_tensor(self, I: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        eigenvalues, eigenvectors = np.linalg.eigh(I)
        return eigenvalues, eigenvectors

    def _build_rotation_matrix(self, eigenvalues: np.ndarray, eigenvectors: np.ndarray) -> np.ndarray:
        idx = np.argsort(eigenvalues)
        v_small = eigenvectors[:, idx[0]]
        v_mid   = eigenvectors[:, idx[1]]
        v_large = eigenvectors[:, idx[2]]

        R = np.column_stack((v_large, v_mid, v_small))
        if np.linalg.det(R) < 0:
            R[:, 0] *= -1
        return R

    def _apply_transformation(self, coords: np.ndarray, R: np.ndarray) -> np.ndarray:
        return np.dot(coords, R)

    def _update_atom_coordinates(self, coords_std: np.ndarray) -> None:
        for i, atom in enumerate(self.atoms):
            atom.z, atom.y, atom.x = coords_std[i]

    def summary(self) -> str:
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

        Format:
            <number of atoms>
            <title>
            <atom symbol> <x> <y> <z>
        """
        lines = []
        lines.append(f"{len(self.atoms)}")
        comment = self.title
        lines.append(comment)
        for atom in self.atoms:
            lines.append(f"{atom.symbol:<5}{atom.x:17.12f}{atom.y:17.12f}{atom.z:17.12f}")
        return "\n".join(lines)

    @staticmethod
    def _parse_energy(comment_line: str) -> Optional[float]:
        if not comment_line.strip():
            return None
        labeled_pattern = r"(?i)\b(?:energy|e)\s*[:=]\s*([-+]?\d+(\.\d+)?([eE][+-]?\d+)?)"
        match = re.search(labeled_pattern, comment_line)
        if match:
            return float(match.group(1))

        all_floats = re.findall(
            r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
            comment_line
        )
        if len(all_floats) == 1:
            return float(all_floats[0])
        return None

    @staticmethod
    def _looks_like_just_energy(comment_line: str) -> bool:
        if not comment_line.strip():
            return True
        possible_energy = Molecule._parse_energy(comment_line)
        if possible_energy is None:
            return False
        all_floats = re.findall(
            r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
            comment_line
        )
        if len(all_floats) == 1:
            without_float = re.sub(
                r"[-+]?\d+\.\d+(?:[eE][+-]?\d+)?|[-+]?\d+(?:[eE][+-]?\d+)?",
                "",
                comment_line
            )
            if len(without_float.strip()) < 15:
                return True
        return False

    @staticmethod
    def _parse_metadata(comment_line: str) -> Dict[str, float]:
        """
        Parses labeled numeric metadata fields from the comment line.
        For example, "coord= 25.38", "time=10.5", etc.
        """
        results = {}
        if not comment_line.strip():
            return results

        known_labels = {
            "coord": ["coord", "coordinate"],
            "time": ["t", "time"],
        }
        label_group = "|".join(
            label for variants in known_labels.values() for label in variants
        )
        pattern = (
            rf"(?i)\b({label_group})\s*[:=]\s*"
            r"([-+]?\d+(\.\d+)?([eE][+-]?\d+)?)"
        )

        for match in re.finditer(pattern, comment_line):
            raw_label = match.group(1).lower()
            value_str = match.group(2)
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

    # ----------------------------------------------------------------------
    # NEW: Demonstration: Convert to TCs, embed back, and dump new XYZ
    # ----------------------------------------------------------------------
    print("\n--- Converting to Total Connections, then re-embedding in 3D ---")
    dist_mat = mol.to_total_connections()
    symbols = [atom.symbol for atom in mol.atoms]

    # Create a new Molecule from the distance matrix
    mol_reconstructed = Molecule.from_total_connections(dist_mat, symbols,
                                                        title="Reconstructed from TCs",
                                                        allow_inconsistent=False)
    print("\nReconstructed geometry from TCs:")
    print(mol_reconstructed.to_xyz())
