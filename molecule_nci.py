#!/usr/bin/env python3
# molecule_nci.py
#
# A subclass "MoleculeWithNCIs" that detects non-covalent interactions
# (H-bonds, sigma-hole interactions, steric clashes), prints results
# in a formatted table, uses 1-based atom numbering, can handle multiple XYZ files
# from the command line, and supports debug output.

import sys
import math
import numpy as np

from molecule import Molecule
from elements_table import Elements

# This list will store tuples of (priority, function)
_NCI_METHODS = []

def register_nci(priority=0):
    """
    Decorator for registering an NCI detection function with a given priority.
    Lower 'priority' values run first.
    """
    def decorator(func):
        _NCI_METHODS.append((priority, func))
        _NCI_METHODS.sort(key=lambda x: x[0])
        return func
    return decorator

# A unified set of acceptable "sigma-hole" acceptors for halogen/chalcogen
VALID_SIGMA_HOLE_ACCEPTORS = {"N", "O", "F", "S", "CL", "BR", "I"}

class MoleculeWithNCIs(Molecule):
    """
    Subclass of Molecule that adds Non-Covalent Interaction (NCI) detection.
    """
    def __init__(self, title: str = ""):
        super().__init__(title=title)
        # Dictionary mapping (atom_i, atom_j) -> list of interaction records
        # Each record is a dict with keys like: type, distance, angle, angle_atoms, intra_or_inter
        self.ncis = {}

    def detect_all_ncis(self, run_steric_clashes: bool = False, debug: bool = False):
        """
        Calls each registered NCI detection method in order of ascending priority.
        If run_steric_clashes=False, the steric clash detection method is skipped.
        If debug=True, prints some debug info.
        """
        if debug:
            print("[DEBUG] detect_all_ncis: run_steric_clashes =", run_steric_clashes)
        for priority, detect_func in _NCI_METHODS:
            # Skip steric clash detection if requested.
            if detect_func.__name__ == "detect_steric_clashes" and (not run_steric_clashes):
                if debug:
                    print("[DEBUG] Skipping steric_clashes.")
                continue
            if debug:
                print(f"[DEBUG] Calling {detect_func.__name__} (priority={priority})")
            detect_func(self, debug=debug)

    def _vectorized_angle(self, donor_idx: int, central_idx: int, candidate_indices: np.ndarray) -> np.ndarray:
        """
        Vectorized calculation of angles at the central atom for given donor and candidate acceptor indices.
        
        Returns:
            angles in degrees as a numpy array.
        """
        # Get coordinates
        r_d = np.array([self.atoms[donor_idx].x, self.atoms[donor_idx].y, self.atoms[donor_idx].z])
        r_c = np.array([self.atoms[central_idx].x, self.atoms[central_idx].y, self.atoms[central_idx].z])
        # Build candidate coordinate array
        r_candidates = np.array([[self.atoms[j].x, self.atoms[j].y, self.atoms[j].z] for j in candidate_indices])
        
        # Vectors from central to donor and candidates
        v_cd = r_d - r_c  # shape (3,)
        v_ca = r_candidates - r_c  # shape (n, 3)
        
        # Normalize
        norm_cd = np.linalg.norm(v_cd)
        norm_ca = np.linalg.norm(v_ca, axis=1)
        # Dot products
        dots = np.dot(v_ca, v_cd)
        # Compute cosine values (avoid division by zero)
        cos_angles = dots / (norm_ca * norm_cd + 1e-8)
        # Clip for numerical safety and compute angles in degrees
        cos_angles = np.clip(cos_angles, -1.0, 1.0)
        angles = np.degrees(np.arccos(cos_angles))
        return angles

    @register_nci(priority=0)
    def detect_hydrogen_bonds(self,
                              max_dist=2.5,
                              angle_cutoff=160.0,
                              donors=("N", "O", "F"),
                              acceptors=("N", "O", "F"),
                              debug=False):
        """
        Simple hydrogen-bond detection:
        1) Atom i must be 'H' covalently bonded to a 'donor' (N, O, F by default).
        2) Atom j must be in 'acceptors' (N, O, F).
        3) dist(H--acceptor) <= max_dist
        4) angle(donor-H-acceptor) >= angle_cutoff
        """
        n = len(self.atoms)
        if debug:
            print("[DEBUG] detect_hydrogen_bonds")
        
        # Precompute full distance matrix once
        dist_matrix = self.compute_distance_matrix()
        
        # Identify indices for hydrogen atoms
        for i in range(n):
            if self.atoms[i].symbol.upper() != "H":
                continue

            # Determine donor: from bond matrix row of hydrogen
            donor_candidates = np.where(self.bond_matrix[i] > 0)[0]
            donor_idx = None
            for d in donor_candidates:
                if self.atoms[d].symbol.upper() in [d.upper() for d in donors]:
                    donor_idx = d
                    break
            if donor_idx is None:
                continue

            # Identify candidate acceptors: not hydrogen itself, not the donor, and with acceptable element.
            candidate_indices = []
            for j in range(n):
                if j == i or j == donor_idx:
                    continue
                if self.atoms[j].symbol.upper() not in [a.upper() for a in acceptors]:
                    continue
                candidate_indices.append(j)
            if not candidate_indices:
                continue
            candidate_indices = np.array(candidate_indices, dtype=int)

            # Vectorized filtering by distance
            dists = dist_matrix[i, candidate_indices]
            valid_mask = dists <= max_dist
            if not np.any(valid_mask):
                continue
            valid_candidates = candidate_indices[valid_mask]
            valid_dists = dists[valid_mask]

            # Compute angles vectorized: angle(donor, H, candidate)
            angles = self._vectorized_angle(donor_idx, i, valid_candidates)
            # Further filter based on angle cutoff
            for idx, ang in zip(valid_candidates, angles):
                if ang >= angle_cutoff:
                    pair = (min(i, idx), max(i, idx))
                    if pair not in self.ncis:
                        self.ncis[pair] = []
                    self.ncis[pair].append({
                        "type": "H-bond",
                        "distance": dist_matrix[i, idx],
                        "angle": ang,
                        "angle_atoms": (donor_idx, i, idx),
                        "intra_or_inter": self._intra_or_inter(donor_idx, idx)
                    })
                    if debug:
                        print(f"[DEBUG] H-bond: i={i}, donor={donor_idx}, j={idx}, "
                              f"dist={dist_matrix[i, idx]:.2f}, angle={ang:.1f}")

    @register_nci(priority=1)
    def detect_sigma_hole_bonds(self, debug=False):
        """
        Halogen + Chalcogen detection using one acceptor set: VALID_SIGMA_HOLE_ACCEPTORS.
        We check angle_min=120 (halogen) or 130 (chalcogen),
        distance < 0.9*(vdw_i+vdw_j).
        """
        groups = [
            {"label": "halogen_bond",   "symbols": {"CL", "BR", "I"},  "angle_min": 120.0},
            {"label": "chalcogen_bond", "symbols": {"S", "SE", "TE"},  "angle_min": 130.0},
        ]
        fraction_vdw = 0.9
        n = len(self.atoms)
        dist_matrix = self.compute_distance_matrix()

        for group in groups:
            bond_type = group["label"]
            donor_symbols = {s.upper() for s in group["symbols"]}
            angle_min = group["angle_min"]

            if debug:
                print(f"[DEBUG] Checking {bond_type} with donors={donor_symbols}")

            # Find candidate donor indices for this group
            for i in range(n):
                sym_i = self.atoms[i].symbol.upper()
                if sym_i not in donor_symbols:
                    continue
                vdw_i = Elements.vdw_radius(sym_i)
                # Get neighbors from bond matrix
                neighbors_i = np.where(self.bond_matrix[i] > 0)[0]
                # Identify potential acceptors: not covalently bonded, not self.
                candidate_indices = []
                for j in range(n):
                    if j == i or j in neighbors_i:
                        continue
                    sym_j = self.atoms[j].symbol.upper()
                    if sym_j not in VALID_SIGMA_HOLE_ACCEPTORS:
                        continue
                    candidate_indices.append(j)
                if not candidate_indices:
                    continue
                candidate_indices = np.array(candidate_indices, dtype=int)
                dists = dist_matrix[i, candidate_indices]
                # Get vdW radii for candidate acceptors
                vdw_j = np.array([Elements.vdw_radius(self.atoms[j].symbol.upper()) for j in candidate_indices])
                thresholds = fraction_vdw * (vdw_i + vdw_j)
                valid_mask = dists < thresholds
                if not np.any(valid_mask):
                    continue
                valid_candidates = candidate_indices[valid_mask]
                valid_dists = dists[valid_mask]
                
                # For each valid candidate, check if any neighbor of i gives an angle above cutoff.
                for idx, dist_val in zip(valid_candidates, valid_dists):
                    best_angle = None
                    best_neighbor = None
                    # Check each neighbor for the angle (angle(neighbor, i, candidate))
                    for neighbor in neighbors_i:
                        ax = self.angle(neighbor, i, idx)
                        if ax > angle_min:
                            best_angle = ax
                            best_neighbor = neighbor
                            break
                    if best_angle is not None:
                        pair = (min(i, idx), max(i, idx))
                        if pair not in self.ncis:
                            self.ncis[pair] = []
                        self.ncis[pair].append({
                            "type": bond_type,
                            "distance": dist_val,
                            "angle": best_angle,
                            "angle_atoms": (best_neighbor, i, idx),
                            "intra_or_inter": self._intra_or_inter(i, idx)
                        })
                        if debug:
                            print(f"[DEBUG] {bond_type}: i={i}, j={idx}, "
                                  f"dist={dist_val:.2f}, angle={best_angle:.1f}")

    @register_nci(priority=999)
    def detect_steric_clashes(self,
                              clash_tolerance=0.4,
                              only_hydrogen_clashes=True,
                              ignore_same_neighbor=True,
                              debug=False):
        """
        Check for steric clashes. Default is to look only for
        hydrogen-hydrogen clashes unless 'only_hydrogen_clashes=False'.
        """
        n = len(self.atoms)
        if debug:
            print("[DEBUG] detect_steric_clashes")
        dist_matrix = self.compute_distance_matrix()
        # Precompute vdW radii array (upper-case symbols)
        vdw_radii = np.array([Elements.vdw_radius(atom.symbol.upper()) for atom in self.atoms])
        
        # Use upper triangle indices (i<j)
        triu_indices = np.triu_indices(n, k=1)
        i_inds, j_inds = triu_indices
        
        # If only hydrogen clashes, filter pairs where both atoms are H
        if only_hydrogen_clashes:
            mask = np.array([self.atoms[i].symbol.upper() == "H" and self.atoms[j].symbol.upper() == "H" 
                             for i, j in zip(i_inds, j_inds)])
            i_inds = i_inds[mask]
            j_inds = j_inds[mask]
        
        # Skip covalently bonded pairs
        bond_mask = self.bond_matrix[i_inds, j_inds] == 0
        i_inds = i_inds[bond_mask]
        j_inds = j_inds[bond_mask]
        
        # Compute clash threshold for each pair
        thresholds = vdw_radii[i_inds] + vdw_radii[j_inds] - clash_tolerance
        dists = dist_matrix[i_inds, j_inds]
        clash_mask = dists < thresholds
        
        # Process each candidate pair that passes vectorized filters
        for i, j, d in zip(i_inds[clash_mask], j_inds[clash_mask], dists[clash_mask]):
            # Check if they share a neighbor if required
            if ignore_same_neighbor and self._share_a_neighbor(i, j):
                continue
            # Avoid duplicating interactions if another NCI exists already.
            pair = (i, j)
            existing = (self.ncis.get(pair, []) or self.ncis.get((j, i), []))
            skip_clash = any(x["type"] in ["H-bond", "halogen_bond", "chalcogen_bond"]
                             for x in existing)
            if not skip_clash:
                if pair not in self.ncis:
                    self.ncis[pair] = []
                self.ncis[pair].append({
                    "type": "steric_clash",
                    "distance": d,
                    "angle": None,
                    "angle_atoms": None,
                    "intra_or_inter": self._intra_or_inter(i, j)
                })
                if debug:
                    print(f"[DEBUG] steric_clash: i={i}, j={j}, dist={d:.2f}")

    def _share_a_neighbor(self, i: int, j: int) -> bool:
        """
        Returns True if atoms i and j share at least one common neighbor
        in the covalent bond matrix.
        """
        if self.bond_matrix is None:
            return False
        n = len(self.atoms)
        neighbors_i = [k for k in range(n) if k != i and self.bond_matrix[i, k] > 0]
        neighbors_j = [k for k in range(n) if k != j and self.bond_matrix[j, k] > 0]
        return len(set(neighbors_i).intersection(set(neighbors_j))) > 0

    def _intra_or_inter(self, i: int, j: int) -> str:
        """
        Returns "intra" if atoms i and j are in the same fragment,
        "inter" if different fragments, or "unknown" if fragments 
        haven't been identified.
        """
        if not self.fragments:
            return "unknown"
        frag_i = frag_j = None
        for f_id, members in self.fragments.items():
            if i in members:
                frag_i = f_id
            if j in members:
                frag_j = f_id
        if frag_i is None or frag_j is None:
            return "unknown"
        return "intra" if (frag_i == frag_j) else "inter"

    def list_interactions(self, interaction_type=None, fragment_scope=None):
        """
        Return a list of (pair, info_dict) for interactions that match
        a given interaction_type (e.g. "H-bond", "halogen_bond") 
        and/or a fragment_scope ("intra", "inter").
        """
        results = []
        for pair, all_info in self.ncis.items():
            for info in all_info:
                if interaction_type and info["type"] != interaction_type:
                    continue
                if fragment_scope and info["intra_or_inter"] != fragment_scope:
                    continue
                results.append((pair, info))
        return results

    def get_fragment_number(self, atom_idx: int) -> int:
        """
        Return the fragment number (1-based) for the given atom index.
        If not found, return 0.
        """
        if not self.fragments:
            return 0
        for f_id, members in self.fragments.items():
            if atom_idx in members:
                return f_id + 1  # switch to 1-based
        return 0


#
# Helper function(s) outside of the class
#
def get_fragments_for_interaction(mol: MoleculeWithNCIs, i: int, j: int, info: dict):
    """
    Returns a string describing which fragments are donor/acceptor 
    for an *intermolecular* interaction. If it's intramolecular, returns a single fragment label.
    """
    scope = info.get("intra_or_inter", "unknown")
    if scope != "inter":
        # intramolecular
        frag = mol.get_fragment_number(i)
        return f"Frag{frag}"

    # "inter" -> figure out donor/acceptor
    interaction_type = info["type"]
    donor_atom = i
    acceptor_atom = j

    if interaction_type == "H-bond":
        # angle_atoms = (donor_idx, H, acceptor)
        angle_atoms = info.get("angle_atoms", None)
        if angle_atoms and len(angle_atoms) == 3:
            donor_atom = angle_atoms[0]
            acceptor_atom = angle_atoms[2]

    elif interaction_type in ["halogen_bond", "chalcogen_bond"]:
        # angle_atoms = (neighbor, i, j)
        # i is the halogen/chalcogen, j is acceptor
        # so donor = i, acceptor = j
        angle_atoms = info.get("angle_atoms", None)
        if angle_atoms and len(angle_atoms) == 3:
            donor_atom = angle_atoms[1]
            acceptor_atom = angle_atoms[2]

    donor_frag = mol.get_fragment_number(donor_atom)
    acceptor_frag = mol.get_fragment_number(acceptor_atom)
    return f"Frag{donor_frag}->Frag{acceptor_frag}"


def print_interactions_table(mol: MoleculeWithNCIs, interactions):
    """
    Print a tabular summary of all interactions found, with:
      - Type
      - Pair (S1-O5, etc.)
      - Dist(Å)
      - Angle(°)
      - AngleAtoms (e.g. C1–S2–O6)
      - Fragments
      - Scope
    """
    print(f"{'Type':<14} {'Pair':<10} {'Dist(Å)':>8} {'Angle(°)':>9} {'AngleAtoms':<15} {'Fragments':<12} {'Scope':>6}")
    print("-" * 100)

    for (i, j), info in interactions:
        typ   = info["type"]
        dist  = info.get("distance", None)
        angle = info.get("angle", None)
        angle_atoms = info.get("angle_atoms", None)
        scope = info.get("intra_or_inter", "unknown")

        sym_i = mol.atoms[i].symbol
        sym_j = mol.atoms[j].symbol
        pair_str = f"{sym_i}{i+1}-{sym_j}{j+1}"

        dist_str  = f"{dist:.2f}"   if dist  is not None else "N/A"
        angle_str = f"{angle:.1f}" if angle is not None else "N/A"

        angle_atoms_str = "N/A"
        if angle_atoms and len(angle_atoms) == 3:
            a1, a2, a3 = angle_atoms
            sym_a1 = mol.atoms[a1].symbol
            sym_a2 = mol.atoms[a2].symbol
            sym_a3 = mol.atoms[a3].symbol
            angle_atoms_str = f"{sym_a1}{a1+1}-{sym_a2}{a2+1}-{sym_a3}{a3+1}"

        # label fragments
        fragments_str = get_fragments_for_interaction(mol, i, j, info)

        print(f"{typ:<14} {pair_str:<10} {dist_str:>8} {angle_str:>9} "
              f"{angle_atoms_str:<15} {fragments_str:<12} {scope:>6}")


def main():
    """
    Example main function to process one or more .xyz files, detect NCIs, and print results.
    """
    args = sys.argv[1:]
    debug_mode = False
    if "--debug" in args:
        debug_mode = True
        args.remove("--debug")

    if len(args) < 1:
        print("Usage: python molecule_nci.py <xyz_file1> [xyz_file2 ...] [--debug]")
        sys.exit(1)

    for xyz_file in args:
        print(f"\n=== Processing {xyz_file} ===")
        mol = MoleculeWithNCIs.from_xyz_data(file_name=xyz_file)

        # We now call detect_bonds from the *parent* class
        mol.detect_bonds(tolerance=0.3)  
        mol.find_fragments()

        # NCI detection
        mol.detect_all_ncis(run_steric_clashes=False, debug=debug_mode)

        print("\n--- Molecule Summary ---")
        print(mol.summary())

        all_ncis = mol.list_interactions()
        if not all_ncis:
            print("\nNo NCIs found.")
        else:
            print("\n--- Non-Covalent Interactions Detected ---")
            print_interactions_table(mol, all_ncis)


if __name__ == "__main__":
    main()
