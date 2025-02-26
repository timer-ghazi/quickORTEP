#!/usr/bin/env python3
"""
trajectory_manager.py

This module provides the _TrajectoryManager class, which encapsulates the logic
for converting a trajectory frame (from a MoleculeWithNCIs) into an ORTEP_Molecule.
It handles:
  - Clamping the frame index.
  - Creating ORTEP_Atom objects.
  - Rebuilding bonds from a bond matrix (if available).
  - Adding non-covalent interactions.
  - Optionally filtering out hydrogen atoms attached to carbons.
"""

from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from bonds import CovalentBond, NCIBond

class _TrajectoryManager:
    def __init__(self, trajectory):
        """
        Initialize the trajectory manager with the given trajectory.

        Parameters:
            trajectory: An object that provides get_frame() and has _raw_frames.
        """
        self.trajectory = trajectory

    def convert_frame(self, frame_index, show_hydrogens=True):
        """
        Convert a trajectory frame into an ORTEP_Molecule.
        
        Parameters:
            frame_index (int): The index of the frame to convert.
            show_hydrogens (bool): Whether to include hydrogen atoms attached to carbons.
        
        Returns:
            new_ortep_mol (ORTEP_Molecule): The molecule converted from the trajectory frame.
        """
        # Clamp the frame index to a valid range.
        total_frames = len(self.trajectory._raw_frames)
        frame_index = max(0, min(frame_index, total_frames - 1))

        # Extract the frame from the trajectory.
        new_mol_nci = self.trajectory.get_frame(frame_index)

        # Create a new ORTEP_Molecule and add atoms.
        new_ortep_mol = ORTEP_Molecule()
        for atom in new_mol_nci.atoms:
            a = ORTEP_Atom(atom.symbol, atom.x, atom.y, atom.z)
            new_ortep_mol.add_atom(a)

        # Rebuild bonds from the bond matrix if available.
        n_atoms = len(new_mol_nci.atoms)
        bondmat = getattr(new_mol_nci, "bond_matrix", None)
        if bondmat is not None:
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    if bondmat[i, j] == 1:
                        b = CovalentBond(new_ortep_mol.atoms[i], new_ortep_mol.atoms[j])
                        new_ortep_mol.add_bond(b)

        # Add non-covalent interactions.
        for (i, j), interactions in new_mol_nci.ncis.items():
            nci_bond = NCIBond(new_ortep_mol.atoms[i], new_ortep_mol.atoms[j])
            new_ortep_mol.add_bond(nci_bond)

        # Apply hydrogen filtering if requested.
        if not show_hydrogens:
            self._apply_hydrogen_filter(new_ortep_mol)

        return new_ortep_mol

    def _apply_hydrogen_filter(self, mol):
        """
        Remove hydrogen atoms attached to carbons and their associated bonds.
        After removal, reassign atom indices.
        
        Parameters:
            mol (ORTEP_Molecule): The molecule to filter.
        """
        atoms_to_remove = []
        bonds_to_remove = []
        for bond in mol.bonds:
            if (bond.atom1.symbol == "H" and bond.atom2.symbol == "C"):
                if bond.atom1 not in atoms_to_remove:
                    atoms_to_remove.append(bond.atom1)
                bonds_to_remove.append(bond)
            elif (bond.atom2.symbol == "H" and bond.atom1.symbol == "C"):
                if bond.atom2 not in atoms_to_remove:
                    atoms_to_remove.append(bond.atom2)
                bonds_to_remove.append(bond)
        for bond in bonds_to_remove:
            if bond in mol.bonds:
                mol.bonds.remove(bond)
        for atom in atoms_to_remove:
            if atom in mol.atoms:
                mol.atoms.remove(atom)
        # Reassign atom indices after filtering.
        for idx, atom in enumerate(mol.atoms, start=1):
            atom.index = idx
