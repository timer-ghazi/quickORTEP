#!/usr/bin/env python3
"""
selection_manager.py

This module provides the _SelectionManager class to encapsulate the logic
for managing selections of atoms and bonds in the molecule viewer. It maintains
the selected atoms, bonds, and their corresponding identifiers for persistent
selection across trajectory frames.
"""

class _SelectionManager:
    def __init__(self):
        # Lists of selected atoms and bonds.
        self.selected_atoms = []
        self.selected_bonds = []
        # Persistent selection identifiers.
        self.selected_atom_indices = []
        self.selected_bond_ids = []

    def select_atom(self, atom, multi_select=False):
        """
        Select an atom. If multi_select is False, clear any existing selection.
        
        Parameters:
            atom: The atom object to be selected. Expected to have attributes
                  'selected' (boolean) and 'index' (int).
            multi_select (bool): Whether to add to the current selection.
        """
        if not multi_select:
            self.clear_selection()
        if atom not in self.selected_atoms:
            self.selected_atoms.append(atom)
            atom.selected = True
            if atom.index not in self.selected_atom_indices:
                self.selected_atom_indices.append(atom.index)

    def select_bond(self, bond, multi_select=False):
        """
        Select a bond. If multi_select is False, clear any existing selection.
        
        Parameters:
            bond: The bond object to be selected. Expected to have attributes
                  'selected', 'atom1', and 'atom2'. Each atom must have an 'index'.
            multi_select (bool): Whether to add to the current selection.
        """
        if not multi_select:
            self.clear_selection()
        if bond not in self.selected_bonds:
            self.selected_bonds.append(bond)
            bond.selected = True
            key = (min(bond.atom1.index, bond.atom2.index),
                   max(bond.atom1.index, bond.atom2.index))
            if key not in self.selected_bond_ids:
                self.selected_bond_ids.append(key)

    def clear_selection(self):
        """
        Clear the current selections for both atoms and bonds.
        """
        for atom in self.selected_atoms:
            atom.selected = False
        for bond in self.selected_bonds:
            bond.selected = False
        self.selected_atoms.clear()
        self.selected_bonds.clear()
        self.selected_atom_indices.clear()
        self.selected_bond_ids.clear()
