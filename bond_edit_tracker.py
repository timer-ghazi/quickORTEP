#!/usr/bin/env python3
"""
bond_edit_tracker.py

This module provides the BondEditTracker class, which tracks bond edits made
by the user and applies them across different frames in a trajectory.
"""

from bonds import CovalentBond, NCIBond, TS_Bond, Distance_Bond

class BondEditTracker:
    """
    Keeps track of user-defined bond modifications and applies them to molecules
    across trajectory frames.
    """
    def __init__(self):
        """
        Initialize the bond edit tracker.
        
        The edits dictionary maps atom index pairs to (action, bond_class) tuples:
        - action: 'add', 'remove', or 'change'
        - bond_class: The bond class to use for 'add' or 'change' actions
        """
        self.edits = {}  # Maps (atom1_idx, atom2_idx) -> (action, bond_class)
        self.enabled = True  # Whether to apply edits
    
    def add_bond(self, atom1_idx, atom2_idx, bond_class):
        """
        Record a bond addition.
        
        Parameters:
            atom1_idx, atom2_idx (int): Atom indices
            bond_class: Bond class to use
        """
        key = (min(atom1_idx, atom2_idx), max(atom1_idx, atom2_idx))
        self.edits[key] = ("add", bond_class)
    
    def remove_bond(self, atom1_idx, atom2_idx):
        """
        Record a bond removal.
        
        Parameters:
            atom1_idx, atom2_idx (int): Atom indices
        """
        key = (min(atom1_idx, atom2_idx), max(atom1_idx, atom2_idx))
        self.edits[key] = ("remove", None)
    
    def change_bond_type(self, atom1_idx, atom2_idx, bond_class):
        """
        Record a bond type change.
        
        Parameters:
            atom1_idx, atom2_idx (int): Atom indices
            bond_class: New bond class
        """
        key = (min(atom1_idx, atom2_idx), max(atom1_idx, atom2_idx))
        self.edits[key] = ("change", bond_class)
    
    def clear_edits(self):
        """Clear all recorded bond edits."""
        self.edits.clear()
    
    def get_edit_count(self):
        """Return the number of bond edits recorded."""
        return len(self.edits)
    
    def toggle(self):
        """Toggle whether edits are applied."""
        self.enabled = not self.enabled
        return self.enabled
    
    def apply_edits(self, molecule):
        """
        Apply all stored edits to the given molecule.
        
        Parameters:
            molecule (ORTEP_Molecule): The molecule to modify
        """
        if not self.enabled:
            return
        
        for (atom1_idx, atom2_idx), (action, bond_class) in self.edits.items():
            # Get the atoms by index
            atom1 = next((a for a in molecule.atoms if a.index == atom1_idx), None)
            atom2 = next((a for a in molecule.atoms if a.index == atom2_idx), None)
            
            if not atom1 or not atom2:
                continue  # Skip if atoms don't exist in this frame
            
            existing_bond = molecule.get_bond(atom1, atom2)
            
            if action == "remove" and existing_bond:
                molecule.remove_bond(existing_bond)
            elif action == "add" and not existing_bond:
                molecule.create_bond(atom1, atom2, bond_class)
            elif action == "change" and existing_bond and type(existing_bond) != bond_class:
                molecule.replace_bond(existing_bond, bond_class)