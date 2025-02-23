# bond_manager.py
"""
This module encapsulates the interactive logic for cycling bond types and toggling bonds.
All bond visual properties are defined in config.py and used by the bond classes.
We now use two cycle lists:
  - A unified cycle that includes removal (None) for toggling via uppercase 'B'.
  - A non-removal cycle for cycling bond types with lowercase 'b'.
"""

from bonds import CovalentBond, NCIBond, TS_Bond, Distance_Bond

# Unified cycle (including removal) used for toggle operations.
UNIFIED_BOND_CYCLE = [Distance_Bond, CovalentBond, NCIBond, TS_Bond, None]

# Non-removal cycle used for cycling bond types with lowercase 'b'.
NON_REMOVAL_BOND_CYCLE = [Distance_Bond, CovalentBond, NCIBond, TS_Bond]


def get_next_bond_type(current_bond, cycle_list):
    """
    Given the current bond (or None) and a cycle list, return the next bond type.
    
    Parameters:
      - current_bond: An instance of a bond or None.
      - cycle_list: A list of bond classes and/or None.
    
    Returns:
      - The next bond class in the cycle list.
    """
    current_type = None if current_bond is None else type(current_bond)
    try:
        idx = cycle_list.index(current_type)
    except ValueError:
        idx = 0
    next_idx = (idx + 1) % len(cycle_list)
    return cycle_list[next_idx]


def cycle_existing_bond(bond, molecule, cycle_list=NON_REMOVAL_BOND_CYCLE):
    """
    Cycle the type of an existing bond within the molecule using a non-removal cycle.
    This function never returns None.
    
    Parameters:
      - bond: The current bond instance.
      - molecule: The ORTEP_Molecule instance managing bonds.
      - cycle_list: The list of bond types to cycle through (excluding None).
    
    Returns:
      - The new bond instance after cycling.
    """
    next_bond_class = get_next_bond_type(bond, cycle_list)
    # Replace the bond with the new type.
    new_bond = molecule.replace_bond(bond, next_bond_class)
    return new_bond


def cycle_atom_pair(atom1, atom2, molecule, cycle_list=NON_REMOVAL_BOND_CYCLE):
    """
    Cycle the bond type between two atoms in the molecule using a non-removal cycle.
    If a bond already exists between atom1 and atom2, cycle its type.
    If no bond exists, create one using the first type in the cycle.
    
    Parameters:
      - atom1, atom2: The two ORTEP_Atom objects.
      - molecule: The ORTEP_Molecule instance.
      - cycle_list: The list of bond types to cycle through (excluding None).
    
    Returns:
      - The new or updated bond instance, or None if no bond is created.
    """
    existing_bond = molecule.get_bond(atom1, atom2)
    if existing_bond:
        return cycle_existing_bond(existing_bond, molecule, cycle_list)
    else:
        # Create a new bond using the first type in the cycle.
        first_bond_class = cycle_list[0]
        return molecule.create_bond(atom1, atom2, first_bond_class)


def toggle_bond(atom1, atom2, molecule, default_bond=Distance_Bond):
    """
    Toggle the bond between two atoms.
      - If a bond exists, remove it.
      - If no bond exists, create one using the default_bond type.
    
    Parameters:
      - atom1, atom2: The two ORTEP_Atom objects.
      - molecule: The ORTEP_Molecule instance.
      - default_bond: The bond class to use when creating a new bond.
    
    Returns:
      - The new bond instance if created, or None if the bond was removed.
    """
    existing_bond = molecule.get_bond(atom1, atom2)
    if existing_bond:
        molecule.remove_bond(existing_bond)
        return None
    else:
        return molecule.create_bond(atom1, atom2, default_bond)
