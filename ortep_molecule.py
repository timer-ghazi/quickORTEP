# ortep_molecule.py

class ORTEP_Atom:
    """
    Represents a single atom in 3D space.
    Stores the element symbol and its x, y, z coordinates.
    """
    def __init__(self, symbol, x, y, z, index=None):
        self.symbol = symbol
        self.x = x
        self.y = y
        self.z = z
        self.index = index
        self.selected = False  # Persist selection state

class ORTEP_Molecule:
    """
    A container for atoms and bonds in the ORTEP style.
    
    Note: The old ORTEP_Bond class is deprecated.
    Bonds should now be instances of classes from bonds.py.
    """
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def add_atom(self, atom):
        # Assign a sequential index (starting from 1)
        atom.index = len(self.atoms) + 1
        self.atoms.append(atom)

    def add_bond(self, bond):
        self.bonds.append(bond)

    def get_bond(self, atom1, atom2):
        """
        Return the bond connecting atom1 and atom2, if it exists.
        Order does not matter.
        """
        for bond in self.bonds:
            if ((bond.atom1 == atom1 and bond.atom2 == atom2) or
                (bond.atom1 == atom2 and bond.atom2 == atom1)):
                return bond
        return None

    def create_bond(self, atom1, atom2, bond_class):
        """
        Creates a new bond between atom1 and atom2 using the specified bond_class.
        
        Parameters:
          - atom1, atom2: The two ORTEP_Atom objects to bond.
          - bond_class: A bond class (e.g., COVALENT_BOND, NCI_BOND, TS_BOND, DISTANCE_BOND)
                        used to instantiate the bond.
        
        Returns:
          - The new bond instance if created, or None if a bond already exists.
        """
        if self.get_bond(atom1, atom2) is not None:
            print("Bond already exists between these atoms.")
            return None

        new_bond = bond_class(atom1, atom2)
        self.add_bond(new_bond)
        return new_bond

    def remove_bond(self, bond):
        """
        Safely remove the specified bond from the molecule.
        """
        if bond in self.bonds:
            self.bonds.remove(bond)
        else:
            print("Bond not found in the molecule.")

    def replace_bond(self, old_bond, new_bond_class):
        """
        Replace an existing bond with a new bond of type new_bond_class.
        Transfers selection state from old_bond to the new bond.
        
        Parameters:
          - old_bond: The bond instance to be replaced.
          - new_bond_class: The new bond class to instantiate.
        
        Returns:
          - The new bond instance if replacement is successful, or None otherwise.
        """
        if old_bond not in self.bonds:
            print("Old bond not found in the molecule.")
            return None

        # Retrieve the atoms from the old bond.
        atom1 = old_bond.atom1
        atom2 = old_bond.atom2

        # Remove the old bond.
        self.remove_bond(old_bond)

        # Create a new bond of the desired type.
        new_bond = new_bond_class(atom1, atom2)
        # Transfer selection state (or other needed properties).
        new_bond.selected = old_bond.selected

        self.add_bond(new_bond)
        return new_bond
