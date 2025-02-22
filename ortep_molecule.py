# ortep_molecule.py

class ORTEP_Atom:
    """
    Represents a single atom in 3D space.
    Stores the element symbol and its x, y, z coordinates.
    """
    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.x = x
        self.y = y
        self.z = z

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
        self.atoms.append(atom)

    def add_bond(self, bond):
        self.bonds.append(bond)
