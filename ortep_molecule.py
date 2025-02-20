class ORTEP_Atom:
    """
    Represents a single atom in 3D space.
    We'll store symbol, x, y, z. Other fields (e.g. color, radius) can be
    computed or cached later.
    """
    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.x = x
        self.y = y
        self.z = z

class ORTEP_Bond:
    """
    A simple bond class storing references to the two atoms and a 'bond_type'.
    For now, bond_type will be 'covalent' for all the standard adjacency.
    We can easily extend it for 'NCI', 'transition_state', etc.
    """
    def __init__(self, atom1, atom2, bond_type="covalent"):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_type = bond_type
        # We'll add more fields (like bond order, style, etc.) as needed.

class ORTEP_Molecule:
    """
    A container for atoms and bonds in the new "ORTEP" style.
    We'll build this from MoleculeWithNCIs or from other data sources.
    """
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def add_atom(self, atom):
        self.atoms.append(atom)

    def add_bond(self, bond):
        self.bonds.append(bond)
