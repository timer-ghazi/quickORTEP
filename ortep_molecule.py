# ortep_molecule.py

from vectors import Vector, AxisSystem
from config import AXES_STYLE

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
    
    Also stores vectors for visualization of coordinate axes, normal modes,
    or other vector quantities.
    """
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.vectors = []
        self.axis_system = None  # Will hold coordinate axes if enabled

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
    
    def add_vector(self, vector):
        """
        Add a vector to the molecule.
        
        Parameters:
            vector: A Vector instance
        """
        self.vectors.append(vector)
        return vector
    
    def remove_vector(self, vector):
        """
        Remove a vector from the molecule.
        
        Parameters:
            vector: The Vector instance to remove
        """
        if vector in self.vectors:
            self.vectors.remove(vector)
        else:
            print("Vector not found in the molecule.")
    
    def clear_vectors(self):
        """
        Remove all vectors from the molecule.
        """
        self.vectors.clear()
    
    def create_coordinate_axes(self, origin=None):
        """
        Create or update a coordinate axis system at the specified origin.
        If origin is None, use the geometric center of the molecule.
        
        Parameters:
            origin: (x, y, z) coordinates for the axes origin, or None for molecule center
            
        Returns:
            The created AxisSystem instance
        """
        # If no origin specified, use geometric center of the molecule
        if origin is None:
            if not self.atoms:
                # No atoms, use (0,0,0)
                origin = (0.0, 0.0, 0.0)
            else:
                # Calculate geometric center
                x_sum = sum(atom.x for atom in self.atoms)
                y_sum = sum(atom.y for atom in self.atoms)
                z_sum = sum(atom.z for atom in self.atoms)
                n_atoms = len(self.atoms)
                origin = (x_sum / n_atoms, y_sum / n_atoms, z_sum / n_atoms)
        
        # Create the axis system
        self.axis_system = AxisSystem(origin=origin)
        
        # Clear any existing axis vectors and add the new ones
        self.vectors = [v for v in self.vectors if not hasattr(v, '_is_axis_vector')]
        for vector in self.axis_system.get_vectors():
            vector._is_axis_vector = True  # Mark as an axis vector
            self.add_vector(vector)
        
        return self.axis_system
    
    def toggle_coordinate_axes(self, show):
        """
        Show or hide the coordinate axes.
        
        Parameters:
            show: Boolean, whether to show (True) or hide (False) the axes
        """
        # Remove existing axes if present
        self.vectors = [v for v in self.vectors if not hasattr(v, '_is_axis_vector')]
        
        # Create new axes if showing
        if show:
            if self.axis_system:
                # Reuse existing axis system
                for vector in self.axis_system.get_vectors():
                    vector._is_axis_vector = True
                    self.add_vector(vector)
            else:
                # Create new axis system
                self.create_coordinate_axes()
        else:
            # Just clear the reference if hiding
            self.axis_system = None