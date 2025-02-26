#!/usr/bin/env python3
"""
main.py

Usage:
  python main.py <xyz_file> [ss_factor=1] [tile_size=128]

This replicates the old ortep-view-v2.py behavior using the refactored architecture,
which now includes a multi-line HUD for displaying interactive information and trajectory navigation.
"""

import sys
from molecule_nci import MoleculeWithNCIs
from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from bonds import CovalentBond, NCIBond
from ortep_viewer import MoleculeViewer
from trajectory import Trajectory

def main():
    if len(sys.argv) < 2:
        print("Usage: main.py <xyz_file> [ss_factor=1] [tile_size=128]")
        sys.exit(1)

    xyz_file = sys.argv[1]
    # Default to no anti-aliasing => ss_factor=1
    ss_factor = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    # Default tile_size = 128
    tile_size = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

    # --- Load the full trajectory ---
    traj = Trajectory.from_xyz_file(xyz_file)
    # Extract the first frame to build our initial molecule.
    mol_nci = traj.get_frame(0)

    # --- Build the ORTEP_Molecule from the first frame ---
    ortep_mol = ORTEP_Molecule()

    # Create ORTEP_Atom objects.
    for atom in mol_nci.atoms:
        a = ORTEP_Atom(atom.symbol, atom.x, atom.y, atom.z)
        ortep_mol.add_atom(a)

    # Convert the bond matrix to bonds using the CovalentBond class.
    n_atoms = len(mol_nci.atoms)
    bondmat = getattr(mol_nci, "bond_matrix", None)
    if bondmat is not None:
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if bondmat[i, j] == 1:
                    b = CovalentBond(ortep_mol.atoms[i], ortep_mol.atoms[j])
                    ortep_mol.add_bond(b)

    # Add non-covalent interactions.
    for (i, j), interactions in mol_nci.ncis.items():
        nci_bond = NCIBond(ortep_mol.atoms[i], ortep_mol.atoms[j])
        ortep_mol.add_bond(nci_bond)

    # --- Launch the viewer with trajectory support ---
    viewer = MoleculeViewer(ortep_mol, width=700, height=700,
                              ss_factor=ss_factor, tile_size=tile_size)
    # Pass the trajectory to the viewer.
    viewer.trajectory = traj
    viewer.current_frame = 0
    viewer.total_frames = len(traj._raw_frames)
    
    viewer.fit_molecule_to_window()  # Adjust zoom and centering so the molecule fits
    viewer.run()

if __name__ == "__main__":
    main()
