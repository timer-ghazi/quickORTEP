#!/usr/bin/env python3
"""
quickORTEP.py

Usage:
  python quickORTEP.py <xyz_file> [ss_factor=1] [tile_size=128]

This replicates the old ortep-view-v2.py behavior using the refactored architecture,
which now includes a multi-line HUD for displaying interactive information and trajectory navigation.
"""

import sys
import os
import numpy as np
from molecule_nci import MoleculeWithNCIs
from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from bonds import CovalentBond, NCIBond
from ortep_viewer import MoleculeViewer
from trajectory import Trajectory
from config import DEFAULT_ENERGY_UNIT, ENERGY_UNITS

def main():
    if len(sys.argv) < 2:
        print("Usage: quickORTEP.py <xyz_file> [ss_factor=1] [tile_size=128]")
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
    n_bonds = 0
    if bondmat is not None:
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if bondmat[i, j] == 1:
                    b = CovalentBond(ortep_mol.atoms[i], ortep_mol.atoms[j])
                    ortep_mol.add_bond(b)
                    n_bonds += 1

    # Add non-covalent interactions.
    n_ncis = 0
    for (i, j), interactions in mol_nci.ncis.items():
        nci_bond = NCIBond(ortep_mol.atoms[i], ortep_mol.atoms[j])
        ortep_mol.add_bond(nci_bond)
        n_ncis += 1

    # --- Launch the viewer with trajectory support ---
    viewer = MoleculeViewer(ortep_mol, width=700, height=700,
                              ss_factor=ss_factor, tile_size=tile_size)
    
    # Pass the trajectory to the viewer using the proper method
    viewer.set_trajectory(traj)
    viewer.current_frame = 0  # This might be redundant now since set_trajectory sets total_frames
    
    # Log molecule and trajectory information
    filename = os.path.basename(xyz_file)
    viewer.message_service.log_info(f"Loaded {filename} with {n_atoms} atoms")
    viewer.message_service.log_info(f"Found {n_bonds} covalent bonds and {n_ncis} non-covalent interactions")
    
    # Log trajectory information
    if viewer.total_frames > 1:
        viewer.message_service.log_info(f"Trajectory contains {viewer.total_frames} frames")
        
        # Calculate energy statistics if available using enhanced energy_trajectory
        energies, energy_info = traj.energy_trajectory(
            convert_if_hartrees=True,
            convert_to_unit=DEFAULT_ENERGY_UNIT
        )
        
        valid_energies = energies[~np.isnan(energies)]
        if len(valid_energies) > 0:
            # Get unit symbol
            unit_symbol = ENERGY_UNITS.get(
                energy_info['converted_unit'], 
                {'symbol': energy_info['converted_unit']}
            )['symbol']
            
            min_e = np.min(valid_energies)
            max_e = np.max(valid_energies)
            min_idx = np.nanargmin(energies)
            max_idx = np.nanargmax(energies)
            
            # Check if data has energy method information
            method_info = ""
            if 'energy_data' in traj.metadata:
                energy_types = set()
                for frame_data in traj.metadata['energy_data'].values():
                    if 'type' in frame_data:
                        energy_types.add(frame_data['type'])
                if energy_types:
                    method_info = f" ({', '.join(energy_types)})"
            
            # Include normalization info if applicable
            norm_info = " (rel. to min)" if energy_info['normalized'] else ""
            
            viewer.message_service.log_info(
                f"Energy range{method_info}: {min_e:.4f} (frame {min_idx}) to {max_e:.4f} {unit_symbol}{norm_info}"
            )
    
    # Log rendering settings
    if ss_factor > 1:
        viewer.message_service.log_info(f"Anti-aliasing enabled (factor: {ss_factor}, tile size: {tile_size})")
    
    viewer.fit_molecule_to_window()  # Adjust zoom and centering so the molecule fits
    viewer.run()

if __name__ == "__main__":
    main()
