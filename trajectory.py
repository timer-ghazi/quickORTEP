#!/usr/bin/env python3
"""
trajectory.py
-------------
A module for handling molecular trajectories from XYZ files.
Supports both single structures and multi-frame trajectories.
Uses lazy processing: each frame is stored as a raw XYZ block and is only
converted into a fully processed MoleculeWithNCIs when needed.
"""

import os
import sys
import numpy as np
from molecule_nci import MoleculeWithNCIs  # Assumes molecule_nci.py (and molecule.py) are in the path

class Trajectory:
    def __init__(self, raw_frames):
        """
        Initialize a Trajectory instance.

        Parameters:
            raw_frames (list): List of raw XYZ data blocks.
                               Each block is a list of lines (strings).
        """
        self._raw_frames = raw_frames              # List of raw frame data (each a list of lines)
        self._molecule_cache = {}                  # Cache: frame_index -> MoleculeWithNCIs instance
        self._frame_energies = [None] * len(raw_frames)  # Will be filled on first full processing

    @classmethod
    def from_xyz_file(cls, file_name):
        """
        Load a trajectory from an XYZ file.
        The file can contain one or multiple concatenated XYZ frames.

        This method splits the file into frames and stores each raw frame.
        Energy extraction (and full parsing) is deferred until needed.

        Parameters:
            file_name (str): Path to the XYZ file.

        Returns:
            Trajectory: An instance of Trajectory.
        """
        if not os.path.exists(file_name):
            raise FileNotFoundError(f"File not found: {file_name}")

        with open(file_name, 'r') as f:
            data = f.read()

        # Split the file into non-empty lines
        lines = [line.strip() for line in data.splitlines() if line.strip()]
        frames = []
        i = 0
        n_lines = len(lines)

        # Try to interpret the file as a concatenated XYZ file.
        # Each frame is assumed to start with a line that is an integer atom count.
        while i < n_lines:
            try:
                num_atoms = int(lines[i])
            except ValueError:
                # If the very first line is not a number, assume the whole file is one frame.
                frames = [lines]
                break

            # A standard XYZ frame has: number line, comment line, then num_atoms lines.
            frame_end = i + num_atoms + 2
            if frame_end > n_lines:
                # Not enough lines: treat the remainder as one frame.
                frames.append(lines[i:])
                break
            frame_data = lines[i:frame_end]
            frames.append(frame_data)
            i = frame_end

        return cls(raw_frames=frames)

    def is_single_structure(self):
        """
        Return True if the trajectory contains only one frame.
        """
        return len(self._raw_frames) == 1

    def get_frame(self, frame_index):
        """
        Return the processed MoleculeWithNCIs instance for the given frame.
        Uses lazy evaluation with caching.

        Parameters:
            frame_index (int): Index of the desired frame.

        Returns:
            MoleculeWithNCIs: Processed molecule for the frame.
        """
        if frame_index not in self._molecule_cache:
            raw_data = self._raw_frames[frame_index]
            # Create the molecule using the robust parsing provided by read_xyz_data()
            mol = MoleculeWithNCIs.from_xyz_data(data=raw_data, frame_number=frame_index)
            # Run the expensive processing steps
            mol.to_standard_orientation()
            mol.detect_bonds()
            mol.find_fragments()
            mol.detect_all_ncis()
            self._molecule_cache[frame_index] = mol
            # Store the energy (if parsed) so we can print an energy table without reprocessing later.
            self._frame_energies[frame_index] = mol.Energy
        return self._molecule_cache[frame_index]

    def print_energy_table(self):
        """
        Print a table of frame numbers vs. energy.
        For each frame, the molecule is fully processed (if not already) so that
        its Energy attribute is available.
        """
        print("Frame\tEnergy")
        print("----------------")
        for i in range(len(self._raw_frames)):
            mol = self.get_frame(i)
            energy_str = f"{mol.Energy:.4f}" if mol.Energy is not None else "N/A"
            print(f"{i:3d}\t{energy_str}")

    def print_distance_table(self, atom1, atom2):
        """
        Print a table of the distance between two specified atoms for each frame.

        Parameters:
            atom1 (int): Index of the first atom.
            atom2 (int): Index of the second atom.
        """
        print(f"Frame\tDistance (Atom {atom1} - Atom {atom2})")
        print("-----------------------------------------")
        for i in range(len(self._raw_frames)):
            try:
                mol = self.get_frame(i)
                dist = mol.distance(atom1, atom2)
                print(f"{i:3d}\t{dist:8.3f} Å")
            except Exception as e:
                print(f"{i:3d}\tError: {e}")

    def distance_trajectory(self, atom1, atom2):
        """
        Return a NumPy array of distances between the specified atoms across all frames.
        
        Parameters:
            atom1 (int): Index of the first atom.
            atom2 (int): Index of the second atom.
            
        Returns:
            np.ndarray: Array of distances.
        """
        distances = []
        for i in range(len(self._raw_frames)):
            mol = self.get_frame(i)
            distances.append(mol.distance(atom1, atom2))
        return np.array(distances)

    def angle_trajectory(self, i, j, k, degrees=True):
        """
        Return a NumPy array of angles (atoms i, j, k) across all frames.
        
        Parameters:
            i, j, k (int): Atom indices.
            degrees (bool): Whether to return angles in degrees (default True).
            
        Returns:
            np.ndarray: Array of angles.
        """
        angles = []
        for frame in range(len(self._raw_frames)):
            mol = self.get_frame(frame)
            angles.append(mol.angle(i, j, k, degrees=degrees))
        return np.array(angles)

    def dihedral_trajectory(self, i, j, k, l, degrees=True):
        """
        Return a NumPy array of dihedral angles (atoms i, j, k, l) across all frames.
        
        Parameters:
            i, j, k, l (int): Atom indices.
            degrees (bool): Whether to return angles in degrees (default True).
            
        Returns:
            np.ndarray: Array of dihedral angles.
        """
        dihedrals = []
        for frame in range(len(self._raw_frames)):
            mol = self.get_frame(frame)
            dihedrals.append(mol.dihedral(i, j, k, l, degrees=degrees))
        return np.array(dihedrals)

    def summary(self):
        """
        Print summary information about the trajectory.
        """
        total_frames = len(self._raw_frames)
        print("Trajectory Summary")
        print("------------------")
        print(f"Total frames: {total_frames}")
        if self.is_single_structure():
            print("Single structure loaded.")
        else:
            print("Multi-frame trajectory.")
        # Force processing of frames to obtain energies (if available)
        energies = []
        for i in range(total_frames):
            mol = self.get_frame(i)
            if mol.Energy is not None:
                energies.append(mol.Energy)
        if energies:
            print(f"Energy range: {min(energies):.4f} to {max(energies):.4f}")
        else:
            print("No energy information available.")

def main():
    if len(sys.argv) < 2:
        print("Usage: python trajectory.py <xyz_file>")
        sys.exit(1)

    xyz_file = sys.argv[1]
    # Load the trajectory from the XYZ file.
    traj = Trajectory.from_xyz_file(xyz_file)

    # Print overall summary.
    print("\n--- Trajectory Summary ---")
    traj.summary()

    # Print an energy vs. frame table.
    print("\n--- Energy vs. Frame ---")
    traj.print_energy_table()

    # Print a distance vs. frame table for atoms 0 and 1.
    print("\n--- Distance vs. Frame (Atom 0 - Atom 1) ---")
    traj.print_distance_table(atom1=0, atom2=1)

    # Demonstrate trajectory-level geometric analysis:
    distances = traj.distance_trajectory(0, 1)
    print("\nDistance array (Atom 0 - Atom 1):")
    print(distances)

    # If there are at least 3 atoms, print an angle trajectory.
    try:
        angles = traj.angle_trajectory(0, 1, 2)
        print("\nAngle array (Atoms 0, 1, 2):")
        print(angles)
    except Exception as e:
        print("\nCould not compute angle (0, 1, 2):", e)

    # If there are at least 4 atoms, print a dihedral trajectory.
    try:
        dihedrals = traj.dihedral_trajectory(0, 1, 2, 3)
        print("\nDihedral array (Atoms 0, 1, 2, 3):")
        print(dihedrals)
    except Exception as e:
        print("\nCould not compute dihedral (0, 1, 2, 3):", e)

if __name__ == "__main__":
    main()
