#!/usr/bin/env python3
"""
trajectory.py
-------------
A module for handling molecular trajectories from various file formats.
Supports both single structures and multi-frame trajectories.
Uses lazy processing: each frame is stored as raw data and is only
converted into a fully processed MoleculeWithNCIs when needed.
"""

import os
import sys
import numpy as np
from abc import ABC, abstractmethod
from molecule_nci import MoleculeWithNCIs  # Assumes molecule_nci.py (and molecule.py) are in the path


class TrajectoryParser(ABC):
    """
    Abstract base class for trajectory parsers.
    Each parser class should implement methods to check file compatibility
    and parse file content into frames.
    """
    
    @classmethod
    @abstractmethod
    def can_parse(cls, file_path):
        """
        Determines if this parser can handle the given file
        
        Parameters:
            file_path (str): Path to the file to check
            
        Returns:
            bool: True if this parser can handle the file
        """
        pass
    
    @classmethod
    @abstractmethod
    def parse(cls, file_path):
        """
        Parse file into raw_frames and metadata
        
        Parameters:
            file_path (str): Path to the file to parse
            
        Returns:
            tuple: (raw_frames, metadata) where:
                  - raw_frames is a list of frame data blocks
                  - metadata is a dict with file-level metadata
        """
        pass
    
    @staticmethod
    def read_file_content(file_path):
        """Helper method to read file content"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
            
        with open(file_path, 'r') as f:
            return f.read()


class XYZTrajectoryParser(TrajectoryParser):
    """Parser for concatenated XYZ files"""
    
    @classmethod
    def can_parse(cls, file_path):
        """Check if file is in XYZ format"""
        # Simple check: file extension and basic structure
        if not file_path.lower().endswith('.xyz'):
            return False
            
        try:
            with open(file_path, 'r') as f:
                # Try to read first few lines
                first_line = next(f, '').strip()
                # Check if first line is a number (atom count)
                try:
                    num_atoms = int(first_line)
                    return True
                except ValueError:
                    return False
        except:
            return False
    
    @classmethod
    def parse(cls, file_path):
        """
        Parse XYZ file into raw frames and metadata.
        
        Returns:
            tuple: (raw_frames, metadata)
                  - raw_frames is a list of lists, each containing lines for one frame
                  - metadata contains file info like name and format
        """
        content = cls.read_file_content(file_path)
        
        # Extract metadata
        base_name = os.path.basename(file_path)
        if base_name.endswith('.xyz'):
            base_name = base_name[:-4]
            
        metadata = {
            'file_name': base_name,
            'format': 'xyz',
            'path': file_path
        }
        
        # Split the file into non-empty lines
        lines = [line.strip() for line in content.splitlines() if line.strip()]
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
        
        return frames, metadata


class Trajectory:
    """
    Class for handling molecular trajectories from various file formats.
    Uses lazy evaluation with caching for memory efficiency.
    """
    
    # Register all available parsers
    _parsers = [XYZTrajectoryParser]
    
    def __init__(self, raw_frames, metadata=None):
        """
        Initialize a Trajectory instance.

        Parameters:
            raw_frames (list): List of raw frame data blocks.
                              Each block format depends on the source.
            metadata (dict): Optional metadata about the trajectory
        """
        self._raw_frames = raw_frames              # List of raw frame data
        self._molecule_cache = {}                  # Cache: frame_index -> MoleculeWithNCIs instance
        self._frame_energies = [None] * len(raw_frames)  # Will be filled on processing
        self.metadata = metadata or {}             # Trajectory-level metadata

    @classmethod
    def register_parser(cls, parser_class):
        """
        Register a new parser class for additional file format support.
        
        Parameters:
            parser_class: Class that implements the TrajectoryParser interface
        """
        if parser_class not in cls._parsers:
            cls._parsers.append(parser_class)

    @classmethod
    def from_xyz_file(cls, file_name):
        """
        Load a trajectory from an XYZ file.
        Maintained for backward compatibility - calls from_file() internally.
        
        Parameters:
            file_name (str): Path to the XYZ file.
            
        Returns:
            Trajectory: An instance of Trajectory.
        """
        return cls.from_file(file_name)
        
    @classmethod
    def from_file(cls, file_path):
        """
        Load a trajectory from a file using the appropriate parser.
        
        Parameters:
            file_path (str): Path to the trajectory file.
            
        Returns:
            Trajectory: An instance of Trajectory.
            
        Raises:
            ValueError: If no parser can handle the given file format.
        """
        for parser_class in cls._parsers:
            if parser_class.can_parse(file_path):
                raw_frames, metadata = parser_class.parse(file_path)
                traj = cls(raw_frames=raw_frames, metadata=metadata)
                
                # Preload energies for all frames (for sorting, etc.)
                for i in range(len(raw_frames)):
                    traj.get_frame(i)
                return traj
                
        raise ValueError(f"Unsupported file format: {file_path}")

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

    def filter_frames(self, criteria_func):
        """
        Return a new Trajectory with frames that match the criteria.
        
        Parameters:
            criteria_func: Function that takes a MoleculeWithNCIs instance and returns True/False
            
        Returns:
            Trajectory: A new trajectory containing only matching frames
        """
        filtered_indices = []
        for i in range(len(self._raw_frames)):
            mol = self.get_frame(i)
            if criteria_func(mol):
                filtered_indices.append(i)
                
        filtered_frames = [self._raw_frames[i] for i in filtered_indices]
        return Trajectory(filtered_frames, metadata=self.metadata)

    def extract_frames(self, indices):
        """
        Extract specific frames by index.
        
        Parameters:
            indices (list): List of frame indices to extract
            
        Returns:
            Trajectory: A new trajectory containing only the specified frames
        """
        extracted_frames = [self._raw_frames[i] for i in indices if 0 <= i < len(self._raw_frames)]
        return Trajectory(extracted_frames, metadata=self.metadata)

    def merge(self, other_trajectory):
        """
        Merge with another trajectory.
        
        Parameters:
            other_trajectory (Trajectory): Another trajectory to merge with
            
        Returns:
            Trajectory: A new trajectory containing frames from both
        """
        new_frames = self._raw_frames + other_trajectory._raw_frames
        # Merge metadata (prioritize self's metadata for conflicts)
        merged_metadata = {**other_trajectory.metadata, **self.metadata}
        
        return Trajectory(new_frames, metadata=merged_metadata)

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
            
        # Print metadata
        if self.metadata:
            print("\nMetadata:")
            for key, value in self.metadata.items():
                print(f"  {key}: {value}")
            
        # Force processing of frames to obtain energies (if available)
        energies = []
        for i in range(total_frames):
            mol = self.get_frame(i)
            if mol.Energy is not None:
                energies.append(mol.Energy)
        if energies:
            print(f"\nEnergy range: {min(energies):.4f} to {max(energies):.4f}")
        else:
            print("\nNo energy information available.")


def main():
    if len(sys.argv) < 2:
        print("Usage: python trajectory.py <trajectory_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    
    # Load the trajectory from the file
    try:
        traj = Trajectory.from_file(file_path)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Print overall summary
    print("\n--- Trajectory Summary ---")
    traj.summary()

    # Print an energy vs. frame table
    print("\n--- Energy vs. Frame ---")
    traj.print_energy_table()

    # Print a distance vs. frame table for atoms 0 and 1
    print("\n--- Distance vs. Frame (Atom 0 - Atom 1) ---")
    traj.print_distance_table(atom1=0, atom2=1)

    # Demonstrate trajectory-level geometric analysis:
    distances = traj.distance_trajectory(0, 1)
    print("\nDistance array (Atom 0 - Atom 1):")
    print(distances)

    # If there are at least 3 atoms, print an angle trajectory
    try:
        angles = traj.angle_trajectory(0, 1, 2)
        print("\nAngle array (Atoms 0, 1, 2):")
        print(angles)
    except Exception as e:
        print("\nCould not compute angle (0, 1, 2):", e)

    # If there are at least 4 atoms, print a dihedral trajectory
    try:
        dihedrals = traj.dihedral_trajectory(0, 1, 2, 3)
        print("\nDihedral array (Atoms 0, 1, 2, 3):")
        print(dihedrals)
    except Exception as e:
        print("\nCould not compute dihedral (0, 1, 2, 3):", e)

    # Filter demonstration (only show if multiple frames are present)
    if not traj.is_single_structure():
        print("\n--- Filter Demonstration ---")
        # Filter to frames with energy below the median (if energies are available)
        energies = [traj.get_frame(i).Energy for i in range(len(traj._raw_frames))]
        if all(e is not None for e in energies):
            median_energy = sorted(energies)[len(energies) // 2]
            filtered = traj.filter_frames(lambda mol: mol.Energy < median_energy)
            print(f"Filtered trajectory has {len(filtered._raw_frames)} frames (energy < {median_energy:.4f})")


if __name__ == "__main__":
    main()