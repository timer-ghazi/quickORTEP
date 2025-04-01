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
import re
from abc import ABC, abstractmethod
from molecule_nci import MoleculeWithNCIs  # Assumes molecule_nci.py (and molecule.py) are in the path
from config import ENERGY_UNITS, DEFAULT_ENERGY_UNIT

# ADDED FOR FREQ DATA
# We'll import parse_gaussian_frequencies so we can attach freq data to metadata:
from normal_modes import parse_gaussian_frequencies


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
                    _ = int(first_line)
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


class GaussianTrajectoryParser(TrajectoryParser):
    """
    Parser for Gaussian log/output files containing geometries and energies.
    Handles both single-point calculations and optimization trajectories.
    """
    
    # List of orientation types in order of preference
    ORIENTATIONS = ["Standard", "Input", "Z-Matrix"]
    
    @classmethod
    def can_parse(cls, file_path):
        """
        Determines if this parser can handle the given file
        
        Parameters:
            file_path (str): Path to the file to check
            
        Returns:
            bool: True if this parser can handle the file
        """
        # Check file extension
        if not file_path.lower().endswith(('.log', '.out')):
            return False
            
        # Check content for Gaussian markers
        try:
            with open(file_path, 'r') as f:
                # Read first few lines to identify Gaussian output
                for _ in range(20):  # Check first 20 lines
                    line = next(f, '')
                    if "Gaussian" in line or "GAUSSIAN" in line:
                        return True
            return False
        except:
            return False
    
    @classmethod
    def parse(cls, file_path):
        """
        Parse Gaussian output file into raw frames and metadata
        
        Parameters:
            file_path (str): Path to the file to parse
            
        Returns:
            tuple: (raw_frames, metadata) where:
                  - raw_frames is a list of frame data blocks in XYZ format
                  - metadata is a dict with file-level metadata
        """
        # Extract base name for metadata
        base_name = os.path.basename(file_path)
        if base_name.lower().endswith(('.log', '.out')):
            base_name = base_name[:-4]
            
        metadata = {
            'file_name': base_name,
            'format': 'gaussian',
            'path': file_path,
            'energy_unit': 'hartree'  # Explicitly mark Gaussian energies as Hartrees
        }
        
        # Parse energies and geometries
        energies = cls._parse_opt_energies(file_path)
        geometries = cls._parse_geometries(file_path)
        
        # Convert to XYZ format frames
        raw_frames = []
        
        for step in sorted(energies.keys()):
            energy_value, energy_type = energies[step]
            
            # Get the preferred geometry for this step (0-indexed, so step-1)
            orientation, geom = cls._get_preferred_geometry(geometries, step-1)
            
            if geom:
                # Convert to XYZ format
                xyz_frame = cls._convert_to_xyz(geom, energy_value, energy_type, orientation, step)
                raw_frames.append(xyz_frame)
                
                # Store energy info in metadata for later reference
                if 'energy_data' not in metadata:
                    metadata['energy_data'] = {}
                metadata['energy_data'][len(raw_frames)-1] = {
                    'value': energy_value,
                    'type': energy_type,
                    'orientation': orientation,
                    'step': step
                }
        
        # ADDED FOR FREQ DATA:
        # Attempt to parse vibrational frequencies (if present).
        freq_data = parse_gaussian_frequencies(file_path)
        if freq_data is not None:
            # Attach the entire FrequencyData object to metadata
            metadata["frequency_data"] = freq_data
            # We assume freq data belongs to the final geometry frame (last in raw_frames)
            if raw_frames:
                metadata["freq_data_frame_index"] = len(raw_frames) - 1
        
        return raw_frames, metadata

    @staticmethod
    def _parse_opt_energies(logfile_path):
        """
        Parse geometry-optimization energies from a Gaussian log file.
        Returns a dictionary: step_number -> (final energy, label).

        The label might be 'SCF', 'MP2', 'Double-Hybrid', or 'External'.
        Overwrites SCF energies with correlated or external energies
        for the same step if present.
        """
        final_energies = {}    # step_number -> (energy_value, label)

        scf_done_pat = re.compile(
            r'^ *SCF Done:\s+E\(\w+\)\s*=\s*([-\d.]+(?:[DEe][+\-]?\d+)?)'
        )
        mp2_pat = re.compile(
            r'(?:EUMP2\s*=\s*|E\(MP2\)\s*=\s*)([-\d.]+(?:[DEe][+\-]?\d+)?)'
        )
        dh_correlation_pat = re.compile(r'\bE\([\w\d]+\)\s*=\s*([-\d.]+(?:[DEe][+\-]?\d+)?)')
        ext_pat = re.compile(r'^(?:Energy=|Recovered energy=)\s*([-\d.]+(?:[DEe][+\-]?\d+)?)')

        current_step = 0
        with open(logfile_path, 'r') as f:
            for raw_line in f:
                line = raw_line.strip()

                # SCF line => store (unless overridden)
                scf_match = scf_done_pat.match(line)
                if scf_match:
                    current_step += 1
                    val = float(scf_match.group(1).replace('D','E'))
                    final_energies[current_step] = (val, 'SCF')
                    continue

                # MP2 line => override if found
                mp2_match = mp2_pat.search(line)
                if mp2_match and current_step > 0:
                    val = float(mp2_match.group(1).replace('D','E'))
                    final_energies[current_step] = (val, 'MP2')
                    continue

                # Double-hybrid final correlation => E(...)= ...
                if 'E2(' in line and current_step > 0:
                    matches = dh_correlation_pat.findall(line)
                    if matches:
                        val_str = matches[-1]
                        val = float(val_str.replace('D','E'))
                        final_energies[current_step] = (val, 'Double-Hybrid')
                        continue

                # External lines => override
                ext_match = ext_pat.match(line)
                if ext_match:
                    current_step += 1
                    val = float(ext_match.group(1).replace('D','E'))
                    final_energies[current_step] = (val, 'External')
                    continue

        return final_energies

    @classmethod
    def _parse_geometries(cls, logfile_path):
        """
        Parse geometries from a Gaussian log file.
        Returns a dictionary: orientation_type -> list of geometry blocks
        """
        data = {key: [] for key in cls.ORIENTATIONS}
        
        with open(logfile_path, 'r') as file:
            lines = iter(file)
            for line in lines:
                words = line.split()
                if len(words) >= 2 and words[0] in cls.ORIENTATIONS and words[1] == "orientation:":
                    orientation_type = words[0]

                    # Skip header lines
                    for _ in range(4):
                        next(lines, None)

                    collected_lines = []
                    for data_line in lines:
                        entries = data_line.split()
                        if len(entries) == 1:
                            break
                        if len(entries) >= 6:
                            atom_num, x, y, z = entries[1], entries[3], entries[4], entries[5]
                            collected_lines.append((int(atom_num), float(x), float(y), float(z)))

                    data[orientation_type].append(collected_lines)

        return data

    @classmethod
    def _get_preferred_geometry(cls, geometries, index):
        """
        Get the preferred geometry from the available orientations.
        Returns (orientation, geometry) or (None, None) if no geometry found
        """
        for orientation in cls.ORIENTATIONS:
            if geometries[orientation] and index < len(geometries[orientation]):
                return orientation, geometries[orientation][index]
        return None, None

    @staticmethod
    def _convert_to_xyz(geometry, energy_value, energy_type, orientation, step):
        """
        Convert a Gaussian geometry block to XYZ format
        
        Parameters:
            geometry: List of tuples (atomic_number, x, y, z)
            energy_value: The energy value
            energy_type: The energy type (SCF, MP2, etc.)
            orientation: The orientation type (Standard, Input, Z-Matrix)
            step: The calculation step number
            
        Returns:
            List of strings in XYZ format
        """
        from elements_table import Elements  # Import here to avoid circular imports
        
        num_atoms = len(geometry)
        
        xyz_frame = [
            str(num_atoms),
            f"Step {step} | {orientation} orientation | Energy= {energy_value} | Type: {energy_type}"
        ]
        
        for atom_data in geometry:
            atomic_num, x, y, z = atom_data
            symbol = Elements.symbol(atomic_num)
            xyz_frame.append(f"{symbol:<5}{x:17.12f}{y:17.12f}{z:17.12f}")
        
        return xyz_frame


class Trajectory:
    """
    Class for handling molecular trajectories from various file formats.
    Uses lazy evaluation with caching for memory efficiency.
    """
    
    # Register all available parsers
    _parsers = [XYZTrajectoryParser, GaussianTrajectoryParser]
    
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
        
    def update_frame(self, frame_index, molecule):
        """
        Update a specific frame with a modified molecule.
        
        Parameters:
            frame_index (int): Index of the frame to update
            molecule: The molecule object with updated data
            
        Raises:
            ValueError: If frame_index is out of range
        """
        if frame_index < 0 or frame_index >= len(self._raw_frames):
            raise ValueError(f"Invalid frame index: {frame_index}")
            
        # Update the molecule cache
        self._molecule_cache[frame_index] = molecule
        
        # Update the raw frame data
        xyz_data = molecule.to_xyz().splitlines()
        self._raw_frames[frame_index] = xyz_data
        
        # Update energy if available
        self._frame_energies[frame_index] = molecule.Energy

    def _detect_hartree_units(self):
        """
        Detect if energies are in Hartrees based on magnitude and source.
        
        Returns:
            bool: True if energies are likely to be in Hartrees
        """
        # Check if we know the format from metadata
        if 'format' in self.metadata:
            # Gaussian calculations are typically in Hartrees
            if self.metadata['format'] == 'gaussian':
                return True
        
        # Check if energy_unit is explicitly set in metadata
        if 'energy_unit' in self.metadata:
            return self.metadata['energy_unit'].lower() == 'hartree'
        
        # Check if there are energy values with explicit unit info
        if 'energy_data' in self.metadata:
            for frame_data in self.metadata['energy_data'].values():
                if 'type' in frame_data:
                    # Gaussian energy types (SCF, MP2, etc.) are in Hartrees
                    if frame_data['type'] in ['SCF', 'MP2', 'Double-Hybrid', 'External']:
                        return True
        
        # Check the magnitude of energies
        energies = np.array([e for e in self._frame_energies if e is not None])
        if len(energies) > 0:
            # Hartree energies for molecular systems are typically negative and of order 10^0 to 10^3
            mean_energy = np.mean(energies)
            if -10000 < mean_energy < 0:
                return True
        
        return False

    def energy_trajectory(self, skip_none=False, convert_if_hartrees=True, convert_to_unit=None):
        """
        Return a NumPy array of energies across all frames with optional unit conversion.
        
        Parameters:
            skip_none (bool): If True, frames with None energies will be excluded.
            convert_if_hartrees (bool): If True, will convert energies if detected as Hartrees.
            convert_to_unit (str): Target unit for conversion (hartree, kcal_mol, kj_mol, ev).
                          If None, uses DEFAULT_ENERGY_UNIT from config.
        
        Returns:
            tuple: (np.ndarray, dict) where:
                - np.ndarray is the array of energy values
                - dict contains metadata about the conversion
        """
        # Force loading of all frames to ensure energies are populated
        for i in range(len(self._raw_frames)):
            self.get_frame(i)
        
        # Get raw energies
        if skip_none:
            energies = np.array([e for e in self._frame_energies if e is not None])
        else:
            energies = np.array([e if e is not None else np.nan for e in self._frame_energies])
        
        if len(energies) == 0 or np.isnan(energies).all():
            return energies, {'original_unit': 'unknown', 'converted_unit': 'unknown', 'normalized': False}
        
        conversion_info = {
            'original_unit': 'unknown',
            'converted_unit': 'unknown',
            'normalized': False,
            'min_energy': None,
            'conversion_factor': 1.0
        }
        
        is_hartrees = self._detect_hartree_units()
        
        if is_hartrees:
            conversion_info['original_unit'] = 'hartree'
        elif 'energy_unit' in self.metadata:
            conversion_info['original_unit'] = self.metadata['energy_unit']
        
        valid_energies = energies[~np.isnan(energies)]
        min_energy = np.min(valid_energies)
        conversion_info['min_energy'] = min_energy
        
        if convert_if_hartrees and is_hartrees:
            target_unit = convert_to_unit if convert_to_unit else DEFAULT_ENERGY_UNIT
            if target_unit in ENERGY_UNITS:
                conversion_factor = ENERGY_UNITS['hartree'][f'to_{target_unit}']
                mask = ~np.isnan(energies)
                energies_converted = energies.copy()
                energies_converted[mask] = energies[mask] * conversion_factor
                min_energy_converted = min_energy * conversion_factor
                energies_normalized = energies_converted.copy()
                energies_normalized[mask] = energies_converted[mask] - min_energy_converted
                
                conversion_info['converted_unit'] = target_unit
                conversion_info['conversion_factor'] = conversion_factor
                conversion_info['normalized'] = True
                
                return energies_normalized, conversion_info
        
        conversion_info['converted_unit'] = conversion_info['original_unit']
        return energies, conversion_info

    def coordinate_trajectory(self, skip_none=False):
        """
        Return a NumPy array of coordinate (or time) values across all frames,
        extracted from the molecule metadata.
        
        If the metadata field 'coord' or 'time' is present in the molecule metadata,
        that value will be used for the x-axis in plots. If both are present, 'coord'
        is preferred. If none is found, frame indices are returned.
        
        Parameters:
            skip_none (bool): If True, frames with None values will be excluded.
        
        Returns:
            tuple: (np.ndarray, dict) where:
                - np.ndarray is the array of coordinate values
                - dict contains metadata about the extracted field (e.g., 'field': 'coord')
        """
        # Force loading of all frames to ensure metadata is populated
        for i in range(len(self._raw_frames)):
            self.get_frame(i)
        
        n = len(self._raw_frames)
        
        # Determine which metadata field to use
        field = None
        for i in range(n):
            mol = self.get_frame(i)
            if 'coord' in mol.metadata:
                field = 'coord'
                break
            elif 'time' in mol.metadata:
                field = 'time'
                break
        
        if field is None:
            # No coordinate/time metadata found; return frame indices
            coords = np.array(list(range(n)))
            conversion_info = {'field': None}
            return coords, conversion_info
        
        # Gather coordinate values from each frame
        coords = []
        for i in range(n):
            mol = self.get_frame(i)
            value = mol.metadata.get(field, np.nan)
            coords.append(value)
        coords = np.array(coords)
        
        if skip_none:
            coords = coords[~np.isnan(coords)]
        
        conversion_info = {
            'field': field,
            'min_value': np.nanmin(coords) if not np.isnan(coords).all() else None
        }
        return coords, conversion_info

    def print_energy_table(self, convert_if_hartrees=True, convert_to_unit=None):
        """
        Print a table of frame numbers vs. energy.
        For each frame, the molecule is fully processed (if not already) so that
        its Energy attribute is available.
        """
        energies, energy_info = self.energy_trajectory(
            convert_if_hartrees=convert_if_hartrees,
            convert_to_unit=convert_to_unit
        )
        
        unit_display = ENERGY_UNITS.get(
            energy_info['converted_unit'], 
            {'symbol': energy_info['converted_unit']}
        )['symbol']
        
        normalized_str = " (rel. to min)" if energy_info['normalized'] else ""
        
        print(f"Frame\tEnergy ({unit_display}{normalized_str})\tType\tOrientation")
        print("----------------------------------------")
        
        for i in range(len(self._raw_frames)):
            energy_val = energies[i]
            energy_str = f"{energy_val:.6f}" if not np.isnan(energy_val) else "N/A"
            
            energy_type = "N/A"
            orientation = "N/A"
            if self.metadata and 'energy_data' in self.metadata and i in self.metadata['energy_data']:
                energy_type = self.metadata['energy_data'][i].get('type', 'N/A')
                orientation = self.metadata['energy_data'][i].get('orientation', 'N/A')
                
            print(f"{i:3d}\t{energy_str}\t{energy_type}\t{orientation}")

    def print_distance_table(self, atom1, atom2):
        """
        Print a table of the distance between two specified atoms for each frame.
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

    def summary(self, convert_if_hartrees=True, convert_to_unit=None):
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
                if key != 'energy_data':  # Skip detailed energy data in summary
                    print(f"  {key}: {value}")
            
        # Get energies with unit conversion
        energies, energy_info = self.energy_trajectory(
            convert_if_hartrees=convert_if_hartrees,
            convert_to_unit=convert_to_unit
        )
        
        valid_energies = energies[~np.isnan(energies)]
        if len(valid_energies) > 0:
            unit_display = ENERGY_UNITS.get(
                energy_info['converted_unit'], 
                {'symbol': energy_info['converted_unit']}
            )['symbol']
            
            normalized_str = " (relative to minimum)" if energy_info['normalized'] else ""
            
            print(f"\nEnergy range: {np.min(valid_energies):.6f} to {np.max(valid_energies):.6f} {unit_display}{normalized_str}")
            
            if energy_info['original_unit'] != energy_info['converted_unit'] and energy_info['original_unit'] != 'unknown':
                original_unit = ENERGY_UNITS.get(
                    energy_info['original_unit'], 
                    {'name': energy_info['original_unit']}
                )['name']
                
                converted_unit = ENERGY_UNITS.get(
                    energy_info['converted_unit'], 
                    {'name': energy_info['converted_unit']}
                )['name']
                
                print(f"Energy conversion: {original_unit} → {converted_unit}")
            
            if 'format' in self.metadata and self.metadata['format'] == 'gaussian':
                energy_types = set()
                if 'energy_data' in self.metadata:
                    for data in self.metadata['energy_data'].values():
                        if 'type' in data:
                            energy_types.add(data['type'])
                if energy_types:
                    print(f"Energy types: {', '.join(energy_types)}")
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

    # Print overall summary with unit conversion
    print("\n--- Trajectory Summary ---")
    traj.summary(convert_if_hartrees=True, convert_to_unit=DEFAULT_ENERGY_UNIT)

    # Print an energy vs. frame table with unit conversion
    print("\n--- Energy vs. Frame ---")
    traj.print_energy_table(convert_if_hartrees=True, convert_to_unit=DEFAULT_ENERGY_UNIT)

    # Print a distance vs. frame table for atoms 0 and 1
    print("\n--- Distance vs. Frame (Atom 0 - Atom 1) ---")
    traj.print_distance_table(atom1=0, atom2=1)

    # Demonstrate trajectory-level geometric analysis:
    distances = traj.distance_trajectory(0, 1)
    print("\nDistance array (Atom 0 - Atom 1):")
    print(distances)

    # Get energy information with unit conversion
    energies, energy_info = traj.energy_trajectory(
        convert_if_hartrees=True, 
        convert_to_unit=DEFAULT_ENERGY_UNIT
    )
    print("\nEnergy array:")
    print(energies)
    print(f"Units: {ENERGY_UNITS[energy_info['converted_unit']]['name']}")
    if energy_info['normalized']:
        print("Energies normalized relative to minimum")

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
        energies, _ = traj.energy_trajectory()
        valid_energies = energies[~np.isnan(energies)]
        if len(valid_energies) > 0:
            median_energy = np.median(valid_energies)
            filtered = traj.filter_frames(lambda mol: mol.Energy is not None and mol.Energy < median_energy)
            print(f"Filtered trajectory has {len(filtered._raw_frames)} frames (energy < {median_energy:.6f})")


if __name__ == "__main__":
    main()
