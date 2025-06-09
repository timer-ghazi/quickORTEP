#!/usr/bin/env python3
"""
gaussian_parser.py
------------------
Parser for Gaussian log/output files.
Handles geometry extraction, energy parsing, and vibrational frequency analysis.
"""

import os
import re
from .base_parser import TrajectoryParser
from normal_modes import NormalMode, FrequencyData


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
        # Create an instance to access the frequency parsing method
        parser_instance = cls()
        
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
        freq_data = parser_instance._parse_frequencies(file_path)
        if freq_data is not None:
            # Attach the entire FrequencyData object to metadata
            metadata["frequency_data"] = freq_data
            # We assume freq data belongs to the final geometry frame (last in raw_frames)
            if raw_frames:
                metadata["freq_data_frame_index"] = len(raw_frames) - 1
        
        return raw_frames, metadata

    def _parse_frequencies(self, filename):
        """
        Parse vibrational frequencies from a Gaussian output file.
        
        Parameters:
            filename (str): Path to the Gaussian output file
            
        Returns:
            FrequencyData: Object containing all frequency data, or None if no data found
        """
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except IOError:
            return None

        # Find the number of atoms
        n_atoms = None
        for line in lines:
            match = re.search(r'NAtoms=\s+(\d+)', line)
            if match:
                n_atoms = int(match.group(1))
                break

        if n_atoms is None:
            return None

        freq_data = FrequencyData(n_atoms)

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Look for frequency blocks
            if line.startswith("Frequencies --"):
                try:
                    # Parse frequencies
                    freq_parts = line.split()
                    frequencies = [float(f) for f in freq_parts[2:]]  # Skip "Frequencies" and "--"
                    n_freqs = len(frequencies)
                    
                    if n_freqs == 0:
                        i += 1
                        continue

                    # Parse reduced masses
                    i += 1
                    if i >= len(lines):
                        break
                    red_mass_line = lines[i].strip()
                    if not "Red. masses" in red_mass_line:
                        i += 1
                        continue
                    red_mass_parts = red_mass_line.split()
                    # Find "--" and take everything after it
                    dash_idx = None
                    for idx, part in enumerate(red_mass_parts):
                        if "--" in part:
                            dash_idx = idx + 1
                            break
                    if dash_idx is None:
                        i += 1
                        continue
                    reduced_masses = [float(m) for m in red_mass_parts[dash_idx:]]

                    # Parse force constants
                    i += 1
                    if i >= len(lines):
                        break
                    frc_const_line = lines[i].strip()
                    if not "Frc consts" in frc_const_line:
                        i += 1
                        continue
                    frc_const_parts = frc_const_line.split()
                    # Find "--" and take everything after it
                    dash_idx = None
                    for idx, part in enumerate(frc_const_parts):
                        if "--" in part:
                            dash_idx = idx + 1
                            break
                    if dash_idx is None:
                        i += 1
                        continue
                    force_constants = [float(f) for f in frc_const_parts[dash_idx:]]

                    # Parse IR intensities
                    i += 1
                    if i >= len(lines):
                        break
                    ir_inten_line = lines[i].strip()
                    if not "IR Inten" in ir_inten_line:
                        i += 1
                        continue
                    ir_inten_parts = ir_inten_line.split()
                    # Find "--" and take everything after it
                    dash_idx = None
                    for idx, part in enumerate(ir_inten_parts):
                        if "--" in part:
                            dash_idx = idx + 1
                            break
                    if dash_idx is None:
                        i += 1
                        continue
                    ir_intensities = [float(ir) for ir in ir_inten_parts[dash_idx:]]

                    # Skip header lines for displacement vectors
                    i += 1
                    if i >= len(lines):
                        break

                    # Parse displacement vectors for each atom
                    displacements_by_freq = [[] for _ in range(n_freqs)]
                    
                    for atom_idx in range(n_atoms):
                        i += 1
                        if i >= len(lines):
                            break
                        disp_line = lines[i].strip()
                        parts = disp_line.split()
                        
                        if len(parts) < 2 + 3 * n_freqs:
                            continue
                            
                        # Extract displacement vectors (3 components per frequency)
                        for freq_idx in range(n_freqs):
                            start_idx = 2 + freq_idx * 3  # Skip atom number and symbol
                            if start_idx + 2 < len(parts):
                                dx = float(parts[start_idx])
                                dy = float(parts[start_idx + 1])
                                dz = float(parts[start_idx + 2])
                                displacements_by_freq[freq_idx].append((dx, dy, dz))

                    # Create NormalMode objects
                    for freq_idx in range(n_freqs):
                        if (freq_idx < len(frequencies) and 
                            freq_idx < len(reduced_masses) and 
                            freq_idx < len(force_constants) and 
                            freq_idx < len(ir_intensities)):
                            
                            mode = NormalMode(
                                frequency=frequencies[freq_idx],
                                reduced_mass=reduced_masses[freq_idx],
                                force_constant=force_constants[freq_idx],
                                ir_intensity=ir_intensities[freq_idx],
                                displacements=displacements_by_freq[freq_idx]
                            )
                            freq_data.add_mode(mode)

                except (ValueError, IndexError):
                    # Skip malformed frequency blocks
                    pass

            i += 1

        return freq_data if freq_data.modes else None

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