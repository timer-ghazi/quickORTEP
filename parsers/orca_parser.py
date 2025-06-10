#!/usr/bin/env python3
"""
orca_parser.py
--------------
Parser for ORCA electronic structure calculation output files.
Handles geometry optimization trajectories, energies, and vibrational frequency analysis.
"""

import os
import re
from .base_parser import TrajectoryParser
from normal_modes import NormalMode, FrequencyData


class ORCATrajectoryParser(TrajectoryParser):
    """
    Parser for ORCA output files containing geometries and energies.
    Handles both single-point calculations and optimization trajectories.
    """
    
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
        if not file_path.lower().endswith('.out'):
            return False
            
        # Check content for ORCA markers
        try:
            with open(file_path, 'r') as f:
                # Read first few lines to identify ORCA output
                for _ in range(50):  # Check first 50 lines for ORCA markers
                    line = next(f, '')
                    if not line:
                        break
                    # Look for distinctive ORCA markers
                    if any(marker in line for marker in [
                        "* O   R   C   A *",
                        "ORCA TERMINATED NORMALLY",
                        "This ORCA versions",
                        "Program Version"
                    ]):
                        return True
            return False
        except:
            return False
    
    @classmethod
    def parse(cls, file_path):
        """
        Parse ORCA output file into raw frames and metadata
        
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
        if base_name.lower().endswith('.out'):
            base_name = base_name[:-4]
            
        metadata = {
            'file_name': base_name,
            'format': 'orca',
            'path': file_path,
            'energy_unit': 'hartree'  # ORCA energies are in Hartrees
        }
        
        # Parse energies and geometries
        energies = cls._parse_energies(file_path)
        geometries = cls._parse_geometries(file_path)
        
        # Check if optimization converged
        converged = cls._check_convergence(file_path)
        metadata['converged'] = converged
        
        # Convert to XYZ format frames
        raw_frames = []
        
        for cycle in sorted(energies.keys()):
            energy_value = energies[cycle]
            
            # Get geometry for this cycle (cycles are 1-indexed, list is 0-indexed)
            if cycle <= len(geometries):
                geometry = geometries[cycle - 1]
                
                if geometry:
                    # Convert to XYZ format
                    xyz_frame = cls._convert_to_xyz(geometry, energy_value, cycle)
                    raw_frames.append(xyz_frame)
                    
                    # Store energy info in metadata for later reference
                    if 'energy_data' not in metadata:
                        metadata['energy_data'] = {}
                    metadata['energy_data'][len(raw_frames)-1] = {
                        'value': energy_value,
                        'type': 'DFT',
                        'cycle': cycle
                    }
        
        # Parse vibrational frequencies (if present)
        freq_data = parser_instance._parse_frequencies(file_path)
        if freq_data is not None:
            # Attach the entire FrequencyData object to metadata
            metadata["frequency_data"] = freq_data
            # Assume freq data belongs to the final geometry frame (last in raw_frames)
            if raw_frames:
                metadata["freq_data_frame_index"] = len(raw_frames) - 1
        
        return raw_frames, metadata
    
    @classmethod
    def _parse_energies(cls, file_path):
        """
        Parse optimization cycle energies from ORCA output file.
        Returns a dictionary: cycle_number -> energy_value
        """
        energies = {}
        
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Pattern to match FINAL SINGLE POINT ENERGY lines
        energy_pattern = re.compile(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)')
        
        # Find all energy matches
        energy_matches = energy_pattern.findall(content)
        
        # Assign cycle numbers (1-indexed)
        for i, energy_str in enumerate(energy_matches):
            cycle = i + 1
            energies[cycle] = float(energy_str)
        
        return energies
    
    @classmethod
    def _parse_geometries(cls, file_path):
        """
        Parse geometries from ORCA optimization cycles.
        Returns a list of geometry blocks (one per cycle)
        """
        geometries = []
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Look for geometry optimization cycle headers
            if "GEOMETRY OPTIMIZATION CYCLE" in line:
                # Skip to the cartesian coordinates section
                while i < len(lines) and "CARTESIAN COORDINATES (ANGSTROEM)" not in lines[i]:
                    i += 1
                
                if i >= len(lines):
                    break
                
                # Skip the header line
                i += 1
                
                # Skip dashed line
                if i < len(lines) and "---" in lines[i]:
                    i += 1
                
                # Parse atom coordinates
                geometry = []
                while i < len(lines):
                    coord_line = lines[i].strip()
                    if not coord_line or "---" in coord_line or coord_line.startswith("*"):
                        break
                    
                    parts = coord_line.split()
                    if len(parts) >= 4:
                        try:
                            symbol = parts[0]
                            x = float(parts[1])
                            y = float(parts[2]) 
                            z = float(parts[3])
                            geometry.append((symbol, x, y, z))
                        except (ValueError, IndexError):
                            break
                    else:
                        break
                    i += 1
                
                if geometry:
                    geometries.append(geometry)
            
            i += 1
        
        return geometries
    
    @staticmethod
    def _convert_to_xyz(geometry, energy_value, cycle):
        """
        Convert an ORCA geometry block to XYZ format
        
        Parameters:
            geometry: List of tuples (symbol, x, y, z)
            energy_value: The energy value in Hartree
            cycle: The optimization cycle number
            
        Returns:
            List of strings in XYZ format
        """
        num_atoms = len(geometry)
        
        xyz_frame = [
            str(num_atoms),
            f"Cycle {cycle} | Energy= {energy_value} | Type: DFT"
        ]
        
        for atom_data in geometry:
            symbol, x, y, z = atom_data
            xyz_frame.append(f"{symbol:<5}{x:17.12f}{y:17.12f}{z:17.12f}")
        
        return xyz_frame
    
    @classmethod
    def _check_convergence(cls, file_path):
        """
        Check if the optimization converged
        
        Parameters:
            file_path (str): Path to the ORCA output file
            
        Returns:
            bool: True if optimization converged
        """
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Look for convergence indicators
            convergence_patterns = [
                "THE OPTIMIZATION HAS CONVERGED",
                "OPTIMIZATION CONVERGED",
                "HURRAY"
            ]
            
            for pattern in convergence_patterns:
                if pattern in content:
                    return True
            
            return False
        except:
            return False
    
    def _parse_frequencies(self, file_path):
        """
        Parse vibrational frequencies from an ORCA output file.
        
        Parameters:
            file_path (str): Path to the ORCA output file
            
        Returns:
            FrequencyData: Object containing all frequency data, or None if no data found
        """
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except IOError:
            return None
        
        # Find the number of atoms from the final geometry
        n_atoms = 7  # Default for our test case, but let's try to detect it
        
        # Try to determine number of atoms from geometry sections
        for i, line in enumerate(lines):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                # Skip header and dashed line
                j = i + 2
                atom_count = 0
                while j < len(lines):
                    coord_line = lines[j].strip()
                    if not coord_line or "---" in coord_line or coord_line.startswith("*"):
                        break
                    parts = coord_line.split()
                    if len(parts) >= 4:
                        atom_count += 1
                    j += 1
                if atom_count > 0:
                    n_atoms = atom_count
                    break
        
        freq_data = FrequencyData(n_atoms)
        
        # Find vibrational frequencies section
        freq_start = -1
        for i, line in enumerate(lines):
            if "VIBRATIONAL FREQUENCIES" in line:
                freq_start = i
                break
        
        if freq_start == -1:
            return None
        
        # Parse frequencies
        frequencies = []
        i = freq_start + 1
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("NORMAL MODES"):
                break
            
            # Parse frequency lines: "     6:     127.27 cm**-1"
            if ":" in line and "cm**-1" in line:
                parts = line.split(":")
                if len(parts) >= 2:
                    freq_part = parts[1].strip().replace("cm**-1", "").strip()
                    try:
                        freq_value = float(freq_part)
                        frequencies.append(freq_value)
                    except ValueError:
                        pass
            i += 1
        
        # Find normal modes section
        modes_start = -1
        for i, line in enumerate(lines):
            if "NORMAL MODES" in line:
                modes_start = i
                break
        
        if modes_start == -1:
            return freq_data if frequencies else None
        
        # Parse normal mode displacement vectors
        i = modes_start + 1
        
        # Skip description lines until we find the mode header
        while i < len(lines):
            line = lines[i].strip()
            # Look for mode header like "                  0          1          2          3          4          5"
            if re.match(r'\s*\d+\s+\d+\s+\d+', line):
                break
            i += 1
        
        # Parse mode blocks
        while i < len(lines):
            line = lines[i].strip()
            
            # Check for mode header line (mode numbers with lots of spacing)
            if re.match(r'\s*\d+\s+\d+\s+\d+', line):
                # Parse mode numbers from header
                mode_numbers = [int(x) for x in line.split()]
                n_modes = len(mode_numbers)
                
                if n_modes == 0:
                    i += 1
                    continue
                
                # Initialize displacement arrays for each mode (N atoms × 3 coordinates)
                mode_displacements = [[[0.0, 0.0, 0.0] for _ in range(n_atoms)] for _ in range(n_modes)]
                
                # Parse displacement vectors (3N lines for N atoms)
                i += 1
                coord_count = 0
                
                while i < len(lines) and coord_count < 3 * n_atoms:
                    disp_line = lines[i].strip()
                    
                    # Stop if we hit another mode block or end
                    if re.match(r'\s*\d+\s+\d+\s+\d+', disp_line) or not disp_line:
                        break
                    
                    parts = disp_line.split()
                    if len(parts) >= n_modes + 1:  # coord_idx + mode values
                        try:
                            coord_idx = int(parts[0])  # First column is coordinate index
                            atom_idx = coord_idx // 3  # Which atom (0, 1, 2, ...)
                            coord_component = coord_idx % 3  # Which coordinate (0=x, 1=y, 2=z)
                            
                            # Extract displacement values for each mode
                            for mode_idx in range(n_modes):
                                displacement_val = float(parts[mode_idx + 1])
                                
                                if atom_idx < n_atoms:
                                    mode_displacements[mode_idx][atom_idx][coord_component] = displacement_val
                            
                            coord_count += 1
                        except (ValueError, IndexError):
                            pass
                    
                    i += 1
                
                # Create NormalMode objects for each mode in this block
                for mode_idx, mode_num in enumerate(mode_numbers):
                    if mode_num < len(frequencies):
                        frequency = frequencies[mode_num]
                        
                        # Convert displacement format: [[x,y,z], [x,y,z], ...] -> [(x,y,z), (x,y,z), ...]
                        displacements = []
                        if mode_idx < len(mode_displacements):
                            for atom_disp in mode_displacements[mode_idx]:
                                displacements.append((atom_disp[0], atom_disp[1], atom_disp[2]))
                        
                        # ORCA doesn't provide reduced mass, force constant, or IR intensity directly
                        # Set default values (could be enhanced to parse these if available)
                        mode = NormalMode(
                            frequency=frequency,
                            reduced_mass=1.0,  # Default value
                            force_constant=0.0,  # Default value
                            ir_intensity=0.0,  # Default value
                            displacements=displacements
                        )
                        freq_data.add_mode(mode)
                
                continue
            
            i += 1
        
        return freq_data if freq_data.modes else None