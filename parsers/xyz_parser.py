#!/usr/bin/env python3
"""
xyz_parser.py
-------------
Parser for XYZ trajectory files.
Handles concatenated XYZ files with multiple frames.
"""

import os
from .base_parser import TrajectoryParser


class XYZTrajectoryParser(TrajectoryParser):
    """Parser for concatenated XYZ files"""
    
    @classmethod
    def can_parse(cls, file_path):
        """Check if file is in XYZ format"""
        # Simple check: file extension and basic structure
        allowed_extensions = ('.xyz', '.trj')
        if not any(file_path.lower().endswith(ext) for ext in allowed_extensions):
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