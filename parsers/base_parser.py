#!/usr/bin/env python3
"""
base_parser.py
--------------
Abstract base class for trajectory parsers.
Contains the TrajectoryParser interface that all parsers must implement.
"""

import os
from abc import ABC, abstractmethod


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