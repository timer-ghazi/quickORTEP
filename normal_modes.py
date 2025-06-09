#!/usr/bin/env python3
"""
normal_modes.py
---------------
Classes and utilities for handling vibrational normal modes and frequency data.
Contains general data structures for storing vibrational analysis results.
"""


class NormalMode:
    """
    Represents a single vibrational normal mode.
    Contains frequency, reduced mass, force constant, IR intensity,
    and displacement vectors for all atoms.
    """
    def __init__(self, frequency, reduced_mass, force_constant, ir_intensity, displacements):
        self.frequency = frequency  # in cm^-1
        self.reduced_mass = reduced_mass  # in amu
        self.force_constant = force_constant  # in mdyn/A
        self.ir_intensity = ir_intensity  # in KM/Mole
        self.displacements = displacements  # List of (dx, dy, dz) tuples for each atom


class FrequencyData:
    """
    Container for vibrational frequency data from a quantum chemistry calculation.
    Contains the number of atoms and a list of normal modes.
    """
    def __init__(self, n_atoms):
        self.n_atoms = n_atoms
        self.modes = []  # List of NormalMode objects
    
    def add_mode(self, mode):
        """Add a NormalMode to this frequency data."""
        self.modes.append(mode)