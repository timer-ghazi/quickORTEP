#!/usr/bin/env python3
"""
Parse vibrational frequencies from a Gaussian frequency log file.
Stores them in a FrequencyData object, which has NormalMode entries
for each vibrational mode. Prints out normal mode displacements for
modes with imaginary frequencies so you can debug/inspect them.
"""

import re

class NormalMode:
    """
    Stores one vibrational mode from a Gaussian frequency calculation.
    """
    def __init__(self, frequency, red_mass, force_const, ir_intensity, displacements):
        """
        :param frequency: float (cm^-1). Negative if imaginary.
        :param red_mass: float (amu).
        :param force_const: float (mDyne/A).
        :param ir_intensity: float (KM/mol).
        :param displacements: list of (dx, dy, dz) for each atom.
        """
        self.frequency = frequency
        self.red_mass = red_mass
        self.force_const = force_const
        self.ir_intensity = ir_intensity
        self.displacements = displacements  # list of tuples

    def __repr__(self):
        return (f"<NormalMode freq={self.frequency:.2f} cm^-1, "
                f"IR={self.ir_intensity:.2f} KM/mol>")


class FrequencyData:
    """
    Container holding all normal modes for a single Gaussian freq calculation.
    """
    def __init__(self, n_atoms=0):
        self.n_atoms = n_atoms
        self.modes = []

    def add_mode(self, mode: NormalMode):
        self.modes.append(mode)

    def __repr__(self):
        return (f"<FrequencyData with {len(self.modes)} modes; "
                f"{self.n_atoms} atoms>")


def parse_gaussian_frequencies(filename):
    """
    Read a Gaussian frequency log file and return a FrequencyData object.

    :param filename: Path to the Gaussian .log or .out file.
    :return: FrequencyData object with all the normal modes parsed.
    """
    freq_data = FrequencyData()

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Try to locate the number of atoms from lines like "NAtoms=      7".
    n_atoms_pattern = re.compile(r'NAtoms=\s+(\d+)')
    for line in lines:
        match = n_atoms_pattern.search(line)
        if match:
            freq_data.n_atoms = int(match.group(1))
            break

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        line = lines[i].strip()

        # Look for a line that starts with "Frequencies --"
        if line.startswith("Frequencies --"):
            # 1) Parse frequencies
            freq_tokens = line.split()[2:]  # skip "Frequencies --"
            freqs = [float(x) for x in freq_tokens]

            # 2) Next line: "Red. masses --"
            i += 1
            red_line = lines[i].strip()
            red_tokens = red_line.split()[3:]  # skip "Red.", "masses", "--"
            red_masses = [float(x) for x in red_tokens]

            # 3) Next line: "Frc consts  --"
            i += 1
            frc_line = lines[i].strip()
            frc_tokens = frc_line.split()[3:]  # skip "Frc", "consts", "--"
            frc_consts = [float(x) for x in frc_tokens]

            # 4) Next line: "IR Inten    --"
            i += 1
            ir_line = lines[i].strip()
            ir_tokens = ir_line.split()[3:]  # skip "IR", "Inten", "--"
            ir_ints = [float(x) for x in ir_tokens]

            # 5) Next line should be the header " Atom  AN   X   Y   Z   X   Y   Z ..."
            i += 1
            # Now parse displacements (we need freq_data.n_atoms lines).
            n_atoms = freq_data.n_atoms
            if n_atoms == 0:
                raise ValueError("NAtoms not found in file. Please set it manually.")

            # Prepare a list of displacement-lists for each freq column.
            freq_block_displacements = [[] for _ in freqs]

            # For each of the n_atoms lines, read the displacement data:
            for _ in range(n_atoms):
                i += 1
                disp_line = lines[i].strip().split()
                # skip the first two entries (atom index, atomic number)
                for j in range(len(freqs)):
                    col_start = 2 + j*3
                    dx = float(disp_line[col_start])
                    dy = float(disp_line[col_start+1])
                    dz = float(disp_line[col_start+2])
                    freq_block_displacements[j].append((dx, dy, dz))

            # Create NormalMode objects for each frequency in this block
            for j in range(len(freqs)):
                mode = NormalMode(
                    frequency=freqs[j],
                    red_mass=red_masses[j],
                    force_const=frc_consts[j],
                    ir_intensity=ir_ints[j],
                    displacements=freq_block_displacements[j]
                )
                freq_data.add_mode(mode)
        else:
            i += 1

    return freq_data


# -------------- Example usage (uncomment to run) --------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python parse_freq.py <freq.log>")
        sys.exit(1)

    filename = sys.argv[1]
    data = parse_gaussian_frequencies(filename)
    print(data)

    # Print first few modes to confirm parsing
    print("First few modes:")
    for idx, mode in enumerate(data.modes[:5]):
        print(f"Mode {idx+1}: freq={mode.frequency:.2f} cm^-1, IR={mode.ir_intensity:.2f} KM/mol")

    print("\nChecking for imaginary modes...")
    for idx, mode in enumerate(data.modes, start=1):
        if mode.frequency < 0.0:
            print(f"\nImaginary Mode {idx}: freq={mode.frequency:.2f} cm^-1")
            print("Displacements (x, y, z) for each atom:")
            for atom_i, disp in enumerate(mode.displacements, start=1):
                print(f"  Atom {atom_i}: {disp}")
