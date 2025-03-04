#!/usr/bin/env python3
"""
Parse vibrational frequencies from a Gaussian .log/.out file.
If no frequencies are found, returns None.
Otherwise, returns a FrequencyData object containing NormalMode objects.
Prints out details for imaginary frequencies for debugging.
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
        self.displacements = displacements  # list of (dx, dy, dz)

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
    Read a Gaussian output file and return:
    - FrequencyData object if vibrational data is found,
    - None otherwise (indicating no frequency data present).
    
    :param filename: Path to the Gaussian .log or .out file.
    :return: FrequencyData object or None.
    """
    freq_data = FrequencyData()
    found_any_frequency = False

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Attempt to find NAtoms from lines like "NAtoms=  7"
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
            found_any_frequency = True

            # 1) Parse frequencies
            freq_tokens = line.split()[2:]  # skip ["Frequencies", "--"]
            freqs = [float(x) for x in freq_tokens]

            # 2) Next line: "Red. masses --"
            i += 1
            red_line = lines[i].strip()
            red_tokens = red_line.split()[3:]  # skip ["Red.", "masses", "--"]
            red_masses = [float(x) for x in red_tokens]

            # 3) Next line: "Frc consts  --"
            i += 1
            frc_line = lines[i].strip()
            frc_tokens = frc_line.split()[3:]  # skip ["Frc", "consts", "--"]
            frc_consts = [float(x) for x in frc_tokens]

            # 4) Next line: "IR Inten    --"
            i += 1
            ir_line = lines[i].strip()
            ir_tokens = ir_line.split()[3:]   # skip ["IR", "Inten", "--"]
            ir_ints = [float(x) for x in ir_tokens]

            # 5) Next line is the header " Atom  AN   X   Y   Z   X   Y   Z..."
            i += 1

            # We should have freq_data.n_atoms lines of displacements
            n_atoms = freq_data.n_atoms
            if n_atoms == 0:
                raise ValueError("NAtoms not found in file. Please set it manually.")

            # Prepare containers for each frequency block
            freq_block_displacements = [[] for _ in freqs]

            for _ in range(n_atoms):
                i += 1
                disp_line = lines[i].strip().split()
                # Format: [AtomIndex, AtomicNum, dx(freq1), dy(freq1), dz(freq1), dx(freq2), ... ]
                for j in range(len(freqs)):
                    col_start = 2 + j * 3
                    dx = float(disp_line[col_start])
                    dy = float(disp_line[col_start + 1])
                    dz = float(disp_line[col_start + 2])
                    freq_block_displacements[j].append((dx, dy, dz))

            # Create NormalMode objects for each freq
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

    # If no frequency sections found, return None.
    if not found_any_frequency:
        return None

    return freq_data


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python parse_freq.py <freq.log>")
        sys.exit(1)

    filename = sys.argv[1]
    data = parse_gaussian_frequencies(filename)

    if data is None:
        print("No frequencies found in Gaussian output.")
        sys.exit(0)

    print(data)

    # Print some details:
    print("\nFirst few modes:")
    for idx, mode in enumerate(data.modes[:5], start=1):
        print(f"Mode {idx}: freq={mode.frequency:.2f} cm^-1,"
              f" IR={mode.ir_intensity:.2f} KM/mol")

    # Print imaginary frequencies with their displacements:
    print("\nChecking for imaginary modes...")
    for idx, mode in enumerate(data.modes, start=1):
        if mode.frequency < 0.0:
            print(f"\nImaginary Mode {idx}: freq={mode.frequency:.2f} cm^-1")
            print("Displacements (x, y, z) for each atom:")
            for atom_i, disp in enumerate(mode.displacements, start=1):
                print(f"  Atom {atom_i}: {disp}")
