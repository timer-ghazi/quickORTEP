# quickORTEP

A lightweight, fast molecular visualization tool with interactive bond editing and trajectory analysis capabilities.

![quickORTEP Screenshot](screenshot.png)

## Overview

quickORTEP is a Python-based molecular visualization tool inspired by the classic ORTEP (Oak Ridge Thermal Ellipsoid Plot) program. It provides an interactive, keyboard-driven interface for viewing molecular structures and trajectories, with special attention to bond visualization and manipulation.

## Features

- **Fast X11-based rendering** with optional anti-aliasing (supersampling)
- **Trajectory support** for visualizing molecular dynamics or reaction paths
- **Interactive bond editing** that propagates across trajectory frames
- **Multiple bond types** (covalent, non-covalent interactions, transition state bonds, distance bonds)
- **Real-time energy and bond length plots** for trajectory analysis
- **SVG export** for high-quality publication graphics
- **Keyboard-driven interface** for efficient workflow

## Installation

### Dependencies

- Python 3.6+
- python-xlib
- numpy
- PIL (Pillow)

### Installing from Source

```bash
# Clone the repository
git clone https://github.com/yourusername/quickORTEP.git
cd quickORTEP

# Install dependencies
pip install python-xlib numpy pillow

# Make the main script executable
chmod +x main.py
```

## Basic Usage

### Starting the Viewer

```bash
# Basic usage
python main.py molecule.xyz

# With anti-aliasing (2x supersampling)
python main.py molecule.xyz 2

# With custom tile size for larger displays
python main.py molecule.xyz 2 256
```

### Command Line Arguments

- **Argument 1**: Path to XYZ file (required)
- **Argument 2**: Supersampling factor for anti-aliasing (optional, default=1)
- **Argument 3**: Tile size for rendering (optional, default=128)

## User Interface

### Mouse Controls

- **Left click**: Select atoms and bonds
- **Shift + Left click**: Add to selection / deselect
- **Shift + Double click**: Reset view
- **Middle button/scroll**: Zoom in/out
- **Right drag**: Pan view
- **Left drag**: Rotate molecule (primary axes)
- **Shift + Left drag**: Rotate molecule around Z-axis

### Keyboard Controls

#### View Navigation

- **h/j/k/l**: Rotate around axes (vim-style)
- **u/o**: Rotate around Z-axis
- **H/J/K/L**: Pan the view
- **n/m**: Zoom in/out

#### Atom and Bond Manipulation

- **b**: Cycle through bond types for selected bond
- **B**: Toggle bond (add/remove) between selected atoms
- **d**: Toggle hydrogen display

#### Trajectory Navigation

- **[/]**: Previous/next frame
- **{/}**: First/last frame
- **-/=**: Jump to lowest/highest energy frame

#### Other Controls

- **p**: Toggle between energy and bond length plots
- **Shift+p**: Toggle bond propagation across frames
- **s**: Export current view as SVG
- **q or Escape**: Quit

## Features in Detail

### Bond Types

quickORTEP supports multiple bond types that can be cycled through with the 'b' key:

1. **Distance Bond**: Thin line for tracking distances
2. **Covalent Bond**: Standard bond representation
3. **NCI Bond**: Non-covalent interaction (dashed)
4. **TS Bond**: Transition state bond (special color)

### Trajectory Visualization

When viewing trajectory files (multi-frame XYZ):

- A small energy plot is shown in the bottom-right corner
- When a bond is selected, press 'p' to toggle between energy and bond length plots
- Use bracket keys ([, ]) to navigate through frames
- Use {, } to jump to first or last frames
- Use -, = to jump to lowest or highest energy frames

### Bond Editing and Propagation

Bond edits can be propagated across trajectory frames:

1. Select two atoms (click on first, then Shift+click on second)
2. Press 'B' to create a bond
3. With a bond selected, press 'b' to cycle through bond types
4. Enable bond propagation with Shift+p to apply these changes to all frames

### SVG Export

Press 's' to export the current view as an SVG file (`quickORTEP_export.svg`), which is useful for creating publication-quality figures.

## File Formats

### Supported Input Formats

- **XYZ files**: Standard XYZ format for single molecules or trajectories
- **Gaussian output files**: Log files from Gaussian calculations with energy information

## Advanced Usage

### Anti-aliasing

For better visual quality, use supersampling:

```bash
python main.py molecule.xyz 2  # 2x supersampling
python main.py molecule.xyz 4  # 4x supersampling (slower but higher quality)
```

### Customizing Appearance

Edit `config.py` to customize:
- Atom rendering style
- Bond colors and thicknesses
- Interaction settings
- Graph themes

## Troubleshooting

### X11 Connection Issues

If you encounter X11 connection issues:

```bash
# Check X11 forwarding settings if using SSH
ssh -X user@host

# Set display variable if needed
export DISPLAY=:0
```



