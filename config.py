# config.py

# Global conversion factor:
ANGSTROM_TO_PIXEL = 100

# ================================
# Atom Visualization Settings
# ================================
ATOM_STYLE = {
    # Scales an atom's covalent radius to its visual size.
    "scale": 0.4,  # previously SCALE_ATOM_SPHERE

    # Color of the atom's border.
    "border_color": (0, 0, 0),  # ATOM_BORDER_COLOR

    # Minimum radius in pixels to display an atom.
    "min_radius": 2,  # MIN_ATOM_RADIUS
}

# ================================
# Arc (ORTEP) Drawing Settings
# ================================
ARC_STYLE = {
    # Factor to flatten arcs to simulate sphere curvature.
    "flatten": 0.4,  # ARC_FLATTEN

    # Thickness for the ORTEP meridian arcs.
    "meridian_thickness": 2,  # ARC_MERIDIAN_THICKNESS
}

# ================================
# Bond Visualization Settings
# ================================
COVALENT_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.18,

    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,

    # Color for covalent bonds.
    "color": (0, 0, 0),
}

NCI_BOND = {
    # Base bond thickness in Ångströms (to be modified by thickness factor).
    "thickness": 0.18,

    # Bond segment length in Ångströms.
    "segment_length": 0.20,

    # Color for non-covalent interaction bonds.
    "color": (128, 128, 128),

    # Factor to reduce the bond thickness for NCIs.
    "thickness_factor": 0.4,
}

TS_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.18,

    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,

    # Color for covalent bonds.
    "color": (128, 128, 148),
}

DISTANCE_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.05,

    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,

    # Color for covalent bonds.
    "color": (128, 128, 128),
}

# ================================
# Highlighting Settings
# ================================
HIGHLIGHT = {
    # Color for highlighting selected objects.
    "color": (255, 0, 0),

    # Additional thickness added to bonds when highlighted.
    "thickness_delta": 2,
}

# ================================
# Viewer Interaction Settings
# ================================
VIEWER_INTERACTION = {
    # Rotation increment (in degrees) per key press.
    "rotation_increment_deg": 5,

    # Panning increment (in pixels) for view movement.
    "pan_increment": 10,

    # Zoom factor for key-based zoom adjustments.
    "key_zoom_factor": 1.05,

    # Zoom factor for mouse wheel zoom adjustments.
    "mouse_zoom_factor": 1.1,

    # Sensitivity for mouse-driven rotation.
    "mouse_rotation_sensitivity": 0.5,
}

# graph theme

MINIMAL_THEME = {
    "background_color": (255, 255, 255),
    "axis_color": (0, 0, 0),
    "tick_color": (0, 0, 0),
    "data_line_color": (255, 0, 0),
    "indicator_color": (0, 0, 255),
    "font_color": (0, 0, 0),
    "font_size": 12,
    "tick_length": 5,
    "line_thickness": {
        "axis": 2,
        "data": 2,
        "tick": 2,
        "indicator": 1,
    },
    "margin": {
         "left": 30,
         "right": 10,
         "top": 10,
         "bottom": 20
    }
}

FULL_THEME = {
    "background_color": (255, 255, 255),
    "axis_color": (0, 0, 0),
    "tick_color": (0, 0, 0),
    "data_line_color": (255, 0, 0),
    "indicator_color": (0, 0, 255),
    "font_color": (0, 0, 0),
    "font_size": 12,
    "tick_length": 5,
    "line_thickness": {
        "axis": 2,
        "data": 2,
        "tick": 2,
        "indicator": 1,
    },
    "margin": {
         "left": 50,
         "right": 20,
         "top": 20,
         "bottom": 40
    }
}
