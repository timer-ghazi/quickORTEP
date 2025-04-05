# config.py
# Select active theme here (change this to switch themes)

# ================================
# Theme Definitions
# ================================

# Sci-Fi Holographic Theme
SCIFI_THEME = {
    "background_color": (8, 12, 20),        # Dark blue-black background
    "default_text_color": (0, 255, 170),    # Bright cyan-teal text
    "default_line_color": (0, 179, 136),    # Medium teal for lines

    # Bond colors
 #   "covalent_bond_color": (0, 255, 170),   # Bright cyan for primary bonds
 #   "covalent_bond_color": (68, 221, 204),   # Cyan for normal messages
    "covalent_bond_color": (54, 177, 163),   # Darker cyan

    "nci_bond_color": (0, 179, 136),        # Medium teal for non-covalent
    "ts_bond_color": (68, 221, 170),        # Teal-cyan for TS bonds
    "distance_bond_color": (100, 255, 200), # Light cyan for distance

    # Atom elements - keep standard element colors but add cyan glow effect
    "atom_border_color": (0, 221, 170),     # Cyan borders

    # UI elements
    "hud_text_color": (0, 255, 170),        # Bright cyan
    "message_panel_bg": (10, 21, 32),       # Dark blue panel background
    "message_info_color": (68, 221, 204),   # Cyan for normal messages
    "message_warning_color": (255, 255, 85), # Yellow for warnings
    "message_error_color": (255, 68, 68),   # Red for errors

    # Highlight colors
    "highlight_color": (51, 153, 255),      # Blue highlight

    # Graph colors
    "graph_axis_color": (0, 179, 136),      # Teal axes
    "graph_tick_color": (0, 153, 136),      # Slightly darker teal ticks
    "graph_data_color": (0, 255, 170),      # Bright cyan for data
    "graph_indicator_color": (119, 255, 221), # Light cyan for indicators
    
    # Grid settings
    "show_grid": True,                      # Show grid by default
#    "grid_major_color": (5, 16, 21),       # #051015 almost invisible
#    "grid_major_color": (10, 61, 77),      # #071f28 
    "grid_major_color": (6, 24, 31),        # #06181f
    "grid_minor_color": (5, 16, 21),        # #051015 almost invisible
    
    # Axes visibility
    "show_axes": False,                     # Hide axes by default
}

# Light theme (default, white background)
LIGHT_THEME = {
    "background_color": (255, 255, 255),    # White background
    "default_text_color": (0, 0, 0),        # Black text
    "default_line_color": (0, 0, 0),        # Black lines
    
    # Bond colors
    "covalent_bond_color": (0, 0, 0),       # Black
    "nci_bond_color": (128, 128, 128),      # Gray
    "ts_bond_color": (128, 128, 148),       # Grayish blue
    "distance_bond_color": (128, 128, 128), # Gray
    
    # Atom elements
    "atom_border_color": (0, 0, 0),         # Black
    
    # UI elements
    "hud_text_color": (0, 0, 0),            # Black
    "message_panel_bg": (240, 240, 240),    # Light gray
    "message_info_color": (0, 0, 0),        # Black
    "message_warning_color": (200, 150, 0), # Amber
    "message_error_color": (200, 0, 0),     # Red
    
    # Highlight colors
    "highlight_color": (255, 0, 0),         # Red
    
    # Graph colors
    "graph_axis_color": (0, 0, 0),          # Black
    "graph_tick_color": (0, 0, 0),          # Black
    "graph_data_color": (255, 0, 0),        # Red
    "graph_indicator_color": (0, 0, 255),   # Blue
    
    # Grid settings
    "show_grid": False,                     # Hide grid by default
    "grid_major_color": (200, 200, 200),    # Light gray for major grid lines
    "grid_minor_color": (220, 220, 220),    # Lighter gray for minor grid lines
    
    # Axes visibility
    "show_axes": False,                     # Hide axes by default
}

# Dark theme (black background)
DARK_THEME = {
    "background_color": (0, 0, 0),          # Black background
    "default_text_color": (220, 220, 220),  # Light gray text
    "default_line_color": (220, 220, 220),  # Light gray lines
    
    # Bond colors
    "covalent_bond_color": (220, 220, 220), # Light gray 
    "nci_bond_color": (180, 180, 180),      # Medium gray
    "ts_bond_color": (170, 170, 200),       # Light blue-gray
    "distance_bond_color": (180, 180, 180), # Medium gray
    
    # Atom elements
    "atom_border_color": (220, 220, 220),   # Light gray
    
    # UI elements
    "hud_text_color": (220, 220, 220),      # Light gray
    "message_panel_bg": (40, 40, 40),       # Dark gray
    "message_info_color": (220, 220, 220),  # Light gray
    "message_warning_color": (255, 200, 0), # Bright amber
    "message_error_color": (255, 80, 80),   # Bright red
    
    # Highlight colors
    "highlight_color": (255, 80, 80),       # Bright red
    
    # Graph colors
    "graph_axis_color": (180, 180, 180),    # Medium gray
    "graph_tick_color": (180, 180, 180),    # Medium gray
    "graph_data_color": (255, 100, 100),    # Light red
    "graph_indicator_color": (100, 100, 255),  # Light blue
    
    # Grid settings
    "show_grid": True,                      # Show grid by default
    "grid_major_color": (10, 61, 77),       # #0a3d4d - very dark teal/blue-green
    "grid_minor_color": (15, 71, 87),       # Slightly lighter shade for minor lines
    
    # Axes visibility
    "show_axes": False,                     # Hide axes by default
}

# Print-optimized theme (for SVG export)
PRINT_THEME = {
    "background_color": (255, 255, 255),    # White background
    "default_text_color": (0, 0, 0),        # Black text
    "default_line_color": (0, 0, 0),        # Black lines
    
    # Bond colors
    "covalent_bond_color": (0, 0, 0),       # Black
    "nci_bond_color": (100, 100, 100),      # Darker gray for better contrast in print
    "ts_bond_color": (80, 80, 120),         # Darker blue-gray
    "distance_bond_color": (100, 100, 100), # Darker gray
    
    # Atom elements
    "atom_border_color": (0, 0, 0),         # Black
    
    # UI elements
    "hud_text_color": (0, 0, 0),            # Black
    "message_panel_bg": (245, 245, 245),    # Very light gray
    "message_info_color": (0, 0, 0),        # Black
    "message_warning_color": (180, 90, 0),  # Dark amber
    "message_error_color": (180, 0, 0),     # Dark red
    
    # Highlight colors
    "highlight_color": (200, 0, 0),         # Darker red
    
    # Graph colors
    "graph_axis_color": (0, 0, 0),          # Black
    "graph_tick_color": (0, 0, 0),          # Black
    "graph_data_color": (200, 0, 0),        # Darker red
    "graph_indicator_color": (0, 0, 200),   # Darker blue
    
    # Grid settings
    "show_grid": True,                      # Show grid by default
    "grid_major_color": (220, 220, 220),    # Light gray for major grid lines
    "grid_minor_color": (235, 235, 235),    # Lighter gray for minor grid lines
    
    # Axes visibility
    "show_axes": False,                     # Hide axes by default
}

# Matrix theme (green-on-black cyberpunk aesthetic)
MATRIX_THEME = {
    "background_color": (0, 0, 0),          # Pure black background
    "default_text_color": (0, 215, 0),      # Bright Matrix green text
    "default_line_color": (0, 150, 0),      # Slightly darker green for lines
    
    # Bond colors
    "covalent_bond_color": (0, 200, 0),     # Bright green for primary bonds
    "nci_bond_color": (0, 150, 0),          # Medium green for non-covalent
    "ts_bond_color": (20, 200, 100),        # Green-teal for TS bonds
    "distance_bond_color": (100, 255, 100), # Light green for distance
    
    # Atom elements - keep standard element colors but add green tint
    "atom_border_color": (0, 200, 0),       # Green borders
    
    # UI elements
    "hud_text_color": (0, 255, 0),          # Bright Matrix green
    "message_panel_bg": (10, 30, 10),       # Very dark green-black
    "message_info_color": (0, 255, 0),      # Bright green for normal messages
    "message_warning_color": (255, 255, 0), # Yellow for warnings
    "message_error_color": (255, 0, 0),     # Red for errors
    
    # Highlight colors
    "highlight_color": (0, 255, 160),       # Teal-green highlight
    
    # Graph colors
    "graph_axis_color": (0, 200, 0),        # Green axes
    "graph_tick_color": (0, 180, 0),        # Slightly darker green ticks
    "graph_data_color": (0, 255, 120),      # Bright green-teal for data
    "graph_indicator_color": (180, 255, 180),  # Light green for indicators
    
    # Grid settings
    "show_grid": True,                      # Show grid by default
    "grid_major_color": (0, 30, 0),         # Dark green for major grid lines
    "grid_minor_color": (0, 20, 0),         # Darker green for minor grid lines
    
    # Axes visibility
    "show_axes": False,                     # Hide axes by default
}

# Options: LIGHT_THEME, DARK_THEME, PRINT_THEME
# CURRENT_THEME = LIGHT_THEME  # Set to DARK_THEME for a dark mode
# CURRENT_THEME = DARK_THEME 
CURRENT_THEME = SCIFI_THEME

# Global conversion factor:
ANGSTROM_TO_PIXEL = 100

# ================================
# Energy Unit Settings
# ================================
ENERGY_UNITS = {
    "hartree": {
        "name": "Hartree",
        "symbol": "Ha",
        "to_hartree": 1.0,
        "to_kcal_mol": 627.5095,
        "to_kj_mol": 2625.5,
        "to_ev": 27.2114
    },
    "kcal_mol": {
        "name": "kcal/mol",
        "symbol": "kcal/mol",
        "to_hartree": 1.0 / 627.5095,
        "to_kcal_mol": 1.0,
        "to_kj_mol": 4.184,
        "to_ev": 0.0433641
    },
    "kj_mol": {
        "name": "kJ/mol",
        "symbol": "kJ/mol",
        "to_hartree": 1.0 / 2625.5,
        "to_kcal_mol": 1.0 / 4.184,
        "to_kj_mol": 1.0,
        "to_ev": 0.0103642
    },
    "ev": {
        "name": "eV",
        "symbol": "eV",
        "to_hartree": 1.0 / 27.2114,
        "to_kcal_mol": 1.0 / 0.0433641,
        "to_kj_mol": 1.0 / 0.0103642,
        "to_ev": 1.0
    }
}

# Default energy unit for display
DEFAULT_ENERGY_UNIT = "kcal_mol"

# ================================
# Grid Settings
# ================================
GRID_SETTINGS = {
    "major_spacing": 1.0,      # Major grid lines every 1.0 Å
    "minor_divisions": 2,      # 4 divisions = 0.25 Å between minor lines
    "major_thickness": 2,      # 1 pixel thick
    "minor_thickness": 1,      # 1 pixel thick
}

# ================================
# Canvas Settings
# ================================
CANVAS_SETTINGS = {
    "background_color": CURRENT_THEME["background_color"],
}

# ================================
# Atom Visualization Settings
# ================================
ATOM_STYLE = {
    # Scales an atom's covalent radius to its visual size.
    "scale": 0.4,
    
    # Color of the atom's border.
    "border_color": CURRENT_THEME["atom_border_color"],
    
    # Minimum radius in pixels to display an atom.
    "min_radius": 2,
    
    # Border thickness in Ångströms (not pixels).
    "border_thickness": 0.02,
    
    # Minimum border thickness in pixels (for very low zoom levels).
    "min_border_thickness_px": 1,
    
    # Color palette for atoms
    "color_palette": "Rasmol",
    # "color_palette": "GreenRasmol",
    
    # Fallback color for unknown atoms
    "fallback_color": (200, 200, 200),
    
    # 3D highlight settings
    "highlight": {
        "enabled": True,           # Turn the effect on/off
        "size_ratio": 0.3,         # Size relative to atom radius
        "offset_ratio": 0.4,       # Offset from center (towards upper-left)
        "brightness_factor": 1.5,  # For calculated highlight color
        "use_white_for_dark": True # Use white for dark atoms
    },
    
    # 3D shadow settings
    "shadow": {
        "enabled": True,        # Turn the effect on/off
        "offset_x": 2,          # Horizontal shadow offset (positive = right)
        "offset_y": 3,          # Vertical shadow offset (positive = down)
        "radius_factor": 1.1,   # Shadow size relative to atom radius
        "darkness": 40          # Shadow color darkness (0-255, lower = darker)
    },
    
    # Enhanced 3D lighting settings
    "lighting": {
        "enabled": False,              # Turn enhanced lighting on/off
        "direction": (0.5, -0.5, 1),  # Global light direction vector (x, y, z)
        "ambient": 0.3,               # Ambient light intensity (0.0-1.0)
        "diffuse": 0.6,               # Diffuse reflection intensity (0.0-1.0)
        "specular": 0.4,              # Specular reflection intensity (0.0-1.0)
        "shininess": 16,              # Shininess power factor (higher = smaller highlight)
        "segments": 12,               # Number of segments for lighting calculation
        "gradient_steps": 8           # Number of gradient steps from dark to light
    },
    
    # Fog effect settings
    "fog": {
        "enabled": False,             # Turn fog effect on/off (toggled with enhanced lighting)
        "start_distance": 5.0,        # Distance where fog starts (in Å)
        "end_distance": 20.0,         # Distance where fog is maximum (in Å)
        "max_intensity": 0.7,         # Maximum fog intensity (0.0-1.0)
        "exponent": 1.0,              # Exponent for fog calculation (1.0 = linear, >1.0 = exponential)
    },
}

# ================================
# Arc (ORTEP) Drawing Settings
# ================================
ARC_STYLE = {
    # Factor to flatten arcs to simulate sphere curvature.
    "flatten": 0.4,
    
    # Thickness for the ORTEP meridian arcs in Ångströms (not pixels).
    "meridian_thickness": 0.03,
    
    # Minimum thickness in pixels (for very low zoom levels).
    "min_meridian_thickness_px": 1,
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
    "color": CURRENT_THEME["covalent_bond_color"],
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
    
    # Tapering configuration
    "enable_taper": True,         # Turn tapering on/off
    "min_taper_factor": 0.5,      # Minimum thickness as fraction of normal (thinnest part)
    "taper_range": 0.5,           # Difference between thickest and thinnest parts
}

NCI_BOND = {
    # Base bond thickness in Ångströms (to be modified by thickness factor).
    "thickness": 0.18,
    
    # Bond segment length in Ångströms.
    "segment_length": 0.20,
    
    # Color for non-covalent interaction bonds.
    "color": CURRENT_THEME["nci_bond_color"],
    
    # Factor to reduce the bond thickness for NCIs.
    "thickness_factor": 0.4,
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
    
    # Tapering configuration
    "enable_taper": True,         # Turn tapering on/off
    "min_taper_factor": 0.5,      # Minimum thickness as fraction of normal (thinnest part)
    "taper_range": 0.5,           # Difference between thickest and thinnest parts
}

TS_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.12,
    
    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,
    
    # Color for TS bonds.
    "color": CURRENT_THEME["ts_bond_color"],
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
    
    # Tapering configuration
    "enable_taper": True,         # Turn tapering on/off
    "min_taper_factor": 0.5,      # Minimum thickness as fraction of normal (thinnest part)
    "taper_range": 0.5,           # Difference between thickest and thinnest parts
}

DISTANCE_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.05,
    
    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,
    
    # Color for distance bonds.
    "color": CURRENT_THEME["distance_bond_color"],
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
    
    # Tapering configuration
    "enable_taper": True,         # Turn tapering on/off
    "min_taper_factor": 0.5,      # Minimum thickness as fraction of normal (thinnest part)
    "taper_range": 0.5,           # Difference between thickest and thinnest parts
}

# ================================
# Vector Visualization Settings
# ================================
VECTOR_STYLE = {
    # Vector shaft thickness in Ångströms
    "thickness": 0.05,
    
    # Vector segment length in Ångströms (for Z-ordering)
    "segment_length": 0.22,
    
    # Default color for vectors (can be overridden per vector)
    "color": CURRENT_THEME["default_line_color"],
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
    
    # Arrowhead settings
    "arrowhead": {
        # Length of arrowhead in Ångströms
        "length": 0.3,
        
        # Width of arrowhead as a ratio of its length
        "width_ratio": 0.6,
        
        # Whether arrowhead should be filled (True) or outline only (False)
        "filled": True,
        
        # Color for arrowhead (None = same as vector shaft)
        "color": None
    }
}

# ================================
# Coordinate Axes Settings
# ================================
AXES_STYLE = {
    # Length of each axis in Ångströms
    "length": 2.0,
    
    # Colors for X, Y, Z axes
    "x_color": (255, 0, 0),    # Red
    "y_color": (0, 255, 0),    # Green
    "z_color": (0, 0, 255),    # Blue
    
    # Thickness multiplier compared to VECTOR_STYLE["thickness"]
    "thickness_multiplier": 1.2,
    
    # Show labels (X, Y, Z) at the end of each axis
    "show_labels": True,
    
    # Font size for labels
    "label_font_size": 12,
    
    # Distance of label from axis end in Ångströms
    "label_offset": 0.2
}

# ================================
# Highlighting Settings
# ================================
HIGHLIGHT = {
    # Color for highlighting selected objects.
    "color": CURRENT_THEME["highlight_color"],
    
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
    
    # Key to toggle axes visibility
    "axes_toggle_key": "a",
    
    # Current zoom level (updated dynamically)
    "current_zoom": 100.0,
}

# ================================
# HUD Panel Settings
# ================================
HUD_STYLE = {
    "x": 10,                              # X-coordinate for HUD text
    "y_offset": 75,                       # Vertical offset from bottom
    "line_spacing": 16,                   # Spacing between HUD lines
    "font_size": 14,                      # Font size for HUD text
    "color": CURRENT_THEME["hud_text_color"],  # Color of HUD text
}

# ================================
# Message Panel Settings
# ================================
MESSAGE_PANEL_STYLE = {
    "x": 10,                              # X-coordinate for messages
    "y_offset": 10,                       # Vertical offset from bottom
    "line_spacing": 16,                   # Spacing between messages
    "font_size": 12,                      # Font size for messages
    "bg_color": CURRENT_THEME["message_panel_bg"],  # Background color
    "padding": 5,                         # Padding around text
    "max_messages": 3,                    # Maximum number of messages to show
}

# ================================
# Message Types Settings
# ================================
MESSAGE_TYPES = {
    "info": {
        "prefix": "",
        "color": CURRENT_THEME["message_info_color"],
    },
    "warning": {
        "prefix": "WARNING",
        "color": CURRENT_THEME["message_warning_color"],
    },
    "error": {
        "prefix": "ERROR",
        "color": CURRENT_THEME["message_error_color"],
    }
}

# ================================
# Graph Settings
# ================================
GRAPH_SETTINGS = {
    "thumbnail_width": 150,               # Width of thumbnail graph
    "thumbnail_height": 80,               # Height of thumbnail graph
    "thumbnail_margin": 20,               # Margin from window edges
    "indicator_min_radius": 3,            # Minimum radius for indicator
    "indicator_size_factor": 0.04,        # Factor for indicator size
    "minor_tick_factor": 0.5,             # Factor for minor tick length
}

# ================================
# Graph Themes
# ================================
MINIMAL_THEME = {
    "background_color": CURRENT_THEME["background_color"],
    "axis_color": CURRENT_THEME["graph_axis_color"],
    "tick_color": CURRENT_THEME["graph_tick_color"],
    "data_line_color": CURRENT_THEME["graph_data_color"],
    "indicator_color": CURRENT_THEME["graph_indicator_color"],
    "font_color": CURRENT_THEME["default_text_color"],
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
    "background_color": CURRENT_THEME["background_color"],
    "axis_color": CURRENT_THEME["graph_axis_color"],
    "tick_color": CURRENT_THEME["graph_tick_color"],
    "data_line_color": CURRENT_THEME["graph_data_color"],
    "indicator_color": CURRENT_THEME["graph_indicator_color"],
    "font_color": CURRENT_THEME["default_text_color"],
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
