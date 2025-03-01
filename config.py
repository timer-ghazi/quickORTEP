# config.py
# Select active theme here (change this to switch themes)

# ================================
# Theme Definitions
# ================================

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
    "atom_border_color": (0, 255, 0),       # Green borders
    
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
}

# Options: LIGHT_THEME, DARK_THEME, PRINT_THEME
# CURRENT_THEME = LIGHT_THEME  # Set to DARK_THEME for a dark mode
# CURRENT_THEME = DARK_THEME 
CURRENT_THEME = MATRIX_THEME

# Global conversion factor:
ANGSTROM_TO_PIXEL = 100

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
    
    # Color palette for atoms
    "color_palette": "Rasmol",
    
    # Fallback color for unknown atoms
    "fallback_color": (200, 200, 200),
}

# ================================
# Arc (ORTEP) Drawing Settings
# ================================
ARC_STYLE = {
    # Factor to flatten arcs to simulate sphere curvature.
    "flatten": 0.4,
    
    # Thickness for the ORTEP meridian arcs.
    "meridian_thickness": 2,
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
}

TS_BOND = {
    # Bond thickness in Ångströms.
    "thickness": 0.18,
    
    # Bond segment length in Ångströms used for splitting bonds.
    "segment_length": 0.22,
    
    # Color for TS bonds.
    "color": CURRENT_THEME["ts_bond_color"],
    
    # Minimum thickness in pixels
    "min_thickness_px": 1,
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
