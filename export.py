#!/usr/bin/env python3
"""
export.py

This module provides export functionality for the MoleculeViewer.
It currently supports exporting the rendered view to an SVG file.
"""

from config import PRINT_THEME, ATOM_STYLE, COVALENT_BOND, NCI_BOND, TS_BOND, DISTANCE_BOND, HIGHLIGHT
from x11view.svg_canvas import SVGCanvas

def export_xyz(molecule, filename="export.xyz", message_service=None):
    """
    Export a molecule to an XYZ file.

    Parameters:
        molecule: The Molecule instance to export
        filename (str): The name of the output XYZ file
        message_service: Optional MessageService instance for sending status messages

    Returns:
        The filename of the exported XYZ file
    """
    # Get XYZ format as string
    xyz_content = molecule.to_xyz()

    # Write to file
    with open(filename, 'w') as f:
        f.write(xyz_content)

    # Send a status message if message_service is provided
    if message_service:
        message_service.log_info(f"XYZ export saved to {filename}")

    return filename

def export_svg(canvas, renderer, molecule, view_params, filename="quickORTEP_export.svg", message_service=None):
    """
    Export the current view to an SVG file using the PRINT_THEME.

    Parameters:
        canvas: The canvas object, providing width and height.
        renderer: The molecule renderer.
        molecule: The ORTEP_Molecule instance.
        view_params: The view parameters (scaling, offset, etc.)
        filename (str): The name of the output SVG file.
        message_service: Optional MessageService instance for sending status messages.
        
    Returns:
        The filename of the exported SVG.
    """
    # Create SVG canvas with print theme background color
    svg_canvas = SVGCanvas(
        width=canvas.width, 
        height=canvas.height,
        background_color=PRINT_THEME["background_color"]
    )
    
    # Save original theme-dependent settings
    original_settings = {
        'atom_border_color': ATOM_STYLE["border_color"],
        'covalent_bond_color': COVALENT_BOND["color"],
        'nci_bond_color': NCI_BOND["color"],
        'ts_bond_color': TS_BOND["color"],
        'distance_bond_color': DISTANCE_BOND["color"],
        'highlight_color': HIGHLIGHT["color"]
    }
    
    # Apply print theme settings
    ATOM_STYLE["border_color"] = PRINT_THEME["atom_border_color"]
    COVALENT_BOND["color"] = PRINT_THEME["covalent_bond_color"]
    NCI_BOND["color"] = PRINT_THEME["nci_bond_color"]
    TS_BOND["color"] = PRINT_THEME["ts_bond_color"]
    DISTANCE_BOND["color"] = PRINT_THEME["distance_bond_color"]
    HIGHLIGHT["color"] = PRINT_THEME["highlight_color"]
    
    # For SVG export, create a copy of view_params to customize fog settings
    from config import FOG_STYLE
    import copy
    svg_view_params = copy.copy(view_params)
    
    # Maintain fog_mode and parameters from viewer, but use print theme background color
    # for the fog_color when needed
    if hasattr(view_params, 'fog_mode') and view_params.fog_mode > 0:
        # Override FOG_STYLE["color"] with PRINT_THEME's background color for export
        # but keep the original fog mode, factors, and density
        FOG_STYLE_original_color = FOG_STYLE["color"]
        FOG_STYLE["color"] = PRINT_THEME["background_color"]
    
    # Render with print theme
    view_params.as_float = True
    renderer.draw_molecule(svg_canvas, molecule, view_params)
    view_params.as_float = False
    
    # Restore original FOG_STYLE color if we changed it
    if hasattr(view_params, 'fog_mode') and view_params.fog_mode > 0:
        FOG_STYLE["color"] = FOG_STYLE_original_color
    
    # Restore original settings
    ATOM_STYLE["border_color"] = original_settings['atom_border_color']
    COVALENT_BOND["color"] = original_settings['covalent_bond_color']
    NCI_BOND["color"] = original_settings['nci_bond_color']
    TS_BOND["color"] = original_settings['ts_bond_color']
    DISTANCE_BOND["color"] = original_settings['distance_bond_color']
    HIGHLIGHT["color"] = original_settings['highlight_color']
    
    svg_canvas.flush(filename)
    
    # Send a status message if message_service is provided
    if message_service:
        message_service.log_info(f"SVG export saved to {filename}")
    
    return filename
