#!/usr/bin/env python3
"""
export.py

This module provides export functionality for the MoleculeViewer.
It currently supports exporting the rendered view to an SVG file.
"""

from config import PRINT_THEME
from x11view.svg_canvas import SVGCanvas

def export_svg(canvas, renderer, molecule, view_params, filename="quickORTEP_export.svg"):
    """
    Export the current view to an SVG file.

    Parameters:
        canvas: The canvas object, providing width and height.
        renderer: The molecule renderer.
        molecule: The ORTEP_Molecule instance.
        view_params: The view parameters (scaling, offset, etc.)
        filename (str): The name of the output SVG file.
        
    Returns:
        The filename of the exported SVG.
    """
    # Create SVG canvas with print theme background color
    svg_canvas = SVGCanvas(
        width=canvas.width, 
        height=canvas.height,
        background_color=PRINT_THEME["background_color"]
    )
    
    view_params.as_float = True
    renderer.draw_molecule(svg_canvas, molecule, view_params)
    view_params.as_float = False
    svg_canvas.flush(filename)
    return filename