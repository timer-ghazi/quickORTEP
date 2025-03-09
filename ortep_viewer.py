#!/usr/bin/env python3
"""
ortep_viewer.py

The MoleculeViewer class manages an ORTEP_Molecule, its rendering,
and handles user events. It now delegates trajectory conversion,
selection, event handling, and export functionality to helper modules.
It also supports bond propagation across trajectory frames.

Normal mode visualization is supported through vectors representing
vibrational displacements.
"""

import sys
import time
import math
import numpy as np
from Xlib import XK, X
from x11view.window import X11Window, X11CanvasBasic, X11CanvasSS
from hud import HUDPanel
from ortep_renderer import ORTEP_MoleculeRenderer
from config import CANVAS_SETTINGS, VIEWER_INTERACTION, MINIMAL_THEME
from config import HUD_STYLE, MESSAGE_PANEL_STYLE, MESSAGE_TYPES, GRAPH_SETTINGS
from config import CURRENT_THEME, GRID_SETTINGS, ENERGY_UNITS, DEFAULT_ENERGY_UNIT
from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from message_service import MessageService
from message_panel import MessagePanel
from graph_viewer import GraphViewer
from bond_edit_tracker import BondEditTracker
from vectors import Vector

# Import our helper modules.
from selection_manager import _SelectionManager
from trajectory_manager import _TrajectoryManager
from event_handler import _EventHandler
from export import export_svg, export_xyz
from help_text import HELP_TEXT


# ===============================================================
# Helper Classes
# ===============================================================

class ViewParams:
    def __init__(self, rx=0.0, ry=0.0, rz=0.0,
                 scale=100.0, x_offset=400, y_offset=300):
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.scale = scale
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.as_float = False


class _NormalModeManager:
    """
    Handles normal mode visualization and operations.
    """
    def __init__(self, viewer):
        self.viewer = viewer
        self.show_normal_modes = False
        self.current_normal_mode_index = 0
        self.normal_mode_scale_factor = 1.0
        self.has_normal_modes = False
        self.normal_mode_data = None
        self.normal_mode_frame_index = None
    
    def initialize_from_trajectory(self, trajectory):
        """Initialize normal mode data from trajectory metadata."""
        self.has_normal_modes = False
        self.normal_mode_data = None
        self.normal_mode_frame_index = None
        
        if 'frequency_data' in trajectory.metadata:
            self.normal_mode_data = trajectory.metadata['frequency_data']
            self.normal_mode_frame_index = trajectory.metadata.get('freq_data_frame_index', None)
            if self.normal_mode_data and self.normal_mode_frame_index is not None:
                self.has_normal_modes = True
                self.viewer.message_service.log_info(f"Normal mode data available for frame {self.normal_mode_frame_index}")

    def toggle_normal_modes(self):
        """Toggle display of normal mode vectors."""
        if not self.has_normal_modes:
            self.viewer.message_service.log_info("No normal modes available")
            return

        # Only allow showing normal modes on the frame that has the data
        if self.viewer.current_frame != self.normal_mode_frame_index:
            self.viewer.message_service.log_info(f"Normal modes only available for frame {self.normal_mode_frame_index}")
            return
        
        # Toggle display state
        self.show_normal_modes = not self.show_normal_modes
        
        if self.show_normal_modes:
            # Create vectors for the current normal mode
            self.create_normal_mode_vectors()
            mode = self.normal_mode_data.modes[self.current_normal_mode_index]
            freq_text = f"{abs(mode.frequency):.1f} cm⁻¹"
            if mode.frequency < 0:
                freq_text = f"{freq_text} (imaginary)"
            self.viewer.message_service.log_info(f"Showing normal mode {self.current_normal_mode_index + 1}/{len(self.normal_mode_data.modes)}: {freq_text}")
        else:
            # Clear normal mode vectors
            self.clear_normal_mode_vectors()
            self.viewer.message_service.log_info("Normal mode vectors hidden")
        
        self.viewer.redraw()

    def cycle_normal_mode(self, direction):
        """
        Cycle to the next or previous normal mode.
        
        Parameters:
            direction (int): 1 for next, -1 for previous
        """
        if not self.has_normal_modes or not self.show_normal_modes:
            return
            
        # Update the current normal mode index
        num_modes = len(self.normal_mode_data.modes)
        self.current_normal_mode_index = (self.current_normal_mode_index + direction) % num_modes
        
        # Create vectors for the new current mode
        self.create_normal_mode_vectors()
        
        # Log information about the current mode
        mode = self.normal_mode_data.modes[self.current_normal_mode_index]
        freq_text = f"{abs(mode.frequency):.1f} cm⁻¹"
        if mode.frequency < 0:
            freq_text = f"{freq_text} (imaginary)"
        self.viewer.message_service.log_info(f"Normal mode {self.current_normal_mode_index + 1}/{num_modes}: {freq_text}")
        
        self.viewer.redraw()

    def adjust_normal_mode_scale(self, factor):
        """
        Adjust the scale factor for normal mode vectors.
        
        Parameters:
            factor (float): Multiplier to adjust the scale
        """
        if not self.has_normal_modes or not self.show_normal_modes:
            return
            
        # Update the scale factor
        old_scale = self.normal_mode_scale_factor
        self.normal_mode_scale_factor *= factor
        
        # Clamp to reasonable values
        min_scale = 0.1
        max_scale = 50.0
        self.normal_mode_scale_factor = max(min_scale, min(max_scale, self.normal_mode_scale_factor))
        
        # Recreate vectors with the new scale
        if old_scale != self.normal_mode_scale_factor:
            self.create_normal_mode_vectors()
            self.viewer.message_service.log_info(f"Normal mode scale: {self.normal_mode_scale_factor:.1f}")
            self.viewer.redraw()

    def create_normal_mode_vectors(self):
        """Create vectors representing the current normal mode."""
        # Clear any existing normal mode vectors
        self.clear_normal_mode_vectors()
        
        if not self.has_normal_modes or not self.show_normal_modes:
            return
            
        # Get the current normal mode
        mode = self.normal_mode_data.modes[self.current_normal_mode_index]
        
        # Set colors based on whether the mode is imaginary
        if mode.frequency < 0:
            # Imaginary frequency - use red
            vector_color = (255, 0, 0)
        else:
            # Real frequency - use blue to green gradient based on displacement magnitude
            vector_color = (0, 128, 255)  # Default blue
        
        # Create vectors for each atom's displacement
        displacements = mode.displacements
        
        # Find maximum displacement for scaling
        max_displacement = 0.0
        for dx, dy, dz in displacements:
            magnitude = np.sqrt(dx*dx + dy*dy + dz*dz)
            max_displacement = max(max_displacement, magnitude)
        
        # Base scale - adjusted empirically
        scale_base = 3.0 / max_displacement if max_displacement > 0 else 1.0
        scale = scale_base * self.normal_mode_scale_factor
        
        # Create vectors for each atom
        for i, atom in enumerate(self.viewer.ortep_mol.atoms):
            if i < len(displacements):
                dx, dy, dz = displacements[i]
                
                # Skip negligible displacements
                magnitude = np.sqrt(dx*dx + dy*dy + dz*dz)
                if magnitude < 1e-6:
                    continue
                
                # Create the vector
                start_point = (atom.x, atom.y, atom.z)
                end_point = (
                    atom.x + dx * scale,
                    atom.y + dy * scale,
                    atom.z + dz * scale
                )
                
                # For non-imaginary modes, color by displacement magnitude relative to max
                if mode.frequency >= 0:
                    # Scale from blue to green based on magnitude
                    intensity = min(1.0, magnitude / max_displacement)
                    r = 0
                    g = int(128 + 127 * intensity)
                    b = int(255 - 127 * intensity)
                    vector_color = (r, g, b)
                
                # Create and add the vector
                vector = Vector(start_point=start_point, end_point=end_point, color=vector_color)
                self.viewer.ortep_mol.add_vector(vector)

    def clear_normal_mode_vectors(self):
        """Clear all normal mode vectors."""
        # Remove all vectors that are not coordinate axes
        self.viewer.ortep_mol.vectors = [v for v in self.viewer.ortep_mol.vectors if hasattr(v, '_is_axis_vector')]
        
        # Recreate coordinate axes if they should be shown
        if self.viewer.show_axes and not self.viewer.ortep_mol.axis_system:
            self.viewer.ortep_mol.create_coordinate_axes()

    def update_info_message(self, lines):
        """Add normal mode information to the HUD lines."""
        if self.has_normal_modes:
            if self.viewer.current_frame == self.normal_mode_frame_index:
                if self.show_normal_modes:
                    mode = self.normal_mode_data.modes[self.current_normal_mode_index]
                    freq_text = f"{abs(mode.frequency):.1f} cm⁻¹"
                    if mode.frequency < 0:
                        freq_text = f"{freq_text} (imaginary)"
                    lines.append(f"Normal mode {self.current_normal_mode_index + 1}/{len(self.normal_mode_data.modes)}: {freq_text}")
                    lines.append(f"IR intensity: {mode.ir_intensity:.1f} KM/mol | Scale: {self.normal_mode_scale_factor:.1f}")
                    lines.append("Press ',' / '.' to cycle modes, '<' / '>' to adjust scale")
                else:
                    lines.append("Normal modes available. Press 'v' to view.")
            elif self.normal_mode_frame_index is not None:
                lines.append(f"Normal modes available on frame {self.normal_mode_frame_index}")
        
        return lines


class _GraphManager:
    """
    Manages energy and bond length graphs.
    """
    def __init__(self, viewer):
        self.viewer = viewer
        self.energy_graph = None
        self.graph_mode = "energy"  # Options: "energy", "bond_length"
        self.selected_bond_for_graph = None  # Will store the bond ID (atom1_idx, atom2_idx)

    def ensure_energy_graph(self):
        """Initialize or update the graph based on the current mode (energy or bond length)."""
        # Configure graph position in bottom right corner
        thumb_width = GRAPH_SETTINGS["thumbnail_width"]
        thumb_height = GRAPH_SETTINGS["thumbnail_height"]
        thumb_margin = GRAPH_SETTINGS["thumbnail_margin"]
        thumb_x = self.viewer.canvas.width - thumb_width - thumb_margin
        thumb_y = self.viewer.canvas.height - thumb_height - thumb_margin
        
        # If we're in bond length mode, ensure we have a selected bond
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph is None:
            # If no bond is selected, revert to energy mode
            self.graph_mode = "energy"
        
        if self.graph_mode == "energy":
            # Skip if no trajectory is available
            if self.viewer.trajectory is None:
                return
            
            # Get energy data using the energy_trajectory method with unit conversion
            energies, energy_info = self.viewer.trajectory.energy_trajectory(
                convert_if_hartrees=True,
                convert_to_unit=DEFAULT_ENERGY_UNIT
            )
            
            if len(energies) == 0 or np.isnan(energies).all():
                return
            
            # Check if coordinate or time metadata is available for the x-axis.
            coords, coord_info = self.viewer.trajectory.coordinate_trajectory(skip_none=False)
            if coord_info.get('field') is not None and len(coords) == len(energies):
                x_values = list(coords)
                x_axis_title = "Coord" if coord_info['field'] == "coord" else "Time"
            else:
                x_values = list(range(len(energies)))
                x_axis_title = "Frame"
            
            y_values = [e if not np.isnan(e) else 0.0 for e in energies]
            title = "Energy"
            y_axis_title = ""
        else:
            # Bond length graph
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            atom1 = next((a for a in self.viewer.ortep_mol.atoms if a.index == atom1_idx), None)
            atom2 = next((a for a in self.viewer.ortep_mol.atoms if a.index == atom2_idx), None)
            
            # Get bond length data
            x_values, y_values = self.calculate_bond_length_trajectory(atom1_idx, atom2_idx)
            
            # Prepare title with atom symbols
            if atom1 and atom2:
                title = f"{atom1.symbol}{atom1_idx}-{atom2.symbol}{atom2_idx}"
            else:
                title = "Bond Length"
            x_axis_title = "Frame"
            y_axis_title = ""
        
        # Create a custom minimal theme
        custom_minimal_theme = MINIMAL_THEME.copy()
        custom_minimal_theme["margin"] = {"left": 10, "right": 10, "top": 10, "bottom": 10}
        
        if self.energy_graph is None:
            # Create the graph
            self.energy_graph = GraphViewer(
                canvas=self.viewer.canvas,
                xdata=x_values,
                ydata=y_values,
                mode="minimal",
                current_frame=self.viewer.current_frame,
                region_x=thumb_x,
                region_y=thumb_y,
                region_width=thumb_width,
                region_height=thumb_height,
                x_axis_title=x_axis_title,
                y_axis_title=y_axis_title,
                title=title,
                theme=custom_minimal_theme
            )
        else:
            # Update existing graph position
            self.energy_graph.region_x = thumb_x
            self.energy_graph.region_y = thumb_y
            
            # Update data if the mode changed
            if self.energy_graph.title != title:
                self.energy_graph.update_data(
                    xdata=x_values,
                    ydata=y_values,
                    title=title,
                    y_axis_title=y_axis_title,
                    x_axis_title=x_axis_title
                )

    def calculate_bond_length_trajectory(self, atom1_idx, atom2_idx):
        """
        Calculate bond length between two atoms across all frames in the trajectory.
        
        Parameters:
            atom1_idx, atom2_idx (int): Atom indices
            
        Returns:
            tuple: (x_values, y_values) where x_values are frame numbers and 
                   y_values are bond lengths in Ångströms
        """
        if self.viewer.trajectory is None:
            return [], []
        
        x_values = list(range(self.viewer.total_frames))
        y_values = []
        
        for frame_idx in range(self.viewer.total_frames):
            # Get the molecule for this frame (without changing the current frame)
            mol = self.viewer.trajectory.get_frame(frame_idx)
            # Calculate distance between the atoms (0-indexed in molecule_nci)
            distance = mol.distance(atom1_idx-1, atom2_idx-1)
            y_values.append(distance)
        
        return x_values, y_values

    def toggle_graph_mode(self):
        """Toggle between energy and bond length visualization modes."""
        if self.viewer.selected_bonds and len(self.viewer.selected_bonds) == 1:
            bond = self.viewer.selected_bonds[0]
            atom1_idx = bond.atom1.index
            atom2_idx = bond.atom2.index
            
            if self.graph_mode == "energy":
                # Switch to bond length mode
                self.graph_mode = "bond_length"
                self.selected_bond_for_graph = (atom1_idx, atom2_idx)
                self.viewer.message_service.log_info(
                    f"Showing bond length for {bond.atom1.symbol}{atom1_idx}-{bond.atom2.symbol}{atom2_idx}"
                )
            else:
                # Switch back to energy mode
                self.graph_mode = "energy"
                self.viewer.message_service.log_info("Showing energy plot")
            
            # Refresh the graph
            self.energy_graph = None  # Force recreation
            self.ensure_energy_graph()
        elif self.graph_mode == "bond_length":
            # Switch back to energy mode if in bond length mode but no bond selected
            self.graph_mode = "energy"
            self.viewer.message_service.log_info("Showing energy plot")
            self.energy_graph = None  # Force recreation
            self.ensure_energy_graph()
        else:
            self.viewer.message_service.log_info("Select a bond first to toggle length graph")

    def update_info_message(self, lines):
        """Add graph information to the HUD lines."""
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph:
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            atom1 = next((a for a in self.viewer.ortep_mol.atoms if a.index == atom1_idx), None)
            atom2 = next((a for a in self.viewer.ortep_mol.atoms if a.index == atom2_idx), None)
            if atom1 and atom2:
                lines.append(f"Graph: Bond length {atom1.symbol}{atom1_idx}-{atom2.symbol}{atom2_idx}")
        else:
            # Add unit info to the graph description
            if self.viewer.trajectory:
                _, energy_info = self.viewer.trajectory.energy_trajectory(
                    convert_if_hartrees=True,
                    convert_to_unit=DEFAULT_ENERGY_UNIT
                )
                unit_symbol = ENERGY_UNITS.get(
                    energy_info['converted_unit'], 
                    {'symbol': energy_info['converted_unit']}
                )['symbol']
                if energy_info['normalized']:
                    lines.append(f"Graph: Rel. Energy ({unit_symbol})")
                else:
                    lines.append(f"Graph: Energy ({unit_symbol})")
            else:
                lines.append("Graph: Energy")
        
        return lines

    def check_selected_bond_exists(self):
        """Check if the selected bond for the graph still exists in the current frame."""
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph:
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            bond_exists = False
            for bond in self.viewer.ortep_mol.bonds:
                if ((bond.atom1.index == atom1_idx and bond.atom2.index == atom2_idx) or
                    (bond.atom1.index == atom2_idx and bond.atom2.index == atom1_idx)):
                    bond_exists = True
                    break
            
            if not bond_exists:
                # Bond doesn't exist in this frame, switch back to energy mode
                self.graph_mode = "energy"
                self.selected_bond_for_graph = None
                self.energy_graph = None  # Force recreation


class _ViewManager:
    """
    Handles view transformations and rendering operations.
    """
    def __init__(self, viewer):
        self.viewer = viewer
        self.last_canvas_width = viewer.canvas.width
        self.last_canvas_height = viewer.canvas.height
        self.show_grid = CURRENT_THEME.get("show_grid", False)
        self.show_axes = CURRENT_THEME.get("show_axes", False)
        self.show_help = False
        self.show_hydrogens = True

    def check_resize(self):
        """Check if the canvas has been resized."""
        if self.last_canvas_width != self.viewer.canvas.width or self.last_canvas_height != self.viewer.canvas.height:
            self.last_canvas_width = self.viewer.canvas.width
            self.last_canvas_height = self.viewer.canvas.height
            return True
        return False

    def draw_grid(self):
        """
        Draw a background grid with major and minor lines.
        Major lines are drawn every 1 Å, and minor lines every 0.25 Å.
        The grid remains fixed relative to the molecular origin.
        """
        # Skip if grid is not enabled in the current theme
        if not CURRENT_THEME.get("show_grid", False):
            return
            
        # Get grid settings
        major_spacing = GRID_SETTINGS["major_spacing"]  # Å
        minor_divisions = GRID_SETTINGS["minor_divisions"]
        minor_spacing = major_spacing / minor_divisions  # Å
        
        major_color = CURRENT_THEME.get("grid_major_color", (10, 61, 77))
        minor_color = CURRENT_THEME.get("grid_minor_color", (15, 71, 87))
        
        major_thickness = GRID_SETTINGS["major_thickness"]
        minor_thickness = GRID_SETTINGS["minor_thickness"]
        
        # Calculate the grid based on current view parameters
        scale = self.viewer.view_params.scale
        x_offset = self.viewer.view_params.x_offset
        y_offset = self.viewer.view_params.y_offset
        
        # Calculate visible area in Ångströms
        left_a = -x_offset / scale
        right_a = left_a + (self.viewer.canvas.width / scale)
        top_a = -y_offset / scale
        bottom_a = top_a + (self.viewer.canvas.height / scale)
        
        # Round to the nearest major grid line outside the viewport
        left_a_grid = math.floor(left_a / major_spacing) * major_spacing
        right_a_grid = math.ceil(right_a / major_spacing) * major_spacing
        top_a_grid = math.floor(top_a / major_spacing) * major_spacing
        bottom_a_grid = math.ceil(bottom_a / major_spacing) * major_spacing
        
        # Draw major grid lines (vertical)
        x_a = left_a_grid
        while x_a <= right_a_grid:
            x_px = int(x_a * scale + x_offset)
            self.viewer.canvas.draw_line(x_px, 0, x_px, self.viewer.canvas.height, 
                                thickness=major_thickness, color=major_color)
            x_a += major_spacing
        
        # Draw major grid lines (horizontal)
        y_a = top_a_grid
        while y_a <= bottom_a_grid:
            y_px = int(y_a * scale + y_offset)
            self.viewer.canvas.draw_line(0, y_px, self.viewer.canvas.width, y_px, 
                                thickness=major_thickness, color=major_color)
            y_a += major_spacing
        
        # Draw minor grid lines if enabled and visible enough
        if minor_divisions > 1 and scale > 20:  # Only draw minor lines when sufficiently zoomed in
            # Draw minor grid lines (vertical)
            x_a = left_a_grid
            while x_a <= right_a_grid:
                for i in range(1, minor_divisions):
                    minor_x_a = x_a + (i * minor_spacing)
                    minor_x_px = int(minor_x_a * scale + x_offset)
                    self.viewer.canvas.draw_dashed_line(minor_x_px, 0, minor_x_px, self.viewer.canvas.height, 
                                                thickness=minor_thickness, color=minor_color)
                x_a += major_spacing
            
            # Draw minor grid lines (horizontal)
            y_a = top_a_grid
            while y_a <= bottom_a_grid:
                for i in range(1, minor_divisions):
                    minor_y_a = y_a + (i * minor_spacing)
                    minor_y_px = int(minor_y_a * scale + y_offset)
                    self.viewer.canvas.draw_dashed_line(0, minor_y_px, self.viewer.canvas.width, minor_y_px, 
                                                thickness=minor_thickness, color=minor_color)
                y_a += major_spacing

    def draw_help_overlay(self):
        """
        Draw a rectangle plus the cheat-sheet text,
        using the same colors as the HUD.
        """
        # 1) Use the message panel's background color for the overlay rectangle
        overlay_color = MESSAGE_PANEL_STYLE["bg_color"]
    
        # 2) Use the HUD text color for the help text
        text_color = HUD_STYLE["color"]
    
        # Fill the entire window with this overlay color
        self.viewer.canvas.draw_filled_rect(0, 0, self.viewer.canvas.width, self.viewer.canvas.height, color=overlay_color)
    
        # 3) Draw the grid
        self.draw_grid()
    
        # 4) Decide positioning
        x_start = 50
        y_start = 50
        line_spacing = 16
    
        # 5) Draw each line of the cheat sheet
        for i, line in enumerate(HELP_TEXT):
            y_pos = y_start + i * line_spacing
            self.viewer.canvas.draw_text(
                x_start,
                y_pos,
                line,
                color=text_color,  # from HUD
                font_size=14       # pick size to taste
            )

    def toggle_grid(self):
        """Toggle the grid visibility in the current theme."""
        current_setting = CURRENT_THEME.get("show_grid", False)
        CURRENT_THEME["show_grid"] = not current_setting
        
        status = "enabled" if CURRENT_THEME["show_grid"] else "disabled"
        self.viewer.message_service.log_info(f"Grid {status}")
        self.viewer.redraw()
        
    def toggle_axes(self):
        """Toggle the coordinate axes visibility."""
        self.show_axes = not self.show_axes
        self.viewer.show_axes = self.show_axes  # Update viewer's property as well
        
        # Update molecule
        if self.show_axes:
            self.viewer.ortep_mol.create_coordinate_axes()
            self.viewer.message_service.log_info("Coordinate axes enabled")
        else:
            self.viewer.ortep_mol.toggle_coordinate_axes(False)
            self.viewer.message_service.log_info("Coordinate axes disabled")
        
        self.viewer.redraw()

    def toggle_help(self):
        """Toggle the help overlay display."""
        self.show_help = not self.show_help
        self.viewer.show_help = self.show_help  # Update viewer's property as well
        
        if self.show_help:
            # Print help text to console as well
            for line in HELP_TEXT:
                print(line)
        
        self.viewer.redraw()

    def toggle_hydrogens(self):
        """Toggle the display of hydrogen atoms."""
        self.show_hydrogens = not self.show_hydrogens
        self.viewer.show_hydrogens = self.show_hydrogens  # Update viewer's property as well
        
        status = "shown" if self.show_hydrogens else "hidden"
        self.viewer.message_service.log_info(f"Hydrogens {status}")
        self.viewer.set_frame(self.viewer.current_frame)

    def fit_molecule_to_window(self):
        """Fit the molecule to the window by adjusting scale and offset."""
        from geometry_utils import rotate_point
        xs, ys = [], []
        for atom in self.viewer.ortep_mol.atoms:
            x_rot, y_rot, _ = rotate_point(atom.x, atom.y, atom.z,
                                          self.viewer.view_params.rx,
                                          self.viewer.view_params.ry,
                                          self.viewer.view_params.rz)
            xs.append(x_rot)
            ys.append(y_rot)
        if not xs or not ys:
            return
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        margin = 20
        available_width = self.viewer.canvas.width - 2 * margin
        available_height = self.viewer.canvas.height - 2 * margin
        extent_x = max_x - min_x if max_x > min_x else 1.0
        extent_y = max_y - min_y if max_y > min_y else 1.0
        new_scale = min(available_width / extent_x, available_height / extent_y)
        if new_scale < self.viewer.view_params.scale:
            self.viewer.view_params.scale = new_scale
        center_x = (min_x + max_x) / 2.0
        center_y = (min_y + max_y) / 2.0
        self.viewer.view_params.x_offset = self.viewer.canvas.width / 2 - center_x * self.viewer.view_params.scale
        self.viewer.view_params.y_offset = self.viewer.canvas.height / 2 - center_y * self.viewer.view_params.scale
        self.viewer.message_service.log_info(f"Molecule fitted to window (scale: {new_scale:.1f})")
        self.viewer.redraw()

    def reset_view(self):
        """Reset the view to default parameters."""
        self.viewer.view_params.rx = 0.0
        self.viewer.view_params.ry = 0.0
        self.viewer.view_params.rz = 0.0
        self.viewer.view_params.scale = 100.0
        self.viewer.view_params.x_offset = self.viewer.canvas.width // 2
        self.viewer.view_params.y_offset = self.viewer.canvas.height // 2
        self.viewer.message_service.log_info("View reset to default")
        self.viewer.redraw()

    def hit_test(self, x, y):
        """Perform hit testing to determine what object is at the given coordinates."""
        z_objects = self.viewer.renderer.build_render_list(self.viewer.ortep_mol, self.viewer.view_params)
        for obj in sorted(z_objects, key=lambda o: o.z_value, reverse=True):
            if hasattr(obj, 'contains') and obj.contains(x, y):
                return obj
        return None


# ===============================================================
# Main Viewer Class
# ===============================================================

class MoleculeViewer(X11Window):
    """
    A viewer for ORTEP-style molecule visualization with support for interactive
    bond editing, frame navigation, normal mode visualization, and more.
    """
    def __init__(self, ortep_molecule, width=800, height=600,
                 ss_factor=1, tile_size=128):
        # --- Initialize Window and Canvas ---
        # Select the appropriate canvas class
        if ss_factor <= 1:
            canvas_class = X11CanvasBasic
        else:
            canvas_class = lambda w: X11CanvasSS(w, ss_factor=ss_factor,
                                                 tile_size=tile_size)
        
        # Get background color from theme settings
        background_color = CANVAS_SETTINGS["background_color"]
        
        super().__init__(
            width=width,
            height=height,
            title="ORTEP Style Molecule (Refactored)",
            canvas_class=canvas_class,
            background_color=background_color
        )
        
        # --- Initialize Core Attributes ---
        self.ortep_mol = ortep_molecule
        self.view_params = ViewParams(
            rx=0.0, ry=0.0, rz=0.0,
            scale=100.0,
            x_offset=width // 2,
            y_offset=height // 2
        )
        self.renderer = ORTEP_MoleculeRenderer()

        # --- Initialize Component Managers ---
        self._initialize_managers()
        
        # --- Initialize UI Components ---
        self._initialize_ui()
        
        # --- Initialize State Variables ---
        self._initialize_state()
        
        # --- Set up Coordinate Axes if Enabled ---
        if self.show_axes:
            self.ortep_mol.create_coordinate_axes()
            
        # Log initialization
        self.message_service.log_info(f"Viewer initialized ({width}x{height})")

    def _initialize_managers(self):
        """Initialize all manager components."""
        self._selection_manager = _SelectionManager()
        self._event_handler = _EventHandler(self)
        self._normal_mode_manager = _NormalModeManager(self)
        self._graph_manager = _GraphManager(self)
        self._view_manager = _ViewManager(self)
        self._traj_manager = None  # Set when trajectory is loaded
        self.bond_edit_tracker = BondEditTracker()

    def _initialize_ui(self):
        """Initialize UI components like HUD and message panels."""
        # Initialize HUD panel with themed styles
        self.hud_panel = HUDPanel(
            x=HUD_STYLE["x"],
            y_offset=HUD_STYLE["y_offset"],
            line_spacing=HUD_STYLE["line_spacing"],
            font_size=HUD_STYLE["font_size"],
            color=HUD_STYLE["color"]
        )
        
        # Initialize message service with themed settings
        self.message_service = MessageService(max_messages=MESSAGE_PANEL_STYLE.get("max_messages", 3))
        
        # Initialize message panel with themed styles
        self.message_panel = MessagePanel(
            self.message_service, 
            x=MESSAGE_PANEL_STYLE["x"], 
            y_offset=MESSAGE_PANEL_STYLE["y_offset"],
            line_spacing=MESSAGE_PANEL_STYLE["line_spacing"], 
            font_size=MESSAGE_PANEL_STYLE["font_size"],
            bg_color=MESSAGE_PANEL_STYLE["bg_color"],
            padding=MESSAGE_PANEL_STYLE["padding"]
        )

    def _initialize_state(self):
        """Initialize state variables for the viewer."""
        # Mouse state
        self.last_mouse_x = None
        self.last_mouse_y = None
        self.active_button = None

        # Selection-related members
        self.selected_atoms = []
        self.selected_bonds = []
        self.selected_atom_indices = []
        self.selected_bond_ids = []

        # Trajectory related
        self.trajectory = None
        self.current_frame = 0
        self.total_frames = 1
        
        # Display options (some are managed by _ViewManager but needed at main level too)
        self.show_hydrogens = self._view_manager.show_hydrogens
        self.show_axes = self._view_manager.show_axes
        self.show_help = self._view_manager.show_help
        
        # Expose graph properties for compatibility with event_handler
        self.energy_graph = None  # This will be updated when _graph_manager creates it
        self.graph_mode = "energy"  # Initialize to match _graph_manager

    # --- Public API Methods ---

    def set_trajectory(self, trajectory):
        """
        Set the trajectory and initialize trajectory-related components.
        
        Parameters:
            trajectory: The trajectory object to visualize
        """
        self.trajectory = trajectory
        self.total_frames = len(trajectory._raw_frames)
        self._traj_manager = _TrajectoryManager(trajectory)
        
        # Initialize normal mode manager with trajectory data
        self._normal_mode_manager.initialize_from_trajectory(trajectory)
        
        # Reset the energy graph so it will be recreated with the new trajectory data
        self._graph_manager.energy_graph = None
        self.energy_graph = None  # Keep viewer's property in sync
        
        # Log energy unit information
        energies, energy_info = trajectory.energy_trajectory(
            convert_if_hartrees=True,
            convert_to_unit=DEFAULT_ENERGY_UNIT
        )
        
        if energy_info['original_unit'] != 'unknown':
            original_unit = ENERGY_UNITS.get(
                energy_info['original_unit'], 
                {'name': energy_info['original_unit']}
            )['name']
            
            if energy_info['converted_unit'] != 'unknown':
                converted_unit = ENERGY_UNITS.get(
                    energy_info['converted_unit'], 
                    {'name': energy_info['converted_unit']}
                )['name']
                
                if energy_info['normalized']:
                    self.message_service.log_info(
                        f"Energies: {original_unit} converted to {converted_unit} (normalized to min=0)"
                    )
                else:
                    self.message_service.log_info(
                        f"Energies: {original_unit} converted to {converted_unit}"
                    )

    def set_frame(self, frame_index):
        """
        Set the current frame and update molecule visualization.
        
        Parameters:
            frame_index: Index of the frame to display
        """
        if self.trajectory is None:
            return
            
        # Clear normal mode vectors when changing frames
        if self._normal_mode_manager.show_normal_modes:
            self._normal_mode_manager.show_normal_modes = False
            self._normal_mode_manager.clear_normal_mode_vectors()
        
        # Delegate frame conversion to the trajectory manager with bond edits
        self.current_frame = frame_index
        self.ortep_mol = self._traj_manager.convert_frame(
            frame_index, 
            self.show_hydrogens,
            bond_edit_tracker=self.bond_edit_tracker
        )

        # Update the energy graph if it exists
        if self._graph_manager.energy_graph is not None:
            self._graph_manager.energy_graph.update_current_frame(frame_index)
            # Keep viewer's energy_graph property in sync for compatibility with event_handler
            self.energy_graph = self._graph_manager.energy_graph
            self.graph_mode = self._graph_manager.graph_mode

        # Check if the selected bond for the graph still exists
        self._graph_manager.check_selected_bond_exists()

        # Reapply persistent selection
        self._reapply_selections()

        # Recreate coordinate axes if they should be shown
        if self.show_axes:
            self.ortep_mol.create_coordinate_axes()
            
        # Show normal mode prompt if available for this frame
        if self._normal_mode_manager.has_normal_modes and frame_index == self._normal_mode_manager.normal_mode_frame_index:
            self.message_service.log_info(f"Normal modes available for this frame. Press 'v' to view.")

        self.redraw()

    def _reapply_selections(self):
        """Reapply stored selections to the current frame."""
        # Clear current selections
        self.selected_atoms = []
        self.selected_bonds = []
        
        # Reapply atom selections
        for atom in self.ortep_mol.atoms:
            if atom.index in self.selected_atom_indices:
                atom.selected = True
                self.selected_atoms.append(atom)
            else:
                atom.selected = False
                
        # Reapply bond selections
        for bond in self.ortep_mol.bonds:
            key = (min(bond.atom1.index, bond.atom2.index), max(bond.atom1.index, bond.atom2.index))
            if key in self.selected_bond_ids:
                bond.selected = True
                self.selected_bonds.append(bond)
            else:
                bond.selected = False

    def redraw(self):
        """
        Redraw the entire scene, including the molecule, grid, graph, and UI components.
        """
        # Check if window has been resized
        self._view_manager.check_resize()
    
        # Clear the drawing surface
        self.canvas.clear()
    
        # Draw grid if enabled
        self._view_manager.draw_grid()
    
        # Draw the main molecule (atoms, bonds, vectors)
        self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
    
        # Handle the energy/bond-length graph (thumbnail)
        self._graph_manager.ensure_energy_graph()
        if self._graph_manager.energy_graph is not None:
            self._graph_manager.energy_graph.draw_graph()
            # Keep viewer's properties synchronized
            self.energy_graph = self._graph_manager.energy_graph
            self.graph_mode = self._graph_manager.graph_mode
    
        # Update the HUD and message panel
        self._update_ui()
    
        # Show help overlay if needed
        if self.show_help:
            self._view_manager.draw_help_overlay()
    
        # Push all changes to the actual screen
        self.canvas.flush()

    def _update_ui(self):
        """Update the HUD and message panel."""
        # Update the HUD information
        self.update_info_message()
        self.hud_panel.draw(self.canvas)
        
        # Draw recent messages from the message panel
        self.message_panel.draw(self.canvas)

    def update_info_message(self):
        """Update the HUD with current info about selection, frame, energy, etc."""
        lines = []
        if self.selected_bonds:
            b = self.selected_bonds[0]
            bond_type = type(b).__name__.replace("Bond", "")
            from bond_manager import NON_REMOVAL_BOND_CYCLE, get_next_bond_type
            next_bond = get_next_bond_type(b, NON_REMOVAL_BOND_CYCLE)
            next_bond_type = next_bond.__name__.replace("Bond", "") if next_bond else "None"
            lines.append(f"Bond: {b.atom1.symbol}{b.atom1.index}-{b.atom2.symbol}{b.atom2.index} | Type: {bond_type} | Length: {b.length:.4f} Å")
            lines.append(f"Next bond type: {next_bond_type}")
            if len(self.selected_bonds) == 1:
                lines.append("Press 'p' to toggle bond length plot")
        elif self.selected_atoms:
            sel_info = ", ".join(f"{a.symbol}{a.index}" for a in self.selected_atoms)
            lines.append(f"Selected atoms: {sel_info}")
        else:
            lines.append("No object selected.")
            lines.append("Click to select an atom; Shift-click to multi-select.")
        lines.append(f"Zoom: {self.view_params.scale:.1f}")
        # Modify frame display to include coordinate/time metadata if available
        if self.trajectory:
            coords, coord_info = self.trajectory.coordinate_trajectory(skip_none=False)
            if coord_info.get('field') is not None and len(coords) > self.current_frame and not np.isnan(coords[self.current_frame]):
                field_label = "Coord" if coord_info['field'] == "coord" else "Time"
                lines.append(f"Frame: {self.current_frame} / {self.total_frames - 1}, {field_label}: {coords[self.current_frame]:.3f}")
            else:
                lines.append(f"Frame: {self.current_frame} / {self.total_frames - 1}")
        else:
            lines.append(f"Frame: {self.current_frame} / {self.total_frames - 1}")
        
        # Display energy information with units and method
        if self.trajectory:
            energies, energy_info = self.trajectory.energy_trajectory(
                convert_if_hartrees=True,
                convert_to_unit=DEFAULT_ENERGY_UNIT
            )
            
            if len(energies) > self.current_frame and not np.isnan(energies[self.current_frame]):
                energy = energies[self.current_frame]
                
                # Get energy unit symbol
                unit_symbol = ENERGY_UNITS.get(
                    energy_info['converted_unit'], 
                    {'symbol': energy_info['converted_unit']}
                )['symbol']
                
                # Get energy method if available
                energy_method = ""
                if ('energy_data' in self.trajectory.metadata and 
                    self.current_frame in self.trajectory.metadata['energy_data']):
                    method = self.trajectory.metadata['energy_data'][self.current_frame].get('type', '')
                    if method:
                        energy_method = f" ({method})"
                
                # Construct energy display string
                if energy_info['normalized']:
                    lines.append(f"Energy: {energy:.2f} {unit_symbol}{energy_method} (rel. to min)")
                else:
                    lines.append(f"Energy: {energy:.2f} {unit_symbol}{energy_method}")
        
        # Add normal mode information
        lines = self._normal_mode_manager.update_info_message(lines)
        
        # Add display options status
        display_options = []
        display_options.append("H" if self.show_hydrogens else "no H")
        display_options.append("axes" if self.show_axes else "no axes")
        lines.append(f"Display: {', '.join(display_options)}")
        
        # Add graph mode info
        lines = self._graph_manager.update_info_message(lines)
        
        # Add key bindings help
        lines.append("Press '?' for help.")
            
        self.hud_panel.update_lines(lines)

    def fit_molecule_to_window(self):
        """Fit the molecule to the window."""
        self._view_manager.fit_molecule_to_window()

    def reset_view(self):
        """Reset the view to default parameters."""
        self._view_manager.reset_view()

    def toggle_graph_mode(self):
        """Toggle between energy and bond length visualization modes."""
        self._graph_manager.toggle_graph_mode()
        # Update viewer's properties for compatibility
        self.energy_graph = self._graph_manager.energy_graph
        self.graph_mode = self._graph_manager.graph_mode

    def toggle_normal_modes(self):
        """Toggle display of normal mode vectors."""
        self._normal_mode_manager.toggle_normal_modes()

    def cycle_normal_mode(self, direction):
        """Cycle to the next/previous normal mode."""
        self._normal_mode_manager.cycle_normal_mode(direction)

    def adjust_normal_mode_scale(self, factor):
        """Adjust the scale factor for normal mode vectors."""
        self._normal_mode_manager.adjust_normal_mode_scale(factor)

    def toggle_grid(self):
        """Toggle the grid visibility."""
        self._view_manager.toggle_grid()

    def toggle_axes(self):
        """Toggle the coordinate axes visibility."""
        self._view_manager.toggle_axes()

    def toggle_hydrogens(self):
        """Toggle the display of hydrogen atoms."""
        self._view_manager.toggle_hydrogens()

    def clear_bond_edits(self):
        """Clear all bond edits and refresh the current frame."""
        self.bond_edit_tracker.clear_edits()
        self.message_service.log_info("Cleared all bond edits")
        self.set_frame(self.current_frame)

    def hit_test(self, x, y):
        """Perform hit testing to see what's under the cursor."""
        return self._view_manager.hit_test(x, y)

    def dump_xyz(self):
        """
        Export the current frame to an XYZ file, with dynamic filename based on
        the loaded trajectory file and current frame.
        """
        if not self.trajectory:
            self.message_service.log_info("No trajectory data available to export")
            return
    
        # Generate filename based on trajectory metadata
        base_name = self.trajectory.metadata.get('file_name', 'quickORTEP')
    
        # For multi-frame files, include the frame number
        if self.total_frames > 1:
            # Zero-pad the frame number to 3 digits
            frame_num = str(self.current_frame).zfill(3)
            filename = f"{base_name}_{frame_num}.xyz"
        else:
            filename = f"{base_name}.xyz"
    
        # Get the original molecule from the trajectory
        original_mol = self.trajectory.get_frame(self.current_frame)
    
        # Export using the molecule's to_xyz method
        export_xyz(original_mol, filename, self.message_service)

    def dump_svg(self):
        """
        Export the current view to an SVG file, with dynamic filename based on
        the loaded trajectory file and current frame.
        """
        # Generate filename based on trajectory metadata
        if self.trajectory and hasattr(self.trajectory, 'metadata'):
            base_name = self.trajectory.metadata.get('file_name', 'quickORTEP')
            
            # For multi-frame files, include the frame number
            if self.total_frames > 1:
                # Zero-pad the frame number to 3 digits
                frame_num = str(self.current_frame).zfill(3)
                filename = f"{base_name}_{frame_num}.svg"
            else:
                filename = f"{base_name}.svg"
        else:
            filename = "quickORTEP_export.svg"
        
        export_svg(self.canvas, self.renderer, self.ortep_mol, self.view_params, 
                   filename, self.message_service)
    
    # --- Event Handling Methods (delegated to EventHandler) ---
    
    def handle_key(self, evt):
        self._event_handler.handle_key(evt)

    def handle_button_press(self, evt):
        self._event_handler.handle_button_press(evt)

    def handle_motion(self, evt):
        self._event_handler.handle_motion(evt)

    def handle_button_release(self, evt):
        self._event_handler.handle_button_release(evt)
