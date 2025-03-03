#!/usr/bin/env python3
"""
ortep_viewer.py

The MoleculeViewer class manages an ORTEP_Molecule, its rendering,
and handles user events. It now delegates trajectory conversion,
selection, event handling, and export functionality to helper modules.
It also supports bond propagation across trajectory frames.
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

# Import our helper modules.
from selection_manager import _SelectionManager
from trajectory_manager import _TrajectoryManager
from event_handler import _EventHandler
from export import export_svg

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

class MoleculeViewer(X11Window):
    def __init__(self, ortep_molecule, width=800, height=600,
                 ss_factor=1, tile_size=128):
        # Select the appropriate canvas class.
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
            background_color=background_color  # Pass the themed background color
        )
        
        # Store initial canvas dimensions to detect resizing
        self._last_canvas_width = width
        self._last_canvas_height = height

        self.ortep_mol = ortep_molecule
        self.view_params = ViewParams(
            rx=0.0, ry=0.0, rz=0.0,
            scale=100.0,
            x_offset=width // 2,
            y_offset=height // 2
        )
        self.renderer = ORTEP_MoleculeRenderer()

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
        
        self.message_service.log_info(f"Viewer initialized ({width}x{height})")


        # Mouse state.
        self.last_mouse_x = None
        self.last_mouse_y = None
        self.active_button = None

        # Selection-related members.
        self.selected_atoms = []
        self.selected_bonds = []
        self.selected_atom_indices = []
        self.selected_bond_ids = []

        # Trajectory related.
        self.trajectory = None
        self.current_frame = 0
        self.total_frames = 1
        
        # Energy graph for trajectory visualization
        self.energy_graph = None
        
        # Graph mode tracking
        self.graph_mode = "energy"  # Options: "energy", "bond_length"
        self.selected_bond_for_graph = None  # Will store the bond ID for graph (atom1_idx, atom2_idx)

        self.show_hydrogens = True

        # Initialize helper modules.
        self._selection_manager = _SelectionManager()
        self._event_handler = _EventHandler(self)
        self._traj_manager = None  # Will be set when a trajectory is assigned.

        # Initialize bond edit tracker for propagating bond changes across frames
        self.bond_edit_tracker = BondEditTracker()

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
        scale = self.view_params.scale
        x_offset = self.view_params.x_offset
        y_offset = self.view_params.y_offset
        
        # Calculate visible area in Ångströms
        left_a = -x_offset / scale
        right_a = left_a + (self.canvas.width / scale)
        top_a = -y_offset / scale
        bottom_a = top_a + (self.canvas.height / scale)
        
        # Round to the nearest major grid line outside the viewport
        left_a_grid = math.floor(left_a / major_spacing) * major_spacing
        right_a_grid = math.ceil(right_a / major_spacing) * major_spacing
        top_a_grid = math.floor(top_a / major_spacing) * major_spacing
        bottom_a_grid = math.ceil(bottom_a / major_spacing) * major_spacing
        
        # Draw major grid lines (vertical)
        x_a = left_a_grid
        while x_a <= right_a_grid:
            x_px = int(x_a * scale + x_offset)
            self.canvas.draw_line(x_px, 0, x_px, self.canvas.height, 
                                thickness=major_thickness, color=major_color)
            x_a += major_spacing
        
        # Draw major grid lines (horizontal)
        y_a = top_a_grid
        while y_a <= bottom_a_grid:
            y_px = int(y_a * scale + y_offset)
            self.canvas.draw_line(0, y_px, self.canvas.width, y_px, 
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
                    self.canvas.draw_dashed_line(minor_x_px, 0, minor_x_px, self.canvas.height, 
                                            thickness=minor_thickness, color=minor_color)
                x_a += major_spacing
            
            # Draw minor grid lines (horizontal)
            y_a = top_a_grid
            while y_a <= bottom_a_grid:
                for i in range(1, minor_divisions):
                    minor_y_a = y_a + (i * minor_spacing)
                    minor_y_px = int(minor_y_a * scale + y_offset)
                    self.canvas.draw_dashed_line(0, minor_y_px, self.canvas.width, minor_y_px, 
                                            thickness=minor_thickness, color=minor_color)
                y_a += major_spacing

    def set_trajectory(self, trajectory):
        """
        Set the trajectory and initialize the trajectory manager.
        """
        self.trajectory = trajectory
        self.total_frames = len(trajectory._raw_frames)
        self._traj_manager = _TrajectoryManager(trajectory)
        
        # Reset the energy graph so it will be recreated with the new trajectory data
        self.energy_graph = None
        
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

    def calculate_bond_length_trajectory(self, atom1_idx, atom2_idx):
        """
        Calculate bond length between two atoms across all frames in the trajectory.
        
        Parameters:
            atom1_idx, atom2_idx (int): Atom indices
            
        Returns:
            tuple: (x_values, y_values) where x_values are frame numbers and 
                   y_values are bond lengths in Ångströms
        """
        if self.trajectory is None:
            return [], []
        
        x_values = list(range(self.total_frames))
        y_values = []
        
        for frame_idx in range(self.total_frames):
            # Get the molecule for this frame (without changing the current frame)
            mol = self._traj_manager.trajectory.get_frame(frame_idx)
            # Calculate distance between the atoms (0-indexed in molecule_nci)
            distance = mol.distance(atom1_idx-1, atom2_idx-1)
            y_values.append(distance)
        
        return x_values, y_values

    def ensure_energy_graph(self):
        """
        Initialize or update the graph based on the current mode (energy or bond length).
        """
        # Configure graph position in bottom right corner
        thumb_width = GRAPH_SETTINGS["thumbnail_width"]
        thumb_height = GRAPH_SETTINGS["thumbnail_height"]
        thumb_margin = GRAPH_SETTINGS["thumbnail_margin"]
        thumb_x = self.canvas.width - thumb_width - thumb_margin
        thumb_y = self.canvas.height - thumb_height - thumb_margin
        
        # If we're in bond length mode, ensure we have a selected bond
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph is None:
            # If no bond is selected, revert to energy mode
            self.graph_mode = "energy"
        
        # Prepare data based on graph mode
        if self.graph_mode == "energy":
            # Skip if no trajectory is available
            if self.trajectory is None:
                return
            
            # Get energy data using the energy_trajectory method with unit conversion
            energies, energy_info = self.trajectory.energy_trajectory(
                convert_if_hartrees=True,
                convert_to_unit=DEFAULT_ENERGY_UNIT
            )
            
            if len(energies) == 0 or np.isnan(energies).all():
                return
            
            # Prepare graph data
            x_values = list(range(len(energies)))
            y_values = [e if not np.isnan(e) else 0.0 for e in energies]
            
            # Get unit symbol for the title
            unit_symbol = ENERGY_UNITS.get(
                energy_info['converted_unit'], 
                {'symbol': energy_info['converted_unit']}
            )['symbol']
            
            title = "Energy"
            # Create title with unit info
            #---if energy_info['normalized']:
            #---    title = f"Energyunit_symbol})"
            #---else:
            #---    title = f"Energy ({unit_symbol})"
                
            y_axis_title = ""
        else:
            # Bond length graph
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            atom1 = next((a for a in self.ortep_mol.atoms if a.index == atom1_idx), None)
            atom2 = next((a for a in self.ortep_mol.atoms if a.index == atom2_idx), None)
            
            # Get bond length data
            x_values, y_values = self.calculate_bond_length_trajectory(atom1_idx, atom2_idx)
            
            # Prepare title with atom symbols
            if atom1 and atom2:
                title = f"{atom1.symbol}{atom1_idx}-{atom2.symbol}{atom2_idx}"
            else:
                title = "Bond Length"
            y_axis_title = ""
        
        # Create a custom minimal theme
        custom_minimal_theme = MINIMAL_THEME.copy()
        custom_minimal_theme["margin"] = {"left": 10, "right": 10, "top": 10, "bottom": 10}
        
        if self.energy_graph is None:
            # Create the graph
            self.energy_graph = GraphViewer(
                canvas=self.canvas,
                xdata=x_values,
                ydata=y_values,
                mode="minimal",
                current_frame=self.current_frame,
                region_x=thumb_x,
                region_y=thumb_y,
                region_width=thumb_width,
                region_height=thumb_height,
                x_axis_title="",
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
                    y_axis_title=y_axis_title
                )

    def set_frame(self, frame_index):
        """
        Set the current frame and apply bond edits if bond propagation is enabled.
        """
        if self.trajectory is None:
            return
            
        # Delegate frame conversion to the trajectory manager with bond edits
        self.current_frame = frame_index
        self.ortep_mol = self._traj_manager.convert_frame(
            frame_index, 
            self.show_hydrogens,
            bond_edit_tracker=self.bond_edit_tracker
        )

        # Update the energy graph if it exists
        if self.energy_graph is not None:
            self.energy_graph.update_current_frame(frame_index)

        # Reapply persistent selection.
        self.selected_atoms = []
        for atom in self.ortep_mol.atoms:
            if atom.index in self.selected_atom_indices:
                atom.selected = True
                self.selected_atoms.append(atom)
            else:
                atom.selected = False
                
        self.selected_bonds = []
        for bond in self.ortep_mol.bonds:
            key = (min(bond.atom1.index, bond.atom2.index), max(bond.atom1.index, bond.atom2.index))
            if key in self.selected_bond_ids:
                bond.selected = True
                self.selected_bonds.append(bond)
            else:
                bond.selected = False

        # If in bond length mode, check if the selected bond still exists
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph:
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            bond_exists = False
            for bond in self.ortep_mol.bonds:
                if ((bond.atom1.index == atom1_idx and bond.atom2.index == atom2_idx) or
                    (bond.atom1.index == atom2_idx and bond.atom2.index == atom1_idx)):
                    bond_exists = True
                    break
            
            if not bond_exists:
                # Bond doesn't exist in this frame, switch back to energy mode
                self.graph_mode = "energy"
                self.selected_bond_for_graph = None
                self.energy_graph = None  # Force recreation

        self.redraw()

    def update_info_message(self):
        """
        Update the HUD with current info, including bond propagation status and energy units.
        """
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
        
        lines.append(f"Hydrogens: {'shown' if self.show_hydrogens else 'hidden'}")
        
        #---- # Add bond propagation status
        #---- prop_status = "enabled" if self.bond_edit_tracker.enabled else "disabled"
        #---- edit_count = self.bond_edit_tracker.get_edit_count()
        #---- if edit_count > 0:
        #----     lines.append(f"Bond propagation: {prop_status} ({edit_count} edits)")
        #---- else:
        #----     lines.append(f"Bond propagation: {prop_status}")
            
        # Add graph mode info
        if self.graph_mode == "bond_length" and self.selected_bond_for_graph:
            atom1_idx, atom2_idx = self.selected_bond_for_graph
            atom1 = next((a for a in self.ortep_mol.atoms if a.index == atom1_idx), None)
            atom2 = next((a for a in self.ortep_mol.atoms if a.index == atom2_idx), None)
            if atom1 and atom2:
                lines.append(f"Graph: Bond length {atom1.symbol}{atom1_idx}-{atom2.symbol}{atom2_idx}")
        else:
            # Add unit info to the graph description
            if self.trajectory:
                _, energy_info = self.trajectory.energy_trajectory(
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
            
        self.hud_panel.update_lines(lines)

    def toggle_graph_mode(self):
        """
        Toggle between energy and bond length visualization modes.
        """
        if self.selected_bonds and len(self.selected_bonds) == 1:
            bond = self.selected_bonds[0]
            atom1_idx = bond.atom1.index
            atom2_idx = bond.atom2.index
            
            if self.graph_mode == "energy":
                # Switch to bond length mode
                self.graph_mode = "bond_length"
                self.selected_bond_for_graph = (atom1_idx, atom2_idx)
                self.message_service.log_info(
                    f"Showing bond length for {bond.atom1.symbol}{atom1_idx}-{bond.atom2.symbol}{atom2_idx}"
                )
            else:
                # Switch back to energy mode
                self.graph_mode = "energy"
                self.message_service.log_info("Showing energy plot")
            
            # Refresh the graph
            self.energy_graph = None  # Force recreation
            self.ensure_energy_graph()
        elif self.graph_mode == "bond_length":
            # Switch back to energy mode if in bond length mode but no bond selected
            self.graph_mode = "energy"
            self.message_service.log_info("Showing energy plot")
            self.energy_graph = None  # Force recreation
            self.ensure_energy_graph()
        else:
            self.message_service.log_info("Select a bond first to toggle length graph")

    def clear_bond_edits(self):
        """
        Clear all bond edits and refresh the current frame.
        """
        self.bond_edit_tracker.clear_edits()
        self.message_service.log_info("Cleared all bond edits")
        self.set_frame(self.current_frame)

    def toggle_grid(self):
        """
        Toggle the grid visibility in the current theme.
        """
        current_setting = CURRENT_THEME.get("show_grid", False)
        CURRENT_THEME["show_grid"] = not current_setting
        
        status = "enabled" if CURRENT_THEME["show_grid"] else "disabled"
        self.message_service.log_info(f"Grid {status}")
        self.redraw()

    def redraw(self):
        # Check if window has been resized
        if self._last_canvas_width != self.canvas.width or self._last_canvas_height != self.canvas.height:
            # Window was resized, update stored dimensions
            self._last_canvas_width = self.canvas.width
            self._last_canvas_height = self.canvas.height
            # Make sure graph position is updated on next ensure_energy_graph call
            if self.energy_graph is not None:
                # We'll update the position when ensure_energy_graph is called
                pass
        
        self.canvas.clear()
        
        # Draw the grid if enabled in the current theme
        self.draw_grid()
        
        self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
        
        # Ensure the energy graph is initialized if needed and positioned correctly
        self.ensure_energy_graph()
        
        # Draw the energy graph if it exists
        if self.energy_graph is not None:
            self.energy_graph.draw_graph()
            
        self.update_info_message()
        self.hud_panel.draw(self.canvas)
        self.message_panel.draw(self.canvas)
        self.canvas.flush()

    def fit_molecule_to_window(self):
        from geometry_utils import rotate_point
        xs, ys = [], []
        for atom in self.ortep_mol.atoms:
            x_rot, y_rot, _ = rotate_point(atom.x, atom.y, atom.z,
                                           self.view_params.rx,
                                           self.view_params.ry,
                                           self.view_params.rz)
            xs.append(x_rot)
            ys.append(y_rot)
        if not xs or not ys:
            return
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        margin = 20
        available_width = self.canvas.width - 2 * margin
        available_height = self.canvas.height - 2 * margin
        extent_x = max_x - min_x if max_x > min_x else 1.0
        extent_y = max_y - min_y if max_y > min_y else 1.0
        new_scale = min(available_width / extent_x, available_height / extent_y)
        if new_scale < self.view_params.scale:
            self.view_params.scale = new_scale
        center_x = (min_x + max_x) / 2.0
        center_y = (min_y + max_y) / 2.0
        self.view_params.x_offset = self.canvas.width / 2 - center_x * self.view_params.scale
        self.view_params.y_offset = self.canvas.height / 2 - center_y * self.view_params.scale
        self.message_service.log_info(f"Molecule fitted to window (scale: {new_scale:.1f})")
        self.redraw()

    def hit_test(self, x, y):
        z_objects = self.renderer.build_render_list(self.ortep_mol, self.view_params)
        for obj in sorted(z_objects, key=lambda o: o.z_value, reverse=True):
            if hasattr(obj, 'contains') and obj.contains(x, y):
                return obj
        return None

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

    def reset_view(self):
        self.view_params.rx = 0.0
        self.view_params.ry = 0.0
        self.view_params.rz = 0.0
        self.view_params.scale = 100.0
        self.view_params.x_offset = self.canvas.width // 2
        self.view_params.y_offset = self.canvas.height // 2
        self.message_service.log_info("View reset to default")
        self.redraw()

    # Delegate event handling to the _EventHandler.
    def handle_key(self, evt):
        self._event_handler.handle_key(evt)

    def handle_button_press(self, evt):
        self._event_handler.handle_button_press(evt)

    def handle_motion(self, evt):
        self._event_handler.handle_motion(evt)

    def handle_button_release(self, evt):
        self._event_handler.handle_button_release(evt)
