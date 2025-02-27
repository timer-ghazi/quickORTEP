#!/usr/bin/env python3
"""
ortep_viewer.py

The MoleculeViewer class manages an ORTEP_Molecule, its rendering,
and handles user events. It now delegates trajectory conversion,
selection, event handling, and export functionality to helper modules.
"""

import sys
import time
import threading
from Xlib import XK, X
from x11view.window import X11Window, X11CanvasBasic, X11CanvasSS
from hud import HUDPanel
from ortep_renderer import ORTEP_MoleculeRenderer
from config import VIEWER_INTERACTION, MINIMAL_THEME
from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from message_service import MessageService
from message_panel import MessagePanel
from graph_viewer import GraphViewer

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
        super().__init__(
            width=width,
            height=height,
            title="ORTEP Style Molecule (Refactored)",
            canvas_class=canvas_class
        )

        self.ortep_mol = ortep_molecule
        self.view_params = ViewParams(
            rx=0.0, ry=0.0, rz=0.0,
            scale=100.0,
            x_offset=width // 2,
            y_offset=height // 2
        )
        self.renderer = ORTEP_MoleculeRenderer()

        self.hud_panel = HUDPanel(x=10, y_offset=75, line_spacing=16, font_size=14, color=(0, 0, 0))
        self.message_service = MessageService(max_messages=3)
        self.message_panel = MessagePanel(
            self.message_service, 
            x=10, 
            y_offset=10,
            line_spacing=16, 
            font_size=12
        )
        self.message_service.log_info(f"Viewer initialized ({width}x{height})")

        # Drawing lock to ensure thread safety.
        self.draw_lock = threading.Lock()

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

        self.show_hydrogens = True

        # Initialize helper modules.
        self._selection_manager = _SelectionManager()
        self._event_handler = _EventHandler(self)
        self._traj_manager = None  # Will be set when a trajectory is assigned.

    def set_trajectory(self, trajectory):
        """
        Set the trajectory and initialize the trajectory manager.
        """
        self.trajectory = trajectory
        self.total_frames = len(trajectory._raw_frames)
        self._traj_manager = _TrajectoryManager(trajectory)
        
        # Reset the energy graph so it will be recreated with the new trajectory data
        self.energy_graph = None

    def ensure_energy_graph(self):
        """
        Initialize the energy vs. frame graph if a trajectory with energy data is available.
        """
        # Skip if graph already exists or no trajectory is available
        if self.energy_graph is not None or self.trajectory is None:
            return
        
        # Check for energy data
        energies = self.trajectory._frame_energies
        if not energies or all(e is None for e in energies):
            return
        
        # Prepare graph data
        x_values = list(range(len(energies)))
        y_values = [e if e is not None else 0.0 for e in energies]
        
        # Configure graph position in bottom right corner
        thumb_width = 150
        thumb_height = 80
        thumb_margin = 20
        thumb_x = self.canvas.width - thumb_width - thumb_margin
        thumb_y = self.canvas.height - thumb_height - thumb_margin
        
        # Create a custom minimal theme
        custom_minimal_theme = MINIMAL_THEME.copy()
        custom_minimal_theme["margin"] = {"left": 10, "right": 10, "top": 10, "bottom": 10}
        
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
            y_axis_title="",
            title="Energy",
            theme=custom_minimal_theme
        )

    def set_frame(self, frame_index):
        if self.trajectory is None:
            return
        # Delegate frame conversion to the trajectory manager.
        self.current_frame = frame_index
        self.ortep_mol = self._traj_manager.convert_frame(frame_index, self.show_hydrogens)

        # Update the energy graph if it exists
        if self.energy_graph is not None:
            self.energy_graph.update_current_frame(frame_index)

        # Reapply persistent selection.
        for atom in self.ortep_mol.atoms:
            atom.selected = (atom.index in self.selected_atom_indices)
        self.selected_bonds = []
        for bond in self.ortep_mol.bonds:
            key = (min(bond.atom1.index, bond.atom2.index), max(bond.atom1.index, bond.atom2.index))
            if key in self.selected_bond_ids:
                bond.selected = True
                self.selected_bonds.append(bond)
            else:
                bond.selected = False

        self.redraw()

    def update_info_message(self):
        lines = []
        if self.selected_bonds:
            b = self.selected_bonds[0]
            bond_type = type(b).__name__.replace("Bond", "")
            from bond_manager import NON_REMOVAL_BOND_CYCLE, get_next_bond_type
            next_bond = get_next_bond_type(b, NON_REMOVAL_BOND_CYCLE)
            next_bond_type = next_bond.__name__.replace("Bond", "") if next_bond else "None"
            lines.append(f"Bond: {b.atom1.symbol}{b.atom1.index}-{b.atom2.symbol}{b.atom2.index} | Type: {bond_type} | Length: {b.length:.4f} Å")
            lines.append(f"Next bond type: {next_bond_type}")
        elif self.selected_atoms:
            sel_info = ", ".join(f"{a.symbol}{a.index}" for a in self.selected_atoms)
            lines.append(f"Selected atoms: {sel_info}")
        else:
            lines.append("No object selected.")
            lines.append("Click to select an atom; Shift-click to multi-select.")
        lines.append(f"Zoom: {self.view_params.scale:.1f}")
        lines.append(f"Frame: {self.current_frame} / {self.total_frames - 1}")
        if self.trajectory and self.trajectory._frame_energies[self.current_frame] is not None:
            energy = self.trajectory._frame_energies[self.current_frame]
            lines.append(f"Energy: {energy:.4f} a.u.")
        lines.append(f"Hydrogens: {'shown' if self.show_hydrogens else 'hidden'}")
        self.hud_panel.update_lines(lines)

    def redraw(self):
        with self.draw_lock:
            self.canvas.clear()
            self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
            
            # Ensure the energy graph is initialized if needed
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
        filename = "quickORTEP_export.svg"
        export_svg(self.canvas, self.renderer, self.ortep_mol, self.view_params, filename)
        self.message_service.log_info(f"SVG export saved to {filename}")

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