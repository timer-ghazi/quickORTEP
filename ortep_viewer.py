#!/usr/bin/env python3
"""
ortep_viewer.py

The viewer class that manages an ORTEP_Molecule, its renderer, and handles user events.
Now updated to support:
  - Uppercase 'B' toggles a bond between two selected atoms.
  - Lowercase 'b' cycles through bond types (without removal) for a selected bond.
  - A multi-line HUD panel displaying current state information.
  - Trajectory navigation via keys: '[' / ']' to move frame-by-frame, '{' / '}' to jump to the start/end,
    '=' and '-' to jump to the highest/lowest energy frame.
  - Toggle display of hydrogen atoms attached to carbons via key 'd' (persistent across frames).
  - Persistent selection of atoms and bonds across frames.
  - A drawing lock to avoid thread-safety issues with Xlib.
  - A scrolling message panel at the bottom of the window for status updates.
"""

import sys
import time
import threading
from Xlib import XK, X
from x11view.window import X11Window, X11CanvasBasic, X11CanvasSS
from x11view.svg_canvas import SVGCanvas
from ortep_renderer import ORTEP_MoleculeRenderer
from geometry_utils import rotate_point
from config import VIEWER_INTERACTION
from bond_manager import (
    cycle_existing_bond,
    cycle_atom_pair,
    toggle_bond,
)
from hud import HUDPanel
from ortep_molecule import ORTEP_Molecule, ORTEP_Atom
from bonds import CovalentBond, NCIBond
from message_service import MessageService
from message_panel import MessagePanel

class ViewParams:
    """
    Holds rotation angles, scale, and screen offset for the viewer.
    """
    def __init__(self, rx=0.0, ry=0.0, rz=0.0,
                 scale=100.0, x_offset=400, y_offset=300):
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.scale = scale
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.as_float = False   # default

class MoleculeViewer(X11Window):
    """
    The controller/view class that manages an ORTEP_Molecule, its renderer, and handles user events.
    
    - Uppercase 'B' toggles a bond: if two atoms are selected, creates or removes the bond.
    - Lowercase 'b' cycles through bond types (without removal) for a selected bond.
    - Displays a multi-line HUD panel with current state information.
    - Trajectory Navigation: Supports frame-by-frame navigation and energy-based jumps.
    - Toggle display of hydrogen atoms attached to carbons (key 'd').
    - Persistent Selection: Atom and bond selections persist across trajectory frames.
    - Uses a drawing lock to prevent thread-safety issues with Xlib.
    - Displays status messages in a scrolling panel at the bottom of the window.
    """
    def __init__(self, ortep_molecule, width=800, height=600,
                 ss_factor=1, tile_size=128):
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

        # Instantiate the HUDPanel for displaying multi-line info.
        self.hud_panel = HUDPanel(x=10, y_offset=75, line_spacing=16, font_size=14, color=(0, 0, 0))
        
        # Initialize message service and panel
        self.message_service = MessageService(max_messages=3)
        self.message_panel = MessagePanel(
            self.message_service, 
            x=10, 
            y_offset=10,  # Increase bottom margin
            line_spacing=16, 
            font_size=12
        )
        
        # Log initialization message
        self.message_service.log_info(f"Viewer initialized ({width}x{height})")

        # Drawing lock to prevent concurrent Xlib calls.
        self.draw_lock = threading.Lock()

        # Mouse state variables.
        self.active_button = None         # 'left', 'middle', 'right', or 'shift-left'
        self.last_mouse_x = None
        self.last_mouse_y = None
        self.last_click_time = 0.0

        # Lists for storing selected objects.
        self.selected_atoms = []   # List of ORTEP_Atom objects.
        self.selected_bonds = []   # List of underlying bond objects.
        self.click_start_x = None
        self.click_start_y = None

        # ===== Trajectory and Persistent Selection =====
        self.trajectory = None         # Will hold the Trajectory object
        self.current_frame = 0         # Current frame index in the trajectory
        self.total_frames = 1          # Total frames; updated when trajectory is set

        # Persistent selection storage (by unique atom indices and bond key tuples)
        self.selected_atom_indices = []      # e.g. [1, 3, ...]
        self.selected_bond_ids = []          # e.g. [(1,2), (2,3), ...]

        # Hydrogen display toggle: True = show hydrogens; False = hide hydrogens attached to carbons.
        self.show_hydrogens = True

    def set_frame(self, frame_index):
        """
        Update the viewer to display a new frame from the trajectory.
        Also, reapply persistent selection based on stored indices.
        
        Here we convert the molecule from the trajectory (a MoleculeWithNCIs)
        into an ORTEP_Molecule so that the atoms have the necessary attributes.
        Additionally, if hydrogens are to be hidden, filter out hydrogen atoms
        attached to carbons (and their corresponding bonds), then reassign atom indices.
        """
        if self.trajectory is None:
            return

        # Clamp frame_index to valid range.
        frame_index = max(0, min(frame_index, self.total_frames - 1))
        self.current_frame = frame_index
        
        # Log frame change
        energy = None
        if self.trajectory and self.trajectory._frame_energies[self.current_frame] is not None:
            energy = self.trajectory._frame_energies[self.current_frame]
       #     self.message_service.log_info(f"Frame {self.current_frame}/{self.total_frames-1} (E={energy:.4f} a.u.)")
       # else:
       #     self.message_service.log_info(f"Frame {self.current_frame}/{self.total_frames-1}")
        
        # Load the new molecule from the trajectory (MoleculeWithNCIs).
        new_mol_nci = self.trajectory.get_frame(frame_index)
        
        # Convert to ORTEP_Molecule.
        new_ortep_mol = ORTEP_Molecule()
        # Create ORTEP_Atom objects for each atom.
        for atom in new_mol_nci.atoms:
            a = ORTEP_Atom(atom.symbol, atom.x, atom.y, atom.z)
            new_ortep_mol.add_atom(a)
        
        # Recreate bonds from the bond matrix if available.
        n_atoms = len(new_mol_nci.atoms)
        bondmat = getattr(new_mol_nci, "bond_matrix", None)
        if bondmat is not None:
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    if bondmat[i, j] == 1:
                        b = CovalentBond(new_ortep_mol.atoms[i], new_ortep_mol.atoms[j])
                        new_ortep_mol.add_bond(b)
        
        # Add NCIs.
        for (i, j), interactions in new_mol_nci.ncis.items():
            nci_bond = NCIBond(new_ortep_mol.atoms[i], new_ortep_mol.atoms[j])
            new_ortep_mol.add_bond(nci_bond)
        
        # --- Hydrogen Filtering ---
        # If hydrogen atoms attached to carbons should be hidden, remove them and associated bonds.
        if not self.show_hydrogens:
            atoms_to_remove = []
            bonds_to_remove = []
            for bond in new_ortep_mol.bonds:
                if (bond.atom1.symbol == "H" and bond.atom2.symbol == "C"):
                    if bond.atom1 not in atoms_to_remove:
                        atoms_to_remove.append(bond.atom1)
                    bonds_to_remove.append(bond)
                elif (bond.atom2.symbol == "H" and bond.atom1.symbol == "C"):
                    if bond.atom2 not in atoms_to_remove:
                        atoms_to_remove.append(bond.atom2)
                    bonds_to_remove.append(bond)
            for bond in bonds_to_remove:
                if bond in new_ortep_mol.bonds:
                    new_ortep_mol.bonds.remove(bond)
            for atom in atoms_to_remove:
                if atom in new_ortep_mol.atoms:
                    new_ortep_mol.atoms.remove(atom)
            # Reassign atom indices after filtering.
            for idx, atom in enumerate(new_ortep_mol.atoms, start=1):
                atom.index = idx
        
        # Set the current molecule.
        self.ortep_mol = new_ortep_mol

        # Reapply persistent atom selections.
        for atom in self.ortep_mol.atoms:
            if atom.index in self.selected_atom_indices:
                atom.selected = True
            else:
                atom.selected = False

        # Reapply persistent bond selections and update the selected bonds list.
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
        """
        Build the list of HUD lines based on the current state and update the HUDPanel.
        """
        lines = []
        # Persistent selection info.
        if self.selected_bonds:
            b = self.selected_bonds[0]
            bond_type = type(b).__name__.replace("Bond", "")
            # Determine the next bond type in the cycle.
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
        
        # Append view parameters.
      #  lines.append(f"View: rx={self.view_params.rx:.1f}°, ry={self.view_params.ry:.1f}°, rz={self.view_params.rz:.1f}°")
        lines.append(f"Zoom: {self.view_params.scale:.1f}")
        
        # Trajectory Info.
        lines.append(f"Frame: {self.current_frame} / {self.total_frames - 1}")
        if self.trajectory and self.trajectory._frame_energies[self.current_frame] is not None:
            energy = self.trajectory._frame_energies[self.current_frame]
            lines.append(f"Energy: {energy:.4f} a.u.")
        
        # Hydrogen display status.
        lines.append(f"Hydrogens: {'shown' if self.show_hydrogens else 'hidden'}")
        
        self.hud_panel.update_lines(lines)

    def redraw(self):
        """
        Redraw the canvas. All drawing operations are performed inside a lock
        to ensure thread safety with Xlib.
        """
        with self.draw_lock:
            self.canvas.clear()
            self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
            self.update_info_message()
            self.hud_panel.draw(self.canvas)
            self.message_panel.draw(self.canvas)
            self.canvas.flush()

    def fit_molecule_to_window(self):
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
        """
        Build a list of drawable ZObjects and return the first one whose
        'contains()' method reports a hit.
        """
        z_objects = self.renderer.build_render_list(self.ortep_mol, self.view_params)
        for obj in sorted(z_objects, key=lambda o: o.z_value, reverse=True):
            if hasattr(obj, 'contains') and obj.contains(x, y):
                return obj
        return None

    def handle_key(self, evt):
        keysym = self.display.keycode_to_keysym(evt.detail, evt.state)
        keychar = XK.keysym_to_string(keysym)
        if keychar is None:
            return

        if keychar in ("q", "Escape"):
            self.message_service.log_info("Quitting.")
            self.running = False
            sys.exit(0)
        elif keychar == 'h':
            self.view_params.ry -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'l':
            self.view_params.ry += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'j':
            self.view_params.rx += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'k':
            self.view_params.rx -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'u':
            self.view_params.rz -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'o':
            self.view_params.rz += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'H':
            self.view_params.x_offset -= VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'L':
            self.view_params.x_offset += VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'K':
            self.view_params.y_offset -= VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'J':
            self.view_params.y_offset += VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'n':
            self.view_params.scale *= VIEWER_INTERACTION["key_zoom_factor"]
            self.message_service.log_info(f"Zoomed in (scale: {self.view_params.scale:.1f})")
        elif keychar == 'm':
            self.view_params.scale /= VIEWER_INTERACTION["key_zoom_factor"]
            self.message_service.log_info(f"Zoomed out (scale: {self.view_params.scale:.1f})")
        elif keychar == 's':
            self.dump_svg()
        elif keychar == 'd':
            # Toggle hydrogen display.
            self.show_hydrogens = not self.show_hydrogens
            status = "shown" if self.show_hydrogens else "hidden"
            self.message_service.log_info(f"Hydrogens {status}")
            self.set_frame(self.current_frame)
        elif keychar == 'B':
            # Bond Toggle Logic (Uppercase B)
            if len(self.selected_atoms) == 2:
                atom1, atom2 = self.selected_atoms
                toggled = toggle_bond(atom1, atom2, self.ortep_mol)
                if toggled:
                    for a in self.selected_atoms:
                        a.selected = False
                    self.selected_atoms = []
                    toggled.selected = True
                    self.selected_bonds = [toggled]
                    self.selected_bond_ids = [(min(toggled.atom1.index, toggled.atom2.index),
                                                 max(toggled.atom1.index, toggled.atom2.index))]
                    self.message_service.log_info(f"Created bond between {atom1.symbol}{atom1.index} and {atom2.symbol}{atom2.index}")
                else:
                    self.selected_atoms = []
                    self.selected_bonds = []
                    self.selected_bond_ids = []
                    self.message_service.log_info(f"Removed bond between {atom1.symbol}{atom1.index} and {atom2.symbol}{atom2.index}")
        elif keychar == 'b':
            # Bond Cycling Logic (Lowercase b)
            if self.selected_bonds:
                new_bonds = []
                for bond in self.selected_bonds:
                    atom1, atom2 = bond.atom1, bond.atom2
                    old_type = type(bond).__name__.replace("Bond", "")
                    new_bond = cycle_existing_bond(bond, self.ortep_mol)
                    new_bonds.append(new_bond)
                    new_type = type(new_bond).__name__.replace("Bond", "")
                    self.message_service.log_info(f"Changed bond {atom1.symbol}{atom1.index}-{atom2.symbol}{atom2.index} from {old_type} to {new_type}")
                self.selected_bonds = new_bonds
        elif keychar in ('[', ']', '{', '}', '=', '-'):
            # For trajectory navigation, if there's only one frame, ignore the key.
            if self.total_frames <= 1:
                self.message_service.log_info("Only one frame available; ignoring trajectory key press.")
                return
            if keychar == '[':
                self.set_frame(max(self.current_frame - 1, 0))
            elif keychar == ']':
                self.set_frame(min(self.current_frame + 1, self.total_frames - 1))
            elif keychar == '{':
                self.set_frame(0)
                self.message_service.log_info("First frame")
            elif keychar == '}':
                self.set_frame(self.total_frames - 1)
                self.message_service.log_info("Last frame")
            elif keychar == '=':
                # Jump to the highest energy frame.
                energies = self.trajectory._frame_energies if self.trajectory else []
                if energies and any(e is not None for e in energies):
                    max_frame = max(range(len(energies)),
                                    key=lambda i: energies[i] if energies[i] is not None else float('-inf'))
                    self.set_frame(max_frame)
                    self.message_service.log_info(f"Highest energy frame (E={energies[max_frame]:.4f} a.u.)")
            elif keychar == '-':
                # Jump to the lowest energy frame.
                energies = self.trajectory._frame_energies if self.trajectory else []
                if energies and any(e is not None for e in energies):
                    min_frame = min(range(len(energies)),
                                    key=lambda i: energies[i] if energies[i] is not None else float('inf'))
                    self.set_frame(min_frame)
                    self.message_service.log_info(f"Lowest energy frame (E={energies[min_frame]:.4f} a.u.)")

        self.redraw()

    def handle_button_press(self, evt):
        if evt.detail in (4, 5):
            if evt.detail == 4:
                self.view_params.scale *= VIEWER_INTERACTION["mouse_zoom_factor"]
            elif evt.detail == 5:
                self.view_params.scale /= VIEWER_INTERACTION["mouse_zoom_factor"]
            self.redraw()
            return
        self.last_mouse_x = evt.event_x
        self.last_mouse_y = evt.event_y
        if evt.detail == 1:
            self.click_start_x = evt.event_x
            self.click_start_y = evt.event_y
        shift_pressed = bool(evt.state & X.ShiftMask)
        if evt.detail == 1:
            if shift_pressed:
                current_time = time.time()
                if (current_time - self.last_click_time) < 0.3:
                    self.reset_view()
                    return
                else:
                    self.last_click_time = current_time
                    self.active_button = 'shift-left'
            else:
                self.active_button = 'left'
        elif evt.detail == 2:
            self.active_button = 'middle'
        elif evt.detail == 3:
            self.active_button = 'right'

    def handle_motion(self, evt):
        if self.active_button is None:
            return
        dx = evt.event_x - self.last_mouse_x
        dy = evt.event_y - self.last_mouse_y
        self.last_mouse_x = evt.event_x
        self.last_mouse_y = evt.event_y
        sensitivity = VIEWER_INTERACTION["mouse_rotation_sensitivity"]
        if self.active_button == 'left':
            self.view_params.rx += sensitivity * dy
            self.view_params.ry += sensitivity * dx
        elif self.active_button in ('middle', 'shift-left'):
            self.view_params.rz += sensitivity * dx
        elif self.active_button == 'right':
            self.view_params.x_offset += dx
            self.view_params.y_offset += dy
        self.redraw()

    def handle_button_release(self, evt):
        if self.active_button in ('left', 'shift-left'):
            dx = evt.event_x - self.click_start_x
            dy = evt.event_y - self.click_start_y
            if dx * dx + dy * dy < 100:
                clicked_obj = self.hit_test(evt.event_x, evt.event_y)
                shift_pressed = (self.active_button == 'shift-left')
                if clicked_obj is None:
                    for atom in self.selected_atoms:
                        atom.selected = False
                    self.selected_atoms = []
                    for bond in self.selected_bonds:
                        bond.selected = False
                    self.selected_bonds = []
                    self.selected_atom_indices = []
                    self.selected_bond_ids = []
                    # self.message_service.log_info("Selection cleared")
                else:
                    if hasattr(clicked_obj, 'atom') and clicked_obj.atom is not None:
                        for bond in self.selected_bonds:
                            bond.selected = False
                        self.selected_bonds = []
                        underlying_atom = clicked_obj.atom
                        if shift_pressed:
                            if underlying_atom in self.selected_atoms:
                                underlying_atom.selected = False
                                self.selected_atoms.remove(underlying_atom)
                                if underlying_atom.index in self.selected_atom_indices:
                                    self.selected_atom_indices.remove(underlying_atom.index)
                                self.message_service.log_info(f"Deselected {underlying_atom.symbol}{underlying_atom.index}")
                            else:
                                underlying_atom.selected = True
                                self.selected_atoms.append(underlying_atom)
                                if underlying_atom.index not in self.selected_atom_indices:
                                    self.selected_atom_indices.append(underlying_atom.index)
                                self.message_service.log_info(f"Added {underlying_atom.symbol}{underlying_atom.index} to selection")
                        else:
                            for atom in self.selected_atoms:
                                atom.selected = False
                            self.selected_atoms = [underlying_atom]
                            underlying_atom.selected = True
                            self.selected_atom_indices = [underlying_atom.index]
                            self.message_service.log_info(f"Selected {underlying_atom.symbol}{underlying_atom.index}")
                    elif hasattr(clicked_obj, 'bond') and clicked_obj.bond is not None:
                        for atom in self.selected_atoms:
                            atom.selected = False
                        self.selected_atoms = []
                        underlying_bond = clicked_obj.bond
                        if shift_pressed:
                            if underlying_bond in self.selected_bonds:
                                underlying_bond.selected = False
                                self.selected_bonds.remove(underlying_bond)
                                key = (min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                       max(underlying_bond.atom1.index, underlying_bond.atom2.index))
                                if key in self.selected_bond_ids:
                                    self.selected_bond_ids.remove(key)
                                self.message_service.log_info(f"Deselected bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-{underlying_bond.atom2.symbol}{underlying_bond.atom2.index}")
                            else:
                                underlying_bond.selected = True
                                self.selected_bonds.append(underlying_bond)
                                key = (min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                       max(underlying_bond.atom1.index, underlying_bond.atom2.index))
                                if key not in self.selected_bond_ids:
                                    self.selected_bond_ids.append(key)
                                self.message_service.log_info(f"Added bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-{underlying_bond.atom2.symbol}{underlying_bond.atom2.index} to selection")
                        else:
                            for bond in self.selected_bonds:
                                bond.selected = False
                            self.selected_bonds = [underlying_bond]
                            underlying_bond.selected = True
                            self.selected_bond_ids = [(min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                                       max(underlying_bond.atom1.index, underlying_bond.atom2.index))]
                            self.message_service.log_info(f"Selected bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-{underlying_bond.atom2.symbol}{underlying_bond.atom2.index}")
                    else:
                        for atom in self.selected_atoms:
                            atom.selected = False
                        self.selected_atoms = []
                        for bond in self.selected_bonds:
                            bond.selected = False
                        self.selected_bonds = []
                        self.selected_atom_indices = []
                        self.selected_bond_ids = []
                        # self.message_service.log_info("Selection cleared")
                self.redraw()
        self.active_button = None

    def dump_svg(self):
        svg_canvas = SVGCanvas(width=self.canvas.width, height=self.canvas.height)
        self.view_params.as_float = True
        self.renderer.draw_molecule(svg_canvas, self.ortep_mol, self.view_params)
        self.view_params.as_float = False
        filename = "quickORTEP_export.svg"
        svg_canvas.flush(filename)
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
