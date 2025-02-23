#!/usr/bin/env python3
"""
ortep_viewer.py

The viewer class that manages an ORTEP_Molecule, its renderer, and handles user events.
Now updated to support persistent selection for atoms and bonds by storing selection state
in the underlying objects.
"""

import sys
import time
from Xlib import XK, X
from x11view.window import X11Window, X11CanvasBasic, X11CanvasSS
from x11view.svg_canvas import SVGCanvas
from ortep_renderer import ORTEP_MoleculeRenderer
from geometry_utils import rotate_point

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
    The controller/view class that manages an ORTEP_Molecule, its renderer,
    and handles user events (key presses, mouse events, etc.).

    Now supports persistent multi-selection for atoms and bonds.
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

        # Mouse state variables.
        self.active_button = None         # 'left', 'middle', 'right', or 'shift-left'
        self.last_mouse_x = None
        self.last_mouse_y = None
        self.last_click_time = 0.0

        # Use lists to store selected underlying objects.
        self.selected_atoms = []   # List of ORTEP_Atom objects.
        self.selected_bonds = []   # List of underlying bond objects.
        self.info_message = "No object selected yet."
        self.click_start_x = None
        self.click_start_y = None

    def redraw(self):
        self.canvas.clear()
        self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
        # Draw the info text message at the bottom left.
        self.canvas.draw_text(10, self.canvas.height - 20,
                              self.info_message, color=(0, 0, 0), font_size=14)
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
            print("Quitting.")
            self.running = False
            sys.exit(0)
        elif keychar == 'h':
            self.view_params.ry -= 5
        elif keychar == 'l':
            self.view_params.ry += 5
        elif keychar == 'j':
            self.view_params.rx += 5
        elif keychar == 'k':
            self.view_params.rx -= 5
        elif keychar == 'u':
            self.view_params.rz -= 5
        elif keychar == 'o':
            self.view_params.rz += 5
        elif keychar == 'H':
            self.view_params.x_offset -= 10
        elif keychar == 'L':
            self.view_params.x_offset += 10
        elif keychar == 'K':
            self.view_params.y_offset -= 10
        elif keychar == 'J':
            self.view_params.y_offset += 10
        elif keychar == 'n':
            self.view_params.scale *= 1.05
        elif keychar == 'm':
            self.view_params.scale /= 1.05
        elif keychar == 's':
            self.dump_svg()
        else:
            print(f"Ignored key: {keychar}")
        self.redraw()

    def handle_button_press(self, evt):
        if evt.detail in (4, 5):
            if evt.detail == 4:
                self.view_params.scale *= 1.1
            elif evt.detail == 5:
                self.view_params.scale /= 1.1
            self.redraw()
            return
        self.last_mouse_x = evt.event_x
        self.last_mouse_y = evt.event_y
        if evt.detail == 1:
            self.click_start_x = evt.event_x
            self.click_start_y = evt.event_y
        # Check for Shift modifier.
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
        k_x = 0.5
        k_y = 0.5
        k_z = 0.5
        if self.active_button == 'left':
            self.view_params.rx += k_x * dy
            self.view_params.ry += k_y * dx
        elif self.active_button in ('middle', 'shift-left'):
            self.view_params.rz += k_z * dx
        elif self.active_button == 'right':
            self.view_params.x_offset += dx
            self.view_params.y_offset += dy
        self.redraw()

    def handle_button_release(self, evt):
        if self.active_button in ('left', 'shift-left'):
            dx = evt.event_x - self.click_start_x
            dy = evt.event_y - self.click_start_y
            if dx*dx + dy*dy < 100:  # treat as a click if movement is minimal
                clicked_obj = self.hit_test(evt.event_x, evt.event_y)
                # Determine from active_button whether Shift was held.
                shift_pressed = (self.active_button == 'shift-left')
                if clicked_obj is None:
                    # Clear all selections.
                    for atom in self.selected_atoms:
                        atom.selected = False
                    self.selected_atoms = []
                    for bond in self.selected_bonds:
                        bond.selected = False
                    self.selected_bonds = []
                    self.info_message = "No object selected"
                else:
                    # Process atom selection.
                    if hasattr(clicked_obj, 'atom') and clicked_obj.atom is not None:
                        # Clear bond selections.
                        for bond in self.selected_bonds:
                            bond.selected = False
                        self.selected_bonds = []
                        underlying_atom = clicked_obj.atom
                        if shift_pressed:
                            if underlying_atom in self.selected_atoms:
                                underlying_atom.selected = False
                                self.selected_atoms.remove(underlying_atom)
                            else:
                                underlying_atom.selected = True
                                self.selected_atoms.append(underlying_atom)
                        else:
                            for atom in self.selected_atoms:
                                atom.selected = False
                            self.selected_atoms = [underlying_atom]
                            underlying_atom.selected = True
                        if len(self.selected_atoms) == 1:
                            a = self.selected_atoms[0]
                            self.info_message = f"Atom: {a.symbol}{a.index}"
                        else:
                            sel_info = ", ".join(f"{a.symbol}{a.index}" for a in self.selected_atoms)
                            self.info_message = f"Selected atoms: {sel_info}"
                    # Process bond selection.
                    elif hasattr(clicked_obj, 'bond') and clicked_obj.bond is not None:
                        # Clear atom selections.
                        for atom in self.selected_atoms:
                            atom.selected = False
                        self.selected_atoms = []
                        underlying_bond = clicked_obj.bond
                        if shift_pressed:
                            if underlying_bond in self.selected_bonds:
                                underlying_bond.selected = False
                                self.selected_bonds.remove(underlying_bond)
                            else:
                                underlying_bond.selected = True
                                self.selected_bonds.append(underlying_bond)
                        else:
                            for bond in self.selected_bonds:
                                bond.selected = False
                            self.selected_bonds = [underlying_bond]
                            underlying_bond.selected = True
                        if len(self.selected_bonds) == 1:
                            b = underlying_bond
                            self.info_message = (f"Bond: {b.atom1.symbol}{b.atom1.index} - "
                                                 f"{b.atom2.symbol}{b.atom2.index} {b.length:.4f} Ang")
                        else:
                            sel_info = ", ".join(
                                f"{b.atom1.symbol}{b.atom1.index}-{b.atom2.symbol}{b.atom2.index}"
                                for b in self.selected_bonds)
                            self.info_message = f"Selected bonds: {sel_info}"
                    else:
                        # Unrecognized object; clear selections.
                        for atom in self.selected_atoms:
                            atom.selected = False
                        self.selected_atoms = []
                        for bond in self.selected_bonds:
                            bond.selected = False
                        self.selected_bonds = []
                        self.info_message = "No object selected"
                print(self.info_message)
                self.redraw()
        self.active_button = None

    def dump_svg(self):
        svg_canvas = SVGCanvas(width=self.canvas.width, height=self.canvas.height)
        self.view_params.as_float = True
        self.renderer.draw_molecule(svg_canvas, self.ortep_mol, self.view_params)
        self.view_params.as_float = False
        filename = "quickORTEP_export.svg"
        svg_canvas.flush(filename)
        print(f"SVG export saved to {filename}")

    def reset_view(self):
        self.view_params.rx = 0.0
        self.view_params.ry = 0.0
        self.view_params.rz = 0.0
        self.view_params.scale = 100.0
        self.view_params.x_offset = self.canvas.width // 2
        self.view_params.y_offset = self.canvas.height // 2
        self.redraw()
