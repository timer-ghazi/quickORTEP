#!/usr/bin/env python3
"""
event_handler.py

This module provides the _EventHandler class, which encapsulates the logic
for handling keyboard and mouse events. It maps raw X11 events to high-level
commands and delegates the corresponding actions to the MoleculeViewer.
"""

import sys
import time
from Xlib import XK, X
from bond_manager import cycle_existing_bond, toggle_bond
from config import VIEWER_INTERACTION
from x11view.svg_canvas import SVGCanvas

class _EventHandler:
    def __init__(self, viewer):
        """
        Initialize the event handler with a reference to the viewer.
        
        Parameters:
            viewer: An instance of MoleculeViewer that provides access to
                    view parameters, renderer, message service, etc.
        """
        self.viewer = viewer
        self.last_click_time = 0.0
        self.click_start_x = None
        self.click_start_y = None

    def handle_key(self, evt):
        keysym = self.viewer.display.keycode_to_keysym(evt.detail, evt.state)
        keychar = XK.keysym_to_string(keysym)
        if keychar is None:
            return

        if keychar in ("q", "Escape"):
            self.viewer.message_service.log_info("Quitting.")
            self.viewer.running = False
            sys.exit(0)
        elif keychar == 'h':
            self.viewer.view_params.ry -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'l':
            self.viewer.view_params.ry += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'j':
            self.viewer.view_params.rx += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'k':
            self.viewer.view_params.rx -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'u':
            self.viewer.view_params.rz -= VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'o':
            self.viewer.view_params.rz += VIEWER_INTERACTION["rotation_increment_deg"]
        elif keychar == 'H':
            self.viewer.view_params.x_offset -= VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'L':
            self.viewer.view_params.x_offset += VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'K':
            self.viewer.view_params.y_offset -= VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'J':
            self.viewer.view_params.y_offset += VIEWER_INTERACTION["pan_increment"]
        elif keychar == 'n':
            self.viewer.view_params.scale *= VIEWER_INTERACTION["key_zoom_factor"]
            self.viewer.message_service.log_info(
                f"Zoomed in (scale: {self.viewer.view_params.scale:.1f})"
            )
        elif keychar == 'm':
            self.viewer.view_params.scale /= VIEWER_INTERACTION["key_zoom_factor"]
            self.viewer.message_service.log_info(
                f"Zoomed out (scale: {self.viewer.view_params.scale:.1f})"
            )
        elif keychar == 's':
            # Delegate to export functionality.
            self.viewer.dump_svg()
        elif keychar == 'd':
            # Toggle hydrogen display.
            self.viewer.show_hydrogens = not self.viewer.show_hydrogens
            status = "shown" if self.viewer.show_hydrogens else "hidden"
            self.viewer.message_service.log_info(f"Hydrogens {status}")
            self.viewer.set_frame(self.viewer.current_frame)
        elif keychar == 'B':
            # Bond Toggle Logic (Uppercase B)
            if len(self.viewer.selected_atoms) == 2:
                atom1, atom2 = self.viewer.selected_atoms
                toggled = toggle_bond(atom1, atom2, self.viewer.ortep_mol)
                if toggled:
                    for a in self.viewer.selected_atoms:
                        a.selected = False
                    self.viewer.selected_atoms = []
                    toggled.selected = True
                    self.viewer.selected_bonds = [toggled]
                    self.viewer.selected_bond_ids = [
                        (min(toggled.atom1.index, toggled.atom2.index),
                         max(toggled.atom1.index, toggled.atom2.index))
                    ]
                    self.viewer.message_service.log_info(
                        f"Created bond between {atom1.symbol}{atom1.index} and {atom2.symbol}{atom2.index}"
                    )
                else:
                    self.viewer.selected_atoms = []
                    self.viewer.selected_bonds = []
                    self.viewer.selected_bond_ids = []
                    self.viewer.message_service.log_info(
                        f"Removed bond between {atom1.symbol}{atom1.index} and {atom2.symbol}{atom2.index}"
                    )
        elif keychar == 'b':
            # Bond Cycling Logic (Lowercase b)
            if self.viewer.selected_bonds:
                new_bonds = []
                for bond in self.viewer.selected_bonds:
                    atom1, atom2 = bond.atom1, bond.atom2
                    old_type = type(bond).__name__.replace("Bond", "")
                    new_bond = cycle_existing_bond(bond, self.viewer.ortep_mol)
                    new_bonds.append(new_bond)
                    new_type = type(new_bond).__name__.replace("Bond", "")
                    self.viewer.message_service.log_info(
                        f"Changed bond {atom1.symbol}{atom1.index}-{atom2.symbol}{atom2.index} from {old_type} to {new_type}"
                    )
                self.viewer.selected_bonds = new_bonds
        elif keychar in ('[', ']', '{', '}', '=', '-'):
            if self.viewer.total_frames <= 1:
                self.viewer.message_service.log_info("Only one frame available; ignoring trajectory key press.")
                return
            if keychar == '[':
                self.viewer.set_frame(max(self.viewer.current_frame - 1, 0))
            elif keychar == ']':
                self.viewer.set_frame(min(self.viewer.current_frame + 1, self.viewer.total_frames - 1))
            elif keychar == '{':
                self.viewer.set_frame(0)
                self.viewer.message_service.log_info("First frame")
            elif keychar == '}':
                self.viewer.set_frame(self.viewer.total_frames - 1)
                self.viewer.message_service.log_info("Last frame")
            elif keychar == '=':
                energies = self.viewer.trajectory._frame_energies if self.viewer.trajectory else []
                if energies and any(e is not None for e in energies):
                    max_frame = max(
                        range(len(energies)),
                        key=lambda i: energies[i] if energies[i] is not None else float('-inf')
                    )
                    self.viewer.set_frame(max_frame)
                    self.viewer.message_service.log_info(
                        f"Highest energy frame (E={energies[max_frame]:.4f} a.u.)"
                    )
            elif keychar == '-':
                energies = self.viewer.trajectory._frame_energies if self.viewer.trajectory else []
                if energies and any(e is not None for e in energies):
                    min_frame = min(
                        range(len(energies)),
                        key=lambda i: energies[i] if energies[i] is not None else float('inf')
                    )
                    self.viewer.set_frame(min_frame)
                    self.viewer.message_service.log_info(
                        f"Lowest energy frame (E={energies[min_frame]:.4f} a.u.)"
                    )

        self.viewer.redraw()

    def handle_button_press(self, evt):
        # Zooming via mouse wheel.
        if evt.detail in (4, 5):
            if evt.detail == 4:
                self.viewer.view_params.scale *= VIEWER_INTERACTION["mouse_zoom_factor"]
            elif evt.detail == 5:
                self.viewer.view_params.scale /= VIEWER_INTERACTION["mouse_zoom_factor"]
            self.viewer.redraw()
            return

        self.viewer.last_mouse_x = evt.event_x
        self.viewer.last_mouse_y = evt.event_y
        if evt.detail == 1:
            self.click_start_x = evt.event_x
            self.click_start_y = evt.event_y

        shift_pressed = bool(evt.state & X.ShiftMask)
        if evt.detail == 1:
            if shift_pressed:
                current_time = time.time()
                if (current_time - self.last_click_time) < 0.3:
                    self.viewer.reset_view()
                    return
                else:
                    self.last_click_time = current_time
                    self.viewer.active_button = 'shift-left'
            else:
                self.viewer.active_button = 'left'
        elif evt.detail == 2:
            self.viewer.active_button = 'middle'
        elif evt.detail == 3:
            self.viewer.active_button = 'right'

    def handle_motion(self, evt):
        if self.viewer.active_button is None:
            return
        dx = evt.event_x - self.viewer.last_mouse_x
        dy = evt.event_y - self.viewer.last_mouse_y
        self.viewer.last_mouse_x = evt.event_x
        self.viewer.last_mouse_y = evt.event_y
        sensitivity = VIEWER_INTERACTION["mouse_rotation_sensitivity"]
        if self.viewer.active_button == 'left':
            self.viewer.view_params.rx += sensitivity * dy
            self.viewer.view_params.ry += sensitivity * dx
        elif self.viewer.active_button in ('middle', 'shift-left'):
            self.viewer.view_params.rz += sensitivity * dx
        elif self.viewer.active_button == 'right':
            self.viewer.view_params.x_offset += dx
            self.viewer.view_params.y_offset += dy
        self.viewer.redraw()

    def handle_button_release(self, evt):
        if self.viewer.active_button in ('left', 'shift-left'):
            dx = evt.event_x - self.click_start_x
            dy = evt.event_y - self.click_start_y
            if dx * dx + dy * dy < 100:
                clicked_obj = self.viewer.hit_test(evt.event_x, evt.event_y)
                shift_pressed = (self.viewer.active_button == 'shift-left')
                if clicked_obj is None:
                    # Clear selection if nothing was clicked.
                    self.viewer.selected_atoms.clear()
                    self.viewer.selected_bonds.clear()
                    self.viewer.selected_atom_indices.clear()
                    self.viewer.selected_bond_ids.clear()
                else:
                    if hasattr(clicked_obj, 'atom') and clicked_obj.atom is not None:
                        # Clear bond selections if an atom is clicked.
                        self.viewer.selected_bonds.clear()
                        underlying_atom = clicked_obj.atom
                        if shift_pressed:
                            if underlying_atom in self.viewer.selected_atoms:
                                underlying_atom.selected = False
                                self.viewer.selected_atoms.remove(underlying_atom)
                                if underlying_atom.index in self.viewer.selected_atom_indices:
                                    self.viewer.selected_atom_indices.remove(underlying_atom.index)
                                self.viewer.message_service.log_info(
                                    f"Deselected {underlying_atom.symbol}{underlying_atom.index}"
                                )
                            else:
                                underlying_atom.selected = True
                                self.viewer.selected_atoms.append(underlying_atom)
                                if underlying_atom.index not in self.viewer.selected_atom_indices:
                                    self.viewer.selected_atom_indices.append(underlying_atom.index)
                                self.viewer.message_service.log_info(
                                    f"Added {underlying_atom.symbol}{underlying_atom.index} to selection"
                                )
                        else:
                            for atom in self.viewer.selected_atoms:
                                atom.selected = False
                            self.viewer.selected_atoms = [underlying_atom]
                            underlying_atom.selected = True
                            self.viewer.selected_atom_indices = [underlying_atom.index]
                            self.viewer.message_service.log_info(
                                f"Selected {underlying_atom.symbol}{underlying_atom.index}"
                            )
                    elif hasattr(clicked_obj, 'bond') and clicked_obj.bond is not None:
                        # Clear atom selections if a bond is clicked.
                        self.viewer.selected_atoms.clear()
                        underlying_bond = clicked_obj.bond
                        if shift_pressed:
                            if underlying_bond in self.viewer.selected_bonds:
                                underlying_bond.selected = False
                                self.viewer.selected_bonds.remove(underlying_bond)
                                key = (min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                       max(underlying_bond.atom1.index, underlying_bond.atom2.index))
                                if key in self.viewer.selected_bond_ids:
                                    self.viewer.selected_bond_ids.remove(key)
                                self.viewer.message_service.log_info(
                                    f"Deselected bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-"
                                    f"{underlying_bond.atom2.symbol}{underlying_bond.atom2.index}"
                                )
                            else:
                                underlying_bond.selected = True
                                self.viewer.selected_bonds.append(underlying_bond)
                                key = (min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                       max(underlying_bond.atom1.index, underlying_bond.atom2.index))
                                if key not in self.viewer.selected_bond_ids:
                                    self.viewer.selected_bond_ids.append(key)
                                self.viewer.message_service.log_info(
                                    f"Added bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-"
                                    f"{underlying_bond.atom2.symbol}{underlying_bond.atom2.index} to selection"
                                )
                        else:
                            self.viewer.selected_bonds = [underlying_bond]
                            underlying_bond.selected = True
                            self.viewer.selected_bond_ids = [
                                (min(underlying_bond.atom1.index, underlying_bond.atom2.index),
                                 max(underlying_bond.atom1.index, underlying_bond.atom2.index))
                            ]
                            self.viewer.message_service.log_info(
                                f"Selected bond {underlying_bond.atom1.symbol}{underlying_bond.atom1.index}-"
                                f"{underlying_bond.atom2.symbol}{underlying_bond.atom2.index}"
                            )
                    else:
                        # Clear selection if no valid object was clicked.
                        self.viewer.selected_atoms.clear()
                        self.viewer.selected_bonds.clear()
                        self.viewer.selected_atom_indices.clear()
                        self.viewer.selected_bond_ids.clear()
                self.viewer.redraw()
        self.viewer.active_button = None
