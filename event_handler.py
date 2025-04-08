#!/usr/bin/env python3
"""
event_handler.py

This module provides the _EventHandler class, which encapsulates the logic
for handling keyboard and mouse events. It maps raw X11 events to high-level
commands and delegates the corresponding actions to the MoleculeViewer.
"""

import sys
import time
from functools import lru_cache
import numpy as np
from Xlib import XK, X
from Xlib.XK import XK_Left, XK_Right, XK_Up, XK_Down

from bond_manager import cycle_existing_bond, toggle_bond
from config import VIEWER_INTERACTION, ENERGY_UNITS, DEFAULT_ENERGY_UNIT
from x11view.svg_canvas import SVGCanvas


class InputState:
    """Enumeration of possible input states for the event handler."""
    NORMAL = 0      # Default state
    SELECTING = 1   # During selection operations
    ROTATING = 2    # During rotation operations
    PANNING = 3     # During panning operations
    ZOOMING = 4     # During zoom operations


class _EventHandler:
    """
    Handles user input events for the molecule viewer.
    
    This class is responsible for:
    1. Capturing and processing keyboard and mouse events
    2. Mapping raw input events to high-level commands
    3. Delegating actions to the appropriate components of the viewer
    4. Managing selection state for atoms and bonds
    
    The event handler maintains state about:
    - Active mouse buttons
    - Click locations and timing
    - Selection context
    """
    
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
        self.last_motion_time = 0.0  # For throttling motion events
        
        # Initialize state machine
        self.input_state = InputState.NORMAL
        self.state_data = {}  # Store state-specific data
        
        # Set up key command mappings
        self._initialize_key_commands()

    def _initialize_key_commands(self):
        """Initialize the key command mappings."""
        self.key_commands = {
            # Application control
            'q': self._quit_application,
            'Escape': self._quit_application,
            
            # View rotation
            'h': lambda: self._rotate_view('y', -1),
            'l': lambda: self._rotate_view('y', 1),
            'j': lambda: self._rotate_view('x', 1),
            'k': lambda: self._rotate_view('x', -1),
            'u': lambda: self._rotate_view('z', -1),
            'o': lambda: self._rotate_view('z', 1),
            
            # View panning
            'H': lambda: self._pan_view('x', -1),
            'L': lambda: self._pan_view('x', 1),
            'K': lambda: self._pan_view('y', -1),
            'J': lambda: self._pan_view('y', 1),
            
            # Zoom
            'n': self._zoom_in,
            'm': self._zoom_out,
            
            # Export/Save
            's': self._dump_svg,
            'x': self._dump_xyz,
            'w': self._dump_graph_data,
            
            # Display options
            'd': self._toggle_hydrogens,
            'g': self._toggle_grid,
            'a': self._toggle_axes,  # Add axes toggle command
            't': self._standard_orientation, # Standard orientation (current frame)
            'T': self._standard_orientation_all, # Standard orientation (all frames)
            '3': self._toggle_3d_effects,  # Toggle 3D effects (highlights and shadows)
            
            # Bond operations
            'p': self._toggle_graph_mode,
            'P': self._toggle_bond_propagation,  # Shift+p
            'b': self._cycle_selected_bonds,
            'B': self._toggle_selected_atoms_bond,
            'c': self._clear_bond_edits,
            
            # Frame navigation
            '[': lambda: self._change_frame(-1),
            ']': lambda: self._change_frame(1),
            '{': self._goto_first_frame,
            '}': self._goto_last_frame,
            '=': self._goto_highest_energy_frame,
            '-': self._goto_lowest_energy_frame,
            
            # Normal mode visualization
            'v': self._toggle_normal_modes,
            ',': lambda: self._cycle_normal_mode(-1),  # Previous mode
            '.': lambda: self._cycle_normal_mode(1),   # Next mode
            '<': self._decrease_normal_mode_scale,     # Decrease scale
            '>': self._increase_normal_mode_scale,     # Increase scale
            
            # Reset/fit operations
            'f': self._fit_molecule_to_window,
            'r': self._reset_view,

            'R': self._reload_file,

            '?': self._toggle_help,
        }
        
        # Override with configurable keys if present
        axes_toggle_key = VIEWER_INTERACTION.get("axes_toggle_key", "a")
        if axes_toggle_key != "a":
            # Remove the default 'a' mapping
            if 'a' in self.key_commands:
                del self.key_commands['a']
            # Add the configured key for toggling axes
            self.key_commands[axes_toggle_key] = self._toggle_axes

    def handle_key(self, evt):
        """
        Handle keyboard events and dispatch to appropriate handler methods.
        REVISED: Handles arrow keys by checking base keysym first, then shift state.
    
        Parameters:
            evt: The X11 key event to process.
        """
        try:
            # Get the base keysym WITHOUT considering modifiers first
            # This tells us which physical key was pressed (like 'Left Arrow')
            base_keysym = self.viewer.display.keycode_to_keysym(evt.detail, 0)
    
            # Now check the actual state for modifiers like Shift
            shift_pressed = bool(evt.state & X.ShiftMask)
    
            if self.viewer.show_help:
                self._toggle_help()
                return
    
            action_taken = False
    
            # --- Arrow Key Logic (Check base_keysym first) ---
            if base_keysym == XK_Left:
                if shift_pressed:
                    self._pan_view('x', -1) # Pan left
                else:
                    self._rotate_view('y', -1) # Rotate left
                action_taken = True
            elif base_keysym == XK_Right:
                if shift_pressed:
                    self._pan_view('x', 1)  # Pan right
                else:
                    self._rotate_view('y', 1)  # Rotate right
                action_taken = True
            elif base_keysym == XK_Up:
                if shift_pressed:
                    self._pan_view('y', -1) # Pan up
                else:
                    self._rotate_view('x', -1) # Rotate up
                action_taken = True
            elif base_keysym == XK_Down:
                if shift_pressed:
                    self._pan_view('y', 1)  # Pan down
                else:
                    self._rotate_view('x', 1)  # Rotate down
                action_taken = True
    
            # --- Existing Key Command Logic ---
            # Only check other keys if no arrow key action was taken
            if not action_taken:
                # For other keys, we might still need the keysym derived *with* state
                # to handle shifted symbols or letters correctly.
                keysym_with_state = self.viewer.display.keycode_to_keysym(evt.detail, evt.state)
                keychar = XK.keysym_to_string(keysym_with_state) # Use string from stateful keysym
    
                if keychar is not None:
                     # The key_commands dictionary expects the final character ('h', 'H', '?', etc.)
                     key_to_check = keychar
                     if key_to_check in self.key_commands:
                        self.key_commands[key_to_check]()
                        action_taken = True
    
    
            # Redraw if any action was taken
            if action_taken:
                self.viewer.redraw()
    
        except Exception as e:
            print(f"ERROR handling key event: {str(e)}")
            # import traceback
            # traceback.print_exc()
            # self.viewer.message_service.log_error(f"Error handling key event: {str(e)}")
            # self.viewer.message_service.log_debug(traceback.format_exc())        


    def handle_button_press(self, evt):
        """
        Handle mouse button press events.
        
        Parameters:
            evt: The X11 button press event.
        """
        try:
            # Zooming via mouse wheel
            if evt.detail in (4, 5):
                if evt.detail == 4:
                    self.viewer.view_params.scale *= VIEWER_INTERACTION["mouse_zoom_factor"]
                elif evt.detail == 5:
                    self.viewer.view_params.scale /= VIEWER_INTERACTION["mouse_zoom_factor"]
                # Update the current zoom level in the VIEWER_INTERACTION dictionary
                VIEWER_INTERACTION["current_zoom"] = self.viewer.view_params.scale
                self.viewer.redraw()
                return

            self.viewer.last_mouse_x = evt.event_x
            self.viewer.last_mouse_y = evt.event_y
            
            # Update state based on button and modifiers
            if evt.detail == 1:  # Left button
                self.click_start_x = evt.event_x
                self.click_start_y = evt.event_y
                
                shift_pressed = bool(evt.state & X.ShiftMask)
                if shift_pressed:
                    current_time = time.time()
                    if (current_time - self.last_click_time) < 0.3:
                        self._reset_view()
                        return
                    else:
                        self.last_click_time = current_time
                        self.input_state = InputState.SELECTING
                        self.state_data['multi_select'] = True
                        self.viewer.active_button = 'shift-left'
                else:
                    self.input_state = InputState.ROTATING
                    self.viewer.active_button = 'left'
            elif evt.detail == 2:  # Middle button
                self.input_state = InputState.ROTATING
                self.state_data['rotate_z'] = True
                self.viewer.active_button = 'middle'
            elif evt.detail == 3:  # Right button
                self.input_state = InputState.PANNING
                self.viewer.active_button = 'right'
        except Exception as e:
            self.viewer.message_service.log_error(f"Error in button press: {str(e)}")
            self.input_state = InputState.NORMAL
            self.viewer.active_button = None

    def handle_motion(self, evt):
        """
        Handle mouse motion events with throttling for better performance.
        
        Parameters:
            evt: The X11 mouse motion event.
        """
        try:
            if self.viewer.active_button is None:
                return
                
            # Throttle updates for smoother performance (~60 FPS)
            current_time = time.time()
            if current_time - self.last_motion_time < 0.016:
                return
            self.last_motion_time = current_time
            
            # Calculate movement deltas
            dx = evt.event_x - self.viewer.last_mouse_x
            dy = evt.event_y - self.viewer.last_mouse_y
            self.viewer.last_mouse_x = evt.event_x
            self.viewer.last_mouse_y = evt.event_y
            
            # Apply transformations based on state and active button
            sensitivity = VIEWER_INTERACTION["mouse_rotation_sensitivity"]
            
            if self.input_state == InputState.ROTATING:
                if self.state_data.get('rotate_z', False) or self.viewer.active_button == 'middle' or self.viewer.active_button == 'shift-left':
                    self.viewer.view_params.rz += sensitivity * dx
                else:
                    self.viewer.view_params.rx += sensitivity * dy
                    self.viewer.view_params.ry += sensitivity * dx
            elif self.input_state == InputState.PANNING:
                self.viewer.view_params.x_offset += dx
                self.viewer.view_params.y_offset += dy
                
            self.viewer.redraw()
        except Exception as e:
            self.viewer.message_service.log_error(f"Error in motion handler: {str(e)}")

    def handle_button_release(self, evt):
        """
        Handle mouse button release events, primarily for selection.
        
        This method implements the selection logic:
        - A click (small movement between press and release) selects an object
        - With Shift key, selection is added/removed from the current selection
        - Without Shift key, the current selection is replaced
        - Clicking on empty space clears the selection
        
        Parameters:
            evt: The X11 button release event.
        """
        try:
            if self.viewer.active_button in ('left', 'shift-left'):
                dx = evt.event_x - self.click_start_x
                dy = evt.event_y - self.click_start_y
                if dx * dx + dy * dy < 100:  # Small movement = click
                    clicked_obj = self.viewer.hit_test(evt.event_x, evt.event_y)
                    shift_pressed = (self.viewer.active_button == 'shift-left')
                    
                    if clicked_obj is None:
                        # Clear selection if nothing was clicked
                        self._clear_all_selections()
                    else:
                        if hasattr(clicked_obj, 'atom') and clicked_obj.atom is not None:
                            self._handle_atom_selection(clicked_obj.atom, shift_pressed)
                        elif hasattr(clicked_obj, 'bond') and clicked_obj.bond is not None:
                            self._handle_bond_selection(clicked_obj.bond, shift_pressed)
                        elif hasattr(clicked_obj, 'vector') and clicked_obj.vector is not None:
                            # Add vector selection handling
                            self._handle_vector_selection(clicked_obj.vector, shift_pressed)
                        else:
                            # Clear selection if no valid object was clicked
                            self._clear_all_selections()
                    
                    self.viewer.redraw()
            
            # Reset state
            self.input_state = InputState.NORMAL
            self.state_data.clear()
            self.viewer.active_button = None
        except Exception as e:
            self.viewer.message_service.log_error(f"Error in button release: {str(e)}")
            self.input_state = InputState.NORMAL
            self.viewer.active_button = None

    # Command implementation methods

    def _standard_orientation(self):
        """Convert the current frame to standard orientation."""
        self.viewer.standard_orientation()
        self.viewer.redraw()
        
    def _standard_orientation_all(self):
        """Convert all frames in the trajectory to standard orientation."""
        self.viewer.standard_orientation_all()
        self.viewer.redraw()
    
    def _reload_file(self):
        """Reload the current file from disk."""
        self.viewer.reload_file()
        self.viewer.redraw()
    
    def _quit_application(self):
        """Quit the application."""
        self.viewer.message_service.log_info("Quitting.")
        self.viewer.running = False
        sys.exit(0)
    
    def _rotate_view(self, axis, direction):
        """
        Rotate the view around the specified axis.
        
        Parameters:
            axis: The axis to rotate around ('x', 'y', or 'z').
            direction: The direction of rotation (1 or -1).
        """
        increment = VIEWER_INTERACTION["rotation_increment_deg"] * direction
        if axis == 'x':
            self.viewer.view_params.rx += increment
        elif axis == 'y':
            self.viewer.view_params.ry += increment
        elif axis == 'z':
            self.viewer.view_params.rz += increment
    
    def _pan_view(self, axis, direction):
        """
        Pan the view along the specified axis.
        
        Parameters:
            axis: The axis to pan along ('x' or 'y').
            direction: The direction of panning (1 or -1).
        """
        increment = VIEWER_INTERACTION["pan_increment"] * direction
        if axis == 'x':
            self.viewer.view_params.x_offset += increment
        elif axis == 'y':
            self.viewer.view_params.y_offset += increment
    
    def _zoom_in(self):
        """Zoom in the view."""
        self.viewer.view_params.scale *= VIEWER_INTERACTION["key_zoom_factor"]
        # Update the current zoom level in the VIEWER_INTERACTION dictionary
        VIEWER_INTERACTION["current_zoom"] = self.viewer.view_params.scale
        self.viewer.message_service.log_info(
            f"Zoomed in (scale: {self.viewer.view_params.scale:.1f})"
        )
    
    def _zoom_out(self):
        """Zoom out the view."""
        self.viewer.view_params.scale /= VIEWER_INTERACTION["key_zoom_factor"]
        # Update the current zoom level in the VIEWER_INTERACTION dictionary
        VIEWER_INTERACTION["current_zoom"] = self.viewer.view_params.scale
        self.viewer.message_service.log_info(
            f"Zoomed out (scale: {self.viewer.view_params.scale:.1f})"
        )
    
    def _dump_svg(self):
        """Export the current view as SVG."""
        self.viewer.dump_svg()

    def _dump_xyz(self):
        """Export the current frame as XYZ."""
        self.viewer.dump_xyz()

    def _dump_graph_data(self):
        """Dump the current graph data to a file."""
        self.viewer.dump_graph_data()
    
    def _toggle_hydrogens(self):
        """Toggle the display of hydrogen atoms."""
        self.viewer.show_hydrogens = not self.viewer.show_hydrogens
        status = "shown" if self.viewer.show_hydrogens else "hidden"
        self.viewer.message_service.log_info(f"Hydrogens {status}")
        self.viewer.set_frame(self.viewer.current_frame)
    
    def _toggle_grid(self):
        """Toggle the visibility of the grid."""
        self.viewer.toggle_grid()
    
    def _toggle_axes(self):
        """Toggle the visibility of the coordinate axes."""
        self.viewer.toggle_axes()
        
    def _toggle_3d_effects(self):
        """Toggle 3D effects (atom highlights and shadows)."""
        self.viewer.toggle_3d_effects()
    
    def _toggle_normal_modes(self):
        """Toggle the display of normal mode vectors."""
        self.viewer.toggle_normal_modes()
    
    def _cycle_normal_mode(self, direction):
        """
        Cycle to the next or previous normal mode.
        
        Parameters:
            direction: 1 for next, -1 for previous
        """
        self.viewer.cycle_normal_mode(direction)
    
    def _increase_normal_mode_scale(self):
        """Increase the scale factor for normal mode vectors."""
        self.viewer.adjust_normal_mode_scale(1.25)  # Increase by 25%
    
    def _decrease_normal_mode_scale(self):
        """Decrease the scale factor for normal mode vectors."""
        self.viewer.adjust_normal_mode_scale(0.8)  # Decrease by 20%
    
    def _toggle_graph_mode(self):
        """Toggle the graph mode between energy and bond length."""
        self.viewer.toggle_graph_mode()
    
    def _toggle_bond_propagation(self):
        """Toggle the bond propagation feature."""
        enabled = self.viewer.bond_edit_tracker.toggle()
        status = "enabled" if enabled else "disabled"
        self.viewer.message_service.log_info(f"Bond propagation {status}")
        # Re-apply current frame to reflect the change
        self.viewer.set_frame(self.viewer.current_frame)
    
    def _cycle_selected_bonds(self):
        """Cycle through bond types for selected bonds."""
        if self.viewer.selected_bonds:
            new_bonds = []
            for bond in self.viewer.selected_bonds:
                atom1, atom2 = bond.atom1, bond.atom2
                old_type = type(bond).__name__.replace("Bond", "")
                new_bond = cycle_existing_bond(bond, self.viewer.ortep_mol)
                
                # Record the edit in the bond tracker
                self.viewer.bond_edit_tracker.change_bond_type(
                    atom1.index, atom2.index, type(new_bond))
                
                new_bonds.append(new_bond)
                new_type = type(new_bond).__name__.replace("Bond", "")
                self.viewer.message_service.log_info(
                    f"Changed bond {atom1.symbol}{atom1.index}-{atom2.symbol}{atom2.index} from {old_type} to {new_type}"
                )
            self.viewer.selected_bonds = new_bonds
    
    def _toggle_selected_atoms_bond(self):
        """Toggle a bond between two selected atoms."""
        if len(self.viewer.selected_atoms) == 2:
            atom1, atom2 = self.viewer.selected_atoms
            toggled = toggle_bond(atom1, atom2, self.viewer.ortep_mol)
            
            # Record the edit in the bond tracker
            if toggled:
                self.viewer.bond_edit_tracker.add_bond(
                    atom1.index, atom2.index, type(toggled))
                
                self._clear_atom_selections()
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
                self.viewer.bond_edit_tracker.remove_bond(
                    atom1.index, atom2.index)
                
                self._clear_atom_selections()
                self._clear_bond_selections()
                self.viewer.message_service.log_info(
                    f"Removed bond between {atom1.symbol}{atom1.index} and {atom2.symbol}{atom2.index}"
                )
    
    def _clear_bond_edits(self):
        """Clear all bond edits."""
        self.viewer.clear_bond_edits()
    
    def _change_frame(self, delta):
        """
        Change the current frame by the specified delta.
        
        Parameters:
            delta: The frame change (positive or negative).
        """
        if self.viewer.total_frames <= 1:
            self.viewer.message_service.log_info("Only one frame available; ignoring trajectory key press.")
            return
            
        new_frame = self.viewer.current_frame + delta
        new_frame = max(0, min(new_frame, self.viewer.total_frames - 1))
        self.viewer.set_frame(new_frame)
    
    def _goto_first_frame(self):
        """Go to the first frame of the trajectory."""
        if self.viewer.total_frames <= 1:
            self.viewer.message_service.log_info("Only one frame available; ignoring trajectory key press.")
            return
            
        self.viewer.set_frame(0)
        self.viewer.message_service.log_info("First frame")
    
    def _goto_last_frame(self):
        """Go to the last frame of the trajectory."""
        if self.viewer.total_frames <= 1:
            self.viewer.message_service.log_info("Only one frame available; ignoring trajectory key press.")
            return
            
        self.viewer.set_frame(self.viewer.total_frames - 1)
        self.viewer.message_service.log_info("Last frame")
    
    def _goto_highest_energy_frame(self):
        """Go to the frame with the highest energy."""
        if self.viewer.total_frames <= 1 or self.viewer.trajectory is None:
            self.viewer.message_service.log_info("No trajectory available; ignoring energy frame navigation.")
            return
            
        # Get energies with proper unit conversion and unpack the tuple
        energies, energy_info = self.viewer.trajectory.energy_trajectory(
            convert_if_hartrees=True, 
            convert_to_unit=DEFAULT_ENERGY_UNIT
        )
        
        if len(energies) > 0 and not np.isnan(energies).all():
            # Find index of maximum energy (ignoring NaN values)
            max_frame = np.nanargmax(energies)
            self.viewer.set_frame(max_frame)
            
            # Get unit symbol for display
            unit_symbol = ENERGY_UNITS.get(
                energy_info['converted_unit'], 
                {'symbol': 'a.u.'}
            )['symbol']
            
            # Display with proper units
            self.viewer.message_service.log_info(
                f"Highest energy frame (E={energies[max_frame]:.4f} {unit_symbol})"
            )
    
    def _goto_lowest_energy_frame(self):
        """Go to the frame with the lowest energy."""
        if self.viewer.total_frames <= 1 or self.viewer.trajectory is None:
            self.viewer.message_service.log_info("No trajectory available; ignoring energy frame navigation.")
            return
            
        # Get energies with proper unit conversion and unpack the tuple
        energies, energy_info = self.viewer.trajectory.energy_trajectory(
            convert_if_hartrees=True, 
            convert_to_unit=DEFAULT_ENERGY_UNIT
        )
        
        if len(energies) > 0 and not np.isnan(energies).all():
            # Find index of minimum energy (ignoring NaN values)
            min_frame = np.nanargmin(energies)
            self.viewer.set_frame(min_frame)
            
            # Get unit symbol for display
            unit_symbol = ENERGY_UNITS.get(
                energy_info['converted_unit'], 
                {'symbol': 'a.u.'}
            )['symbol']
            
            # Display with proper units
            self.viewer.message_service.log_info(
                f"Lowest energy frame (E={energies[min_frame]:.4f} {unit_symbol})"
            )
    
    def _fit_molecule_to_window(self):
        """Fit the molecule to the window."""
        self.viewer.fit_molecule_to_window()
    
    def _reset_view(self):
        """Reset the view to default parameters."""
        self.viewer.reset_view()
    
    # Selection handling methods
    
    def _handle_atom_selection(self, atom, shift_pressed):
        """
        Handle atom selection and deselection.
        
        Parameters:
            atom: The atom to select/deselect.
            shift_pressed: Whether Shift key is pressed (for multi-selection).
        """
        # Clear bond selections when an atom is clicked
        self._clear_bond_selections()
        
        if shift_pressed:
            # Toggle selection with Shift key
            if atom in self.viewer.selected_atoms:
                # Deselect the atom
                atom.selected = False
                self.viewer.selected_atoms.remove(atom)
                if atom.index in self.viewer.selected_atom_indices:
                    self.viewer.selected_atom_indices.remove(atom.index)
                self.viewer.message_service.log_info(
                    f"Deselected {atom.symbol}{atom.index}"
                )
            else:
                # Add to selection
                atom.selected = True
                self.viewer.selected_atoms.append(atom)
                if atom.index not in self.viewer.selected_atom_indices:
                    self.viewer.selected_atom_indices.append(atom.index)
                self.viewer.message_service.log_info(
                    f"Added {atom.symbol}{atom.index} to selection"
                )
        else:
            # Replace selection without Shift key
            self._clear_atom_selections()
            self.viewer.selected_atoms = [atom]
            atom.selected = True
            self.viewer.selected_atom_indices = [atom.index]
            self.viewer.message_service.log_info(
                f"Selected {atom.symbol}{atom.index}"
            )
    
    def _handle_bond_selection(self, bond, shift_pressed):
        """
        Handle bond selection and deselection.
        
        Parameters:
            bond: The bond to select/deselect.
            shift_pressed: Whether Shift key is pressed (for multi-selection).
        """
        # Clear atom selections when a bond is clicked
        self._clear_atom_selections()
        
        if shift_pressed:
            # Toggle selection with Shift key
            if bond in self.viewer.selected_bonds:
                # Deselect the bond
                bond.selected = False
                self.viewer.selected_bonds.remove(bond)
                key = (min(bond.atom1.index, bond.atom2.index),
                       max(bond.atom1.index, bond.atom2.index))
                if key in self.viewer.selected_bond_ids:
                    self.viewer.selected_bond_ids.remove(key)
                self.viewer.message_service.log_info(
                    f"Deselected bond {bond.atom1.symbol}{bond.atom1.index}-"
                    f"{bond.atom2.symbol}{bond.atom2.index}"
                )
            else:
                # Add to selection
                bond.selected = True
                self.viewer.selected_bonds.append(bond)
                key = (min(bond.atom1.index, bond.atom2.index),
                       max(bond.atom1.index, bond.atom2.index))
                if key not in self.viewer.selected_bond_ids:
                    self.viewer.selected_bond_ids.append(key)
                self.viewer.message_service.log_info(
                    f"Added bond {bond.atom1.symbol}{bond.atom1.index}-"
                    f"{bond.atom2.symbol}{bond.atom2.index} to selection"
                )
        else:
            # Replace selection without Shift key
            self._clear_bond_selections()
            self.viewer.selected_bonds = [bond]
            bond.selected = True
            self.viewer.selected_bond_ids = [
                (min(bond.atom1.index, bond.atom2.index),
                 max(bond.atom1.index, bond.atom2.index))
            ]
            self.viewer.message_service.log_info(
                f"Selected bond {bond.atom1.symbol}{bond.atom1.index}-"
                f"{bond.atom2.symbol}{bond.atom2.index}"
            )
            
            # Suggest bond length graph if in energy mode
            #if self.viewer._graph_manager.energy_graph is not None and self.viewer._graph_manager.graph_mode == "energy":
            if self.viewer.energy_graph is not None and self.viewer.graph_mode == "energy":
                self.viewer.message_service.log_info(
                    "Press 'p' to show bond length plot"
                )
    
    def _handle_vector_selection(self, vector, shift_pressed):
        """
        Handle vector selection and deselection.
        
        Parameters:
            vector: The vector to select/deselect.
            shift_pressed: Whether Shift key is pressed (for multi-selection).
        """
        # Clear other selections when a vector is clicked
        self._clear_atom_selections()
        self._clear_bond_selections()
        
        # Toggle vector selection state
        vector.selected = not vector.selected
        
        # Inform the user
        status = "Selected" if vector.selected else "Deselected"
        # If it's an axis vector, use a more descriptive message
        if hasattr(vector, '_is_axis_vector'):
            axis_names = ["X", "Y", "Z"]
            # Find which axis this is
            for i, axis_vector in enumerate(self.viewer.ortep_mol.axis_system.get_vectors()):
                if vector is axis_vector and i < len(axis_names):
                    self.viewer.message_service.log_info(
                        f"{status} {axis_names[i]} axis"
                    )
                    break
            else:
                self.viewer.message_service.log_info(f"{status} axis vector")
        else:
            self.viewer.message_service.log_info(f"{status} vector")
    
    def _clear_atom_selections(self):
        """Clear all atom selections."""
        for atom in self.viewer.selected_atoms:
            atom.selected = False
        self.viewer.selected_atoms.clear()
        self.viewer.selected_atom_indices.clear()
    
    def _clear_bond_selections(self):
        """Clear all bond selections."""
        for bond in self.viewer.selected_bonds:
            bond.selected = False
        self.viewer.selected_bonds.clear()
        self.viewer.selected_bond_ids.clear()
    
    def _clear_all_selections(self):
        """Clear all selections (atoms and bonds)."""
        self._clear_atom_selections()
        self._clear_bond_selections()
        
        # Clear vector selections (if any)
        for vector in self.viewer.ortep_mol.vectors:
            vector.selected = False

    def _toggle_help(self):
        """
        Toggle the help overlay on/off in the viewer, and also
        print the cheat sheet to stdout if turning it on.
        """
        self.viewer.show_help = not self.viewer.show_help
        if self.viewer.show_help:
            # Print help text to console as well
            from help_text import HELP_TEXT
            for line in HELP_TEXT:
                print(line)

        # Force a redraw so the overlay becomes visible or hidden
        self.viewer.redraw()

