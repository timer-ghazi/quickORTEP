# ortep_viewer.py

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
    The 'controller/view' class that manages an ORTEP_Molecule, an ORTEP_MoleculeRenderer,
    and handles user events (key presses, mouse events, etc.).
    """

    def __init__(self, ortep_molecule, width=800, height=600,
                 ss_factor=1, tile_size=128):
        """
        :param ss_factor: integer supersampling factor (1 = no AA)
        :param tile_size: tile size if using X11CanvasSS
        """
        # Choose the canvas class
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

        # Initialize mouse state variables
        self.active_button = None         # 'left', 'middle', 'right', or 'shift-left'
        self.last_mouse_x = None
        self.last_mouse_y = None
        self.last_click_time = 0.0

    def redraw(self):
        self.canvas.clear()
        self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
        self.canvas.flush()

    def fit_molecule_to_window(self):

        # Collect rotated x,y coordinates of all atoms.
        xs, ys = [], []
        for atom in self.ortep_mol.atoms:
            # Apply current rotation.
            x_rot, y_rot, _ = rotate_point(atom.x, atom.y, atom.z,
                                           self.view_params.rx,
                                           self.view_params.ry,
                                           self.view_params.rz)
            xs.append(x_rot)
            ys.append(y_rot)
        
        if not xs or not ys:
            return  # nothing to do if no atoms
    
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        
        # Define a margin (in pixels)
        margin = 20
        available_width = self.canvas.width - 2 * margin
        available_height = self.canvas.height - 2 * margin
        
        # Avoid division by zero if extents are zero
        extent_x = max_x - min_x if max_x > min_x else 1.0
        extent_y = max_y - min_y if max_y > min_y else 1.0
        
        new_scale = min(available_width / extent_x, available_height / extent_y)
        
        # Update scale and offsets
        self.view_params.scale = new_scale
        center_x = (min_x + max_x) / 2.0
        center_y = (min_y + max_y) / 2.0
        self.view_params.x_offset = self.canvas.width / 2 - center_x * new_scale
        self.view_params.y_offset = self.canvas.height / 2 - center_y * new_scale
        
        self.redraw()

    def handle_key(self, evt):
        # same as in the old code for key events
        keysym = self.display.keycode_to_keysym(evt.detail, evt.state)
        keychar = XK.keysym_to_string(keysym)

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
        # evt.detail: 1=left, 2=middle, 3=right, 4=scroll up, 5=scroll down

        if evt.detail in (4, 5):
            # Handle zoom events
            if evt.detail == 4:
                self.view_params.scale *= 1.1
            elif evt.detail == 5:
                self.view_params.scale /= 1.1
            self.redraw()
            return

        self.last_mouse_x = evt.event_x
        self.last_mouse_y = evt.event_y

        # Check for Shift modifier (using evt.state & X.ShiftMask)
        shift_pressed = bool(evt.state & X.ShiftMask)

        if evt.detail == 1:
            # For left-click, check for double-click reset if Shift is pressed.
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

        # Define rotation factors (in degrees per pixel)
        k_x = 0.5  # vertical movement -> rotation around X
        k_y = 0.5  # horizontal movement -> rotation around Y
        k_z = 0.5  # horizontal movement -> rotation around Z

        if self.active_button == 'left':
            # Left-drag: rotate around X (vertical) and Y (horizontal)
            self.view_params.rx += k_x * dy
            self.view_params.ry += k_y * dx
        elif self.active_button in ('middle', 'shift-left'):
            # Middle or Shift+left drag: rotate around Z using horizontal movement
            self.view_params.rz += k_z * dx
        elif self.active_button == 'right':
            # Right-drag: pan the view (translate)
            self.view_params.x_offset += dx
            self.view_params.y_offset += dy

        self.redraw()

    def handle_button_release(self, evt):
        # Reset the active button when the button is released.
        self.active_button = None

    def dump_svg(self):
        # Create an SVG canvas
        svg_canvas = SVGCanvas(width=self.canvas.width, height=self.canvas.height)
    
        # Temporarily tell the renderer we want float coords
        self.view_params.as_float = True
    
        # Render to the SVG canvas
        self.renderer.draw_molecule(svg_canvas, self.ortep_mol, self.view_params)
    
        # Turn it off again so normal X11 drawing stays integer-based
        self.view_params.as_float = False
    
        # Write the SVG content
        filename = "quickORTEP_export.svg"
        svg_canvas.flush(filename)
        # print(f"SVG export saved to {filename}")


    def reset_view(self):
        # Reset view parameters to their defaults.
        self.view_params.rx = 0.0
        self.view_params.ry = 0.0
        self.view_params.rz = 0.0
        self.view_params.scale = 100.0
        self.view_params.x_offset = self.canvas.width // 2
        self.view_params.y_offset = self.canvas.height // 2
        self.redraw()
