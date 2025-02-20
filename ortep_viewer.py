import sys
from Xlib import XK

from x11view.window import X11Window, X11CanvasBasic, X11CanvasSS
from ortep_renderer import ORTEP_MoleculeRenderer

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

class MoleculeViewer(X11Window):
    """
    The 'controller/view' class that manages an ORTEP_Molecule, an ORTEP_MoleculeRenderer,
    and handles user events (key presses), etc.
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


    def redraw(self):
        self.canvas.clear()
        # remove old bond/atom calls
        # self.renderer._draw_bonds(...)
        # self.renderer._draw_atoms(...)
        # now we do:
        self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
        self.canvas.flush()

#     def redraw(self):
#         self.canvas.clear()
#        #-----  # self.renderer.draw(self.canvas, self.ortep_mol, self.view_params)
#        #-----  # Instead, try the advanced partial-occlusion approach:
#        #-----  self.renderer.draw_segmented(self.canvas, self.ortep_mol, self.view_params)
# 
#         self.renderer.draw_molecule(self.canvas, self.ortep_mol, self.view_params)
#         self.canvas.flush()

    def handle_key(self, evt):
        # same as in the old code
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
        else:
            print(f"Ignored key: {keychar}")

        self.redraw()
