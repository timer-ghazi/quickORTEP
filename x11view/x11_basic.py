# x11_basic.py

"""
x11view.x11_basic

Implements a basic, double-buffered X11 canvas using a Pixmap.
No supersampling or external libraries are required.
"""

from Xlib import X
from Xlib import error as xerror

from .base import X11CanvasBase
from .common import X11CanvasCommon



class X11CanvasBasic(X11CanvasBase, X11CanvasCommon):
    """
    A direct X11-based canvas that uses a Pixmap for double buffering.
    All drawing operations are performed on the offscreen Pixmap,
    then copied to the real window in flush().
    """

    def __init__(self, x11_window):
        """
        Initialize the basic canvas. Allocates an offscreen Pixmap
        to match the current window size.
        """
        # Initialize abstract base + common mixin
        X11CanvasBase.__init__(self, x11_window)
        X11CanvasCommon.__init__(self)

        # Create our offscreen Pixmap
        self.pixmap = self.x11_window.window.create_pixmap(
            depth=self.screen.root_depth,
            width=self.width,
            height=self.height
        )

    def create_or_resize(self, width: int, height: int) -> None:
        """
        Handle the window resize event. Recreate our Pixmap to match
        the new size and clear our GC cache.
        """
        self.width = width
        self.height = height

        self.pixmap = self.x11_window.window.create_pixmap(
            depth=self.screen.root_depth,
            width=self.width,
            height=self.height
        )

        self.gc_cache.clear()

    def clear(self) -> None:
        """
        Fill the offscreen Pixmap with white.
        """
        white_gc = self.get_gc(
            color=(255, 255, 255),
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        self.pixmap.poly_fill_rectangle(white_gc, [(0, 0, self.width, self.height)])

    def flush(self) -> None:
        """
        Copy the offscreen Pixmap to the real X11 window and flush.
        """
        temp_gc = self.x11_window.window.create_gc()

        self.x11_window.window.copy_area(
            gc=temp_gc,
            src_drawable=self.pixmap,
            src_x=0,
            src_y=0,
            width=self.width,
            height=self.height,
            dst_x=0,
            dst_y=0
        )

        self.display.flush()

    # ------------------------------------------------------------------
    # Implement the required drawing methods from X11CanvasBase
    # ------------------------------------------------------------------

    def draw_filled_circle(self,
                           cx: int,
                           cy: int,
                           radius: int,
                           color=(255, 0, 0)
                           ) -> None:
        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        self.pixmap.fill_arc(
            gc,
            cx - radius,
            cy - radius,
            radius * 2,
            radius * 2,
            0,
            360 * 64
        )

    def draw_circle_border(self,
                           cx: int,
                           cy: int,
                           radius: int,
                           color=(0, 0, 0),
                           thickness=2
                           ) -> None:
        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineSolid,
            fill_style=False
        )
        self.pixmap.poly_arc(
            gc,
            [(cx - radius, cy - radius, radius * 2, radius * 2, 0, 360 * 64)]
        )

    def draw_line(self,
                  x1: int,
                  y1: int,
                  x2: int,
                  y2: int,
                  thickness: int = 4,
                  color=(0, 0, 0)
                  ) -> None:
        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineSolid,
            fill_style=False
        )
        self.pixmap.poly_line(
            gc,
            X.CoordModeOrigin,
            [(x1, y1), (x2, y2)]
        )

    def draw_dashed_line(self,
                         x1: int,
                         y1: int,
                         x2: int,
                         y2: int,
                         thickness: int = 2,
                         color=(128, 128, 128)
                         ) -> None:
        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineOnOffDash,
            fill_style=False
        )
        # Customize dash pattern if desired
        gc.set_dashes(0, [8, 4])  # 8px on, 4px off

        self.pixmap.poly_line(
            gc,
            X.CoordModeOrigin,
            [(x1, y1), (x2, y2)]
        )

    def draw_arc(self,
                 cx: int,
                 cy: int,
                 rx: int,
                 ry: int = None,
                 angle_start_deg: float = 0.0,
                 angle_end_deg: float = 360.0,
                 thickness: int = 2,
                 color=(0, 0, 0)
                 ) -> None:
        if ry is None:
            ry = rx

        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineSolid,
            fill_style=False
        )
        start_64 = int(angle_start_deg * 64)
        extent_64 = int((angle_end_deg - angle_start_deg) * 64)

        self.pixmap.poly_arc(
            gc,
            [(cx - rx, cy - ry, rx * 2, ry * 2, start_64, extent_64)]
        )

    def draw_filled_arc(self,
                        cx: int,
                        cy: int,
                        rx: int,
                        ry: int = None,
                        angle_start_deg: float = 0.0,
                        angle_end_deg: float = 360.0,
                        color=(0, 0, 0)
                        ) -> None:
        if ry is None:
            ry = rx

        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        start_64 = int(angle_start_deg * 64)
        extent_64 = int((angle_end_deg - angle_start_deg) * 64)

        self.pixmap.fill_arc(
            gc,
            cx - rx,
            cy - ry,
            rx * 2,
            ry * 2,
            start_64,
            extent_64
        )


    # -----------------------------------------------------------
    # New: Rectangle border
    # -----------------------------------------------------------
    def draw_rect(self, x, y, width, height, color=(0, 0, 0), thickness=2):
        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineSolid,
            fill_style=False
        )
        # Draws the rectangle border in one go
        self.pixmap.rectangle(gc, x, y, width, height)

    # -----------------------------------------------------------
    # New: Filled rectangle
    # -----------------------------------------------------------
    def draw_filled_rect(self, x, y, width, height, color=(0, 0, 0)):
        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        # Fill a single rectangle
        self.pixmap.poly_fill_rectangle(gc, [(x, y, width, height)])

    # -----------------------------------------------------------
    # New: Text drawing
    # -----------------------------------------------------------

    def draw_text(self, x, y, text, color=(0, 0, 0), font_size=12,
                  font_candidates=None):
        """
        Draw text so that (x, y) is the text's baseline in X11.
        If font_candidates is provided, it should be a list of font strings
        (core X11 font names). We try each in order until one works.
        If none work, we use a hardcoded fallback.
        """
        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=False
        )
    
        fallback_font_str = "-misc-fixed-medium-r-normal--13-120-75-75-c-70-iso8859-1"
        if not font_candidates:
            font_candidates = [fallback_font_str]
    
        chosen_font_str = None
        chosen_font_obj = None
    
        for font_str in font_candidates:
            try:
                font_obj = self.display.open_font(font_str)
                chosen_font_str = font_str
                chosen_font_obj = font_obj
                break
            except (xerror.BadName, xerror.BadFont):
                # Could not open this font name
                pass
    
        if chosen_font_obj is None:
            # Attempt the fallback
            try:
                chosen_font_obj = self.display.open_font(fallback_font_str)
            except:
                # If we fail here as well, nothing left to do; proceed with no font
                pass
    
        if chosen_font_obj is not None:
            # **Extract the numeric ID** from the resource object
            gc.change(font=chosen_font_obj.id)
    
        # If we reach here with no loaded font, the text call might do nothing or fail.
        # If you want to guard, you can check chosen_font_obj is not None.
    
        self.pixmap.draw_text(gc, x, y, text)
    
