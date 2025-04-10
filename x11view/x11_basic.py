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

    def __init__(self, x11_window, background_color=(255, 255, 255)):
        """
        Initialize the basic canvas. Allocates an offscreen Pixmap
        to match the current window size.
        
        Parameters:
            x11_window: The parent X11Window instance
            background_color: RGB tuple for the background color (default: white)
        """
        # Initialize abstract base + common mixin
        X11CanvasBase.__init__(self, x11_window)
        X11CanvasCommon.__init__(self)
        
        # Store the background color
        self.background_color = background_color

        # Create our offscreen Pixmap
        self.pixmap = self.x11_window.window.create_pixmap(
            depth=self.screen.root_depth,
            width=self.width,
            height=self.height
        )
        
        # Initialize font cache
        self.font_cache = {}
        
        # Pre-load commonly used fonts
        self.default_font = None
        common_fonts = [
            "6x13",  # Smaller font, common on X11 systems
            "8x13",  # Alternative smaller font
            "-misc-fixed-medium-r-normal--10-100-75-75-c-60-iso8859-1",  # Smaller misc-fixed
            "9x15",  # Larger font as fallback
            "-misc-fixed-medium-r-normal--13-120-75-75-c-70-iso8859-1"  # Final fallback
        ]
        
        # Try to load at least one default font
        for font_name in common_fonts:
            try:
                self.default_font = self.display.open_font(font_name)
                break
            except (xerror.BadName, xerror.BadFont):
                pass

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
        Fill the offscreen Pixmap with the background color.
        """
        bg_gc = self.get_gc(
            color=self.background_color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        self.pixmap.poly_fill_rectangle(bg_gc, [(0, 0, self.width, self.height)])

    def set_background_color(self, color):
        """
        Set the background color of the canvas.
        
        Parameters:
            color: RGB tuple (r, g, b) with values in range 0-255
        """
        self.background_color = color

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
    # Rectangle border
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
    # Filled rectangle
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
    # Triangle border
    # -----------------------------------------------------------
    def draw_triangle(self,
                     x1: int, y1: int,
                     x2: int, y2: int,
                     x3: int, y3: int,
                     color=(0, 0, 0),
                     thickness=2) -> None:
        """
        Draw an unfilled triangle by connecting three points with lines.
        """
        gc = self.get_gc(
            color=color,
            thickness=thickness,
            line_style=X.LineSolid,
            fill_style=False
        )
        
        # Create a closed path using X11's poly_line
        self.pixmap.poly_line(
            gc,
            X.CoordModeOrigin,
            [(x1, y1), (x2, y2), (x3, y3), (x1, y1)]  # Close the path by repeating first point
        )

    # -----------------------------------------------------------
    # Filled triangle
    # -----------------------------------------------------------
    def draw_filled_triangle(self,
                            x1: int, y1: int,
                            x2: int, y2: int,
                            x3: int, y3: int,
                            color=(0, 0, 0)) -> None:
        """
        Draw a filled triangle using X11's fill_poly.
        """
        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        
        # Use fill_poly to create a filled triangle
        # FIXED: Correct parameter order - shape comes before coord_mode in Python-Xlib
        self.pixmap.fill_poly(
            gc,
            X.Convex,           # Shape (Complex, Nonconvex, Convex)
            X.CoordModeOrigin,  # Coordinates are absolute
            [(x1, y1), (x2, y2), (x3, y3)]
        )


    # -----------------------------------------------------------
    # Improved: Text drawing with font caching
    # -----------------------------------------------------------

    def draw_text(self, x, y, text, color=(0, 0, 0), font_size=12,
                  font_candidates=None):
        """
        Draw text so that (x, y) is the text's baseline in X11.
        Uses font caching to avoid repeated X server round-trips.
        If font_candidates is provided, it should be a list of font strings
        (core X11 font names). We try each in order until one works.
        If none work, we use a hardcoded fallback or pre-loaded default.
        """
        gc = self.get_gc(
            color=color,
            thickness=1,
            line_style=X.LineSolid,
            fill_style=False
        )
    
        # If we have a pre-loaded default font and no specific candidates requested,
        # use the default font directly
        if not font_candidates and self.default_font is not None:
            gc.change(font=self.default_font.id)
            self.pixmap.draw_text(gc, x, y, text)
            return
    
        # Otherwise, proceed with font selection but use caching
        fallback_font_str = "-misc-fixed-medium-r-normal--13-120-75-75-c-70-iso8859-1"
        if not font_candidates:
            font_candidates = [fallback_font_str]
    
        # Create a hashable cache key
        cache_key = tuple(font_candidates)  # Font size isn't used in X11 direct mode
    
        # Check if we already have this font cached
        if cache_key in self.font_cache:
            font_obj = self.font_cache[cache_key]
            if font_obj is not None:
                gc.change(font=font_obj.id)
        else:
            # Not cached - need to try loading fonts
            font_obj = None
            for font_str in font_candidates:
                try:
                    font_obj = self.display.open_font(font_str)
                    break
                except (xerror.BadName, xerror.BadFont):
                    # Could not open this font name
                    pass
    
            if font_obj is None and self.default_font is not None:
                # Use pre-loaded default if all candidates failed
                font_obj = self.default_font
            elif font_obj is None:
                # Try the fallback if no default is available
                try:
                    font_obj = self.display.open_font(fallback_font_str)
                except:
                    # If we fail here as well, proceed with whatever X11 defaults to
                    pass
    
            # Cache the result (even if None, to avoid repeated attempts)
            self.font_cache[cache_key] = font_obj
    
            # Set the font if we found one
            if font_obj is not None:
                gc.change(font=font_obj.id)
    
        # Draw the text (will use whatever font we successfully set, or X11's default)
        self.pixmap.draw_text(gc, x, y, text)
        
    def get_text_dimensions(self, text, font_size=12):
        """
        Estimate the dimensions of text in pixels.
        In X11 mode, this is an approximation as we don't query the actual font metrics.
        
        Parameters:
            text: The text string to measure
            font_size: The font size (not directly used in X11 core fonts)
            
        Returns:
            (width, height) tuple in pixels
        """
        # For X11 core fonts, we don't have direct access to font metrics in the Python binding
        # Instead, we'll use a simple approximation based on typical fixed-width fonts
        
        # For our default fonts, approximate character width (most are fixed-width)
        char_width = 8  # Default estimate for a typical fixed-width X11 font
        width = len(text) * char_width
        
        # Height is typically the font size or slightly larger
        height = 13  # Typical height for misc-fixed fonts
        
        return (width, height)
