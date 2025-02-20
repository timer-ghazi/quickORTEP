# common.py

"""
x11view.common

A mixin/utility class that provides commonly used functions
for X11-based canvases, such as GC (graphics context) caching
and color conversion (RGB -> 24-bit pixel).
"""

from Xlib import X
from typing import Dict, Tuple

class X11CanvasCommon:
    """
    A small mixin/utility class that can be inherited by
    X11CanvasBase subclasses to provide GC caching and
    a rgb_to_pixel() helper.
    """

    def __init__(self):
        # Dictionary to cache Graphics Contexts (GCs) by parameters
        # Key: (color, thickness, line_style, fill_style)
        self.gc_cache: Dict[Tuple, "X.GC"] = {}

    def get_gc(self,
               color: Tuple[int, int, int],
               thickness: int,
               line_style: int,
               fill_style: bool = False
               ) -> "X.GC":
        """
        Return a cached GC for the given parameters if available;
        otherwise create a new GC and store it in the cache.

        :param color:      (R, G, B) tuple in 0..255
        :param thickness:  Line width in pixels
        :param line_style: e.g. X.LineSolid or X.LineOnOffDash
        :param fill_style: True if this GC will be used for filled shapes
        :return:           A Graphics Context (GC) object
        """
        key = (color, thickness, line_style, fill_style)

        if key not in self.gc_cache:
            fg_pixel = self.rgb_to_pixel(color)
            if fill_style:
                # For fills, thickness/line_style don't matter
                gc = self.x11_window.window.create_gc(
                    foreground=fg_pixel,
                    background=self.screen.white_pixel
                )
            else:
                gc = self.x11_window.window.create_gc(
                    foreground=fg_pixel,
                    background=self.screen.white_pixel,
                    line_width=thickness,
                    line_style=line_style,
                    cap_style=X.CapRound,
                    join_style=X.JoinRound
                )
            self.gc_cache[key] = gc

        return self.gc_cache[key]

    @staticmethod
    def rgb_to_pixel(rgb: Tuple[int, int, int]) -> int:
        """
        Convert (R,G,B) in [0..255] to a 24-bit integer pixel value.

        :param rgb: (R,G,B) tuple
        :return:     24-bit integer (0xRRGGBB)
        """
        r, g, b = rgb
        return (r << 16) | (g << 8) | b
