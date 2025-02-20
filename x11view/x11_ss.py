# x11_ss.py

"""
x11view.x11_ss

Implements a supersampled X11 canvas using Pillow (PIL).
All drawing calls are directed to a high-resolution PIL image,
which is then downsampled and copied (in tiles) to an X11 Pixmap,
and finally displayed in the window.
"""

import sys
import numpy as np
from PIL import Image, ImageDraw
from Xlib import X

from .base import X11CanvasBase
from .common import X11CanvasCommon


class X11CanvasSS(X11CanvasBase, X11CanvasCommon):
    """
    A supersampled X11 canvas. Drawing is done in a high-resolution
    (ss_width x ss_height) Pillow image, then downsampled to the
    actual window size for smoother edges (antialiasing).

    The final image is copied to the real window using put_image().
    Tiling is used to avoid sending large images all at once if
    tile_size is smaller than the window dimensions.
    """

    def __init__(self,
                 x11_window,
                 ss_factor: int = 2,
                 tile_size: int = 128):
        """
        :param x11_window: The parent X11Window object
        :param ss_factor:  Supersampling factor (e.g., 2 => 2x in each dimension)
        :param tile_size:  Desired tile size in pixels for put_image().
        """
        # Initialize parent classes
        X11CanvasBase.__init__(self, x11_window)
        X11CanvasCommon.__init__(self)

        self.ss_factor = ss_factor

        # 1) Get max request length from python-xlib's .info object
        max_req_length_units = self.display.display.info.max_request_length


        # 2) Convert to bytes, subtract overhead
        overhead = 4096  # arbitrary safety margin
        max_bytes = (max_req_length_units * 4) - overhead

        # 3) Each BGRA pixel is 4 bytes
        max_pixels = max_bytes // 4

        import math
        max_tile_side = int(math.isqrt(max_pixels) * 0.95)
        if max_tile_side < 1:
            max_tile_side = 1

        # 4) Clamp tile_size to prevent exceeding the server's max
        tile_size = min(tile_size, max_tile_side)
        tile_size = max(tile_size, 1)
        self.tile_size = tile_size

        # Calculate supersampled dimensions
        self.ss_width = self.width * self.ss_factor
        self.ss_height = self.height * self.ss_factor

        # Create an offscreen Pixmap for the final (downsampled) image
        self.pixmap = self.x11_window.window.create_pixmap(
            depth=self.screen.root_depth,
            width=self.width,
            height=self.height
        )

        # Create a high-resolution PIL Image and Draw object
        self.ss_image = Image.new("RGB", (self.ss_width, self.ss_height), "white")
        self.ss_draw = ImageDraw.Draw(self.ss_image)

    def create_or_resize(self, width: int, height: int) -> None:
        """
        Called when the window is resized. Recreate the X11 Pixmap
        and the high-resolution PIL image.
        """
        self.width = width
        self.height = height

        # Update supersampled dimensions
        self.ss_width = width * self.ss_factor
        self.ss_height = height * self.ss_factor

        # Clear GC cache
        self.gc_cache.clear()

        # Recreate X11 Pixmap for final display
        self.pixmap = self.x11_window.window.create_pixmap(
            depth=self.screen.root_depth,
            width=width,
            height=height
        )

        # Recreate the high-res PIL image
        self.ss_image = Image.new("RGB", (self.ss_width, self.ss_height), "white")
        self.ss_draw = ImageDraw.Draw(self.ss_image)

    def clear(self) -> None:
        """
        Clear both the high-resolution PIL buffer and the X11 Pixmap.
        """
        # Clear the high-res image
        self.ss_image = Image.new("RGB", (self.ss_width, self.ss_height), "white")
        self.ss_draw = ImageDraw.Draw(self.ss_image)

        # Clear the X11 Pixmap as well
        white_gc = self.get_gc(
            color=(255, 255, 255),
            thickness=1,
            line_style=X.LineSolid,
            fill_style=True
        )
        self.pixmap.poly_fill_rectangle(white_gc, [(0, 0, self.width, self.height)])

    def flush(self) -> None:
        """
        Downsample the high-resolution image, then copy it
        to the X11 Pixmap in tiles, and finally copy that Pixmap
        to the window.
        """
        # 1) Downsample using a high-quality filter
        final_img = self.ss_image.resize(
            (self.width, self.height),
            Image.Resampling.LANCZOS
        )

        # 2) Convert to BGRA
        rgb_array = np.array(final_img)       # shape: (height, width, 3)
        bgr_array = rgb_array[:, :, ::-1]     # reorder channels from RGB to BGR
        alpha_channel = np.full((self.height, self.width, 1), 255, dtype=np.uint8)
        bgra_array = np.concatenate([bgr_array, alpha_channel], axis=2)  # shape: (h, w, 4)

        # 3) Put image in tiles
        gc = self.get_gc(color=(0, 0, 0),
                         thickness=1,
                         line_style=X.LineSolid,
                         fill_style=False)

        for ystart in range(0, self.height, self.tile_size):
            for xstart in range(0, self.width, self.tile_size):
                xend = min(xstart + self.tile_size, self.width)
                yend = min(ystart + self.tile_size, self.height)

                tile_width = xend - xstart
                tile_height = yend - ystart

                tile_data = bgra_array[ystart:yend, xstart:xend]
                tile_bytes = tile_data.tobytes()

                self.pixmap.put_image(
                    gc,
                    xstart,  # where in the Pixmap to place
                    ystart,
                    tile_width,
                    tile_height,
                    X.ZPixmap,
                    self.screen.root_depth,
                    0,
                    tile_bytes
                )

        # 4) Copy the Pixmap to the actual window
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

    # -----------------------------------------------------------
    # Drawing methods (in high-res space). We multiply coords by
    # ss_factor before using Pillow's drawing commands.
    # -----------------------------------------------------------

    def draw_filled_circle(self, cx: int, cy: int, radius: int,
                           color=(255, 0, 0)) -> None:
        ss_cx = cx * self.ss_factor
        ss_cy = cy * self.ss_factor
        ss_radius = radius * self.ss_factor
        self.ss_draw.ellipse(
            [(ss_cx - ss_radius, ss_cy - ss_radius),
             (ss_cx + ss_radius, ss_cy + ss_radius)],
            fill=color
        )

    def draw_circle_border(self, cx: int, cy: int, radius: int,
                           color=(0, 0, 0), thickness=2) -> None:
        ss_cx = cx * self.ss_factor
        ss_cy = cy * self.ss_factor
        ss_radius = radius * self.ss_factor
        ss_thick = thickness * self.ss_factor
        self.ss_draw.ellipse(
            [(ss_cx - ss_radius, ss_cy - ss_radius),
             (ss_cx + ss_radius, ss_cy + ss_radius)],
            outline=color,
            width=ss_thick
        )

    def draw_line(self, x1: int, y1: int, x2: int, y2: int,
                  thickness: int = 4, color=(0, 0, 0)) -> None:
        ss_x1 = x1 * self.ss_factor
        ss_y1 = y1 * self.ss_factor
        ss_x2 = x2 * self.ss_factor
        ss_y2 = y2 * self.ss_factor
        ss_thick = thickness * self.ss_factor
        self.ss_draw.line([(ss_x1, ss_y1), (ss_x2, ss_y2)],
                          fill=color,
                          width=ss_thick)

    def draw_dashed_line(self, x1: int, y1: int, x2: int, y2: int,
                         thickness: int = 2, color=(128, 128, 128)) -> None:
        """
        Pillow doesn't directly support dashed lines, so we approximate
        by drawing small line segments separated by gaps.
        """
        import math
        ss_x1 = x1 * self.ss_factor
        ss_y1 = y1 * self.ss_factor
        ss_x2 = x2 * self.ss_factor
        ss_y2 = y2 * self.ss_factor
        ss_thick = thickness * self.ss_factor

        dash_len = 8 * self.ss_factor
        gap_len = 4 * self.ss_factor

        dx = ss_x2 - ss_x1
        dy = ss_y2 - ss_y1
        seg_length = math.hypot(dx, dy)
        angle = math.atan2(dy, dx)

        step = dash_len + gap_len
        n_steps = int(seg_length // step)
        cur_dist = 0.0

        for _ in range(n_steps + 1):
            start_dist = cur_dist
            end_dist = min(cur_dist + dash_len, seg_length)
            if start_dist > seg_length:
                break
            sx = ss_x1 + math.cos(angle) * start_dist
            sy = ss_y1 + math.sin(angle) * start_dist
            ex = ss_x1 + math.cos(angle) * end_dist
            ey = ss_y1 + math.sin(angle) * end_dist
            self.ss_draw.line([(sx, sy), (ex, ey)],
                              fill=color,
                              width=ss_thick)
            cur_dist += step

    def draw_arc(self, cx, cy, rx, ry=None,
                 angle_start_deg=0.0, angle_end_deg=360.0,
                 thickness=2, color=(0, 0, 0)):
        import math
        if ry is None:
            ry = rx
    
        ss_cx = cx * self.ss_factor
        ss_cy = cy * self.ss_factor
        ss_rx = rx * self.ss_factor
        ss_ry = ry * self.ss_factor
        ss_thick = thickness * self.ss_factor
    
        # We take the difference so it can be negative for clockwise arcs
        total_deg = angle_end_deg - angle_start_deg
    
        # Number of segments for the polyline
        n_points = max(12, int(abs(total_deg))) * 2
    
        start_deg = angle_start_deg
        # We'll sample from start_deg to start_deg + total_deg
        # This may be up or down in angle if total_deg < 0
    
        points = []
        for i in range(n_points + 1):
            frac = i / n_points
            a_deg = start_deg + total_deg * frac
            a_rad = math.radians(a_deg)
            # Xlib style: 0° => 3 o'clock, angles CCW => 90 => 12 o'clock
            # in a coordinate system with y increasing down => we do "cy - sin(...)"
            x_ = ss_cx + ss_rx * math.cos(a_rad)
            y_ = ss_cy - ss_ry * math.sin(a_rad)
            points.append((x_, y_))
    
        # Draw the polyline that approximates the arc
        self.ss_draw.line(points, fill=color, width=ss_thick)

    def draw_filled_arc(self, cx, cy, rx, ry=None,
                        angle_start_deg=0.0, angle_end_deg=360.0,
                        color=(0, 0, 0)):
        if ry is None:
            ry = rx
    
        # Calculate the “Xlib” extent
        extent_deg = angle_end_deg - angle_start_deg
        # If it is negative, add 360 so we do a CCW wedge
        if extent_deg < 0:
            extent_deg += 360.0
    
        pillow_start = angle_start_deg
        pillow_end   = angle_start_deg + extent_deg
    
        ss_cx = cx * self.ss_factor
        ss_cy = cy * self.ss_factor
        ss_rx = rx * self.ss_factor
        ss_ry = ry * self.ss_factor
    
        bbox = [ss_cx - ss_rx, ss_cy - ss_ry,
                ss_cx + ss_rx, ss_cy + ss_ry]
    
        self.ss_draw.pieslice(
            bbox,
            start=pillow_start,
            end=pillow_end,
            fill=color
        )

