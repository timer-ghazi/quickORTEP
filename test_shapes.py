#!/usr/bin/env python3
"""
test_shapes.py

A dedicated test script that showcases a variety of shapes and angles
to visually inspect anti-aliasing. First, it opens a window with the
basic (non-SS) canvas. When that window is closed, it opens another
window with the supersampled (SS) canvas.

Press 'q' or 'Escape' in either window to close it.
"""

from x11view.window import X11Window
from x11view.x11_basic import X11CanvasBasic
from x11view.x11_ss import X11CanvasSS

import math


class ShapeTestWindow(X11Window):
    """
    A custom window that overrides `redraw()` to draw
    a variety of angled shapes for visual inspection.
    """

    def redraw(self):
        """
        Draw a selection of shapes at different angles, arcs, etc.
        """
        self.canvas.clear()

        w = self.canvas.width
        h = self.canvas.height
        cx = w // 2
        cy = h // 2

        # 1) Draw radial lines from center at multiple angles
        radius = min(cx, cy) - 20
        num_spokes = 12
        for i in range(num_spokes):
            angle_deg = (360 / num_spokes) * i
            angle_rad = math.radians(angle_deg)
            x2 = cx + int(radius * math.cos(angle_rad))
            y2 = cy + int(radius * math.sin(angle_rad))
            self.canvas.draw_line(cx, cy, x2, y2, thickness=2, color=(0, 0, 0))

        # 2) Draw arcs at various angles in a ring around the center
        arc_radius = 40
        step_deg = 45
        for angle_start in range(0, 360, step_deg):
            angle_end = angle_start + 30
            arc_cx = cx + int((radius // 2) * math.cos(math.radians(angle_start)))
            arc_cy = cy - int((radius // 2) * math.sin(math.radians(angle_start)))
            # arc_cx = cx + int((radius // 2) * math.cos(math.radians(angle_start)))
            # arc_cy = cy + int((radius // 2) * math.sin(math.radians(angle_start)))
            self.canvas.draw_arc(
                arc_cx,
                arc_cy,
                arc_radius,
                ry=arc_radius // 2,
                angle_start_deg=angle_start,
                angle_end_deg=angle_end,
                thickness=3,
                color=(0, 0, 255)
            )

        # 3) Filled arcs in corners
        corners = [(50, 50), (w - 50, 50), (50, h - 50), (w - 50, h - 50)]
        corner_angles = [(180, 270), (0, 90), (90, 180), (270, 360)]
        for (corner_x, corner_y), (start_a, end_a) in zip(corners, corner_angles):
            self.canvas.draw_filled_arc(
                corner_x,
                corner_y,
                30,
                ry=20,
                angle_start_deg=start_a,
                angle_end_deg=end_a,
                color=(255, 128, 0)
            )

        # 4) A ring of small circles around the center
        circle_count = 8
        circle_radius = 15
        circle_ring_r = 100
        for i in range(circle_count):
            angle_deg = 360 * i / circle_count
            rad = math.radians(angle_deg)
            cx2 = cx + int(circle_ring_r * math.cos(rad))
            cy2 = cy + int(circle_ring_r * math.sin(rad))
            self.canvas.draw_filled_circle(cx2, cy2, circle_radius, color=(0, 200, 0))
            self.canvas.draw_circle_border(cx2, cy2, circle_radius, color=(0, 0, 0), thickness=2)

        # 5) A dashed line diagonal across the window
        self.canvas.draw_dashed_line(0, 0, w, h, thickness=3, color=(128, 128, 128))

        # Flush to show everything
        self.canvas.flush()


def main():
    """
    Run a two-phase test:
      1) Show the shapes in a window with Basic canvas.
      2) After it's closed, show the shapes in a window with SS canvas
         and an artificially large tile_size to test clamping.
    """
    # 1) Basic (non-SS) version
    print("Opening Basic (non-SS) window...")
    basic_win = ShapeTestWindow(width=800, height=600, title="Shape Test - Basic",
                                canvas_class=X11CanvasBasic)
    basic_win.run()  # Blocks until user closes or hits 'q'

    # 2) SS version with an over-the-top tile_size to test clamping
    print("Opening SS window (factor=2) with very large tile_size=99999 ...")
    def ss_factory(x11_window):
        return X11CanvasSS(x11_window, ss_factor=2, tile_size=99999)

    ss_win = ShapeTestWindow(width=800, height=600, title="Shape Test - SS",
                             canvas_class=ss_factory)
    ss_win.run()


if __name__ == "__main__":
    main()
