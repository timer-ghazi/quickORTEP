#!/usr/bin/env python3
"""
test_shapes.py

A dedicated test script that showcases a variety of shapes and angles
in a basic canvas window. After the window is closed, it also saves an 
SVG file ("test_shapes_export.svg") that contains the same shapes, 
drawn with our SVGCanvas backend.

Press 'q' or 'Escape' to close the window.
"""

import math

from x11view.window import X11Window
from x11view.x11_basic import X11CanvasBasic

# Import our SVGCanvas
from x11view.svg_canvas import SVGCanvas


def draw_test_shapes(canvas, w, h):
    """
    Draw the same shapes that we show in the X11 windows, but on the
    specified canvas (which could be X11 or SVG). We pass 'w' and 'h'
    so the code doesn't rely on canvas.width or canvas.height directly.
    """
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
        canvas.draw_line(cx, cy, x2, y2, thickness=2, color=(0, 0, 0))

    # 2) Draw arcs at various angles in a ring around the center
    arc_radius = 40
    step_deg = 45
    for angle_start in range(0, 360, step_deg):
        angle_end = angle_start + 30
        arc_cx = cx + int((radius // 2) * math.cos(math.radians(angle_start)))
        arc_cy = cy - int((radius // 2) * math.sin(math.radians(angle_start)))
        canvas.draw_arc(
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
        canvas.draw_filled_arc(
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
        canvas.draw_filled_circle(cx2, cy2, circle_radius, color=(0, 200, 0))
        canvas.draw_circle_border(cx2, cy2, circle_radius, color=(0, 0, 0), thickness=2)

    # 5) A dashed line diagonal across the window
    canvas.draw_dashed_line(0, 0, w, h, thickness=3, color=(128, 128, 128))


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

        # Draw shapes onto our X11 canvas
        draw_test_shapes(self.canvas, w, h)

        # Flush to show everything
        self.canvas.flush()


def main():
    """
    Run a simple test:
      1) Show the shapes in a window with Basic canvas.
      2) After it's closed, export an SVG file with the same shapes for offline inspection.
    """
    # 1) Basic window
    print("Opening Basic window...")
    basic_win = ShapeTestWindow(width=800, height=600, title="Shape Test")
    basic_win.run()  # Blocks until user closes or hits 'q'

    # 2) Generate an SVG version of the same shapes
    print("Now saving an SVG version of the same shapes to 'test_shapes_export.svg' ...")
    svg_canvas = SVGCanvas(width=800, height=600)
    draw_test_shapes(svg_canvas, 800, 600)
    svg_canvas.flush("test_shapes_export.svg")


if __name__ == "__main__":
    main()
