# x11view/svg_canvas.py

import math

class SVGCanvas:
    """
    SVGCanvas implements a drawing API similar to our X11 canvases,
    but instead of drawing on an X11 window it builds an SVG document.
    """

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.elements = []
        self._init_svg_header()

    def _init_svg_header(self):
        # Create the header and footer for the SVG file.
        self.header = (
            '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
            f'width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n'
        )
        self.footer = '</svg>\n'

    def clear(self):
        """Clear the current SVG elements."""
        self.elements = []

    def create_or_resize(self, width, height):
        """Resize the canvas by updating dimensions and reinitializing the header."""
        self.width = width
        self.height = height
        self._init_svg_header()
        self.clear()

    def flush(self, filename="quickORTEP_export.svg"):
        """
        Finalize the SVG document and write it to a file.
        
        :param filename: Name of the file to save the SVG.
        """
        svg_content = self.header
        for elem in self.elements:
            svg_content += "  " + elem + "\n"
        svg_content += self.footer

        with open(filename, "w") as f:
            f.write(svg_content)
        print(f"SVG export saved to {filename}")

    def _rgb_to_hex(self, color):
        """Convert an (R,G,B) tuple to a hex string (e.g. #ff0000)."""
        return '#{:02x}{:02x}{:02x}'.format(*color)

    def draw_filled_circle(self, cx, cy, radius, color=(255, 0, 0)):
        fill = self._rgb_to_hex(color)
        elem = f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{fill}" stroke="none" />'
        self.elements.append(elem)

    def draw_circle_border(self, cx, cy, radius, color=(0, 0, 0), thickness=2):
        stroke = self._rgb_to_hex(color)
        elem = (f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" '
                f'stroke="{stroke}" stroke-width="{thickness}" />')
        self.elements.append(elem)

    def draw_line(self, x1, y1, x2, y2, thickness=4, color=(0, 0, 0)):
        stroke = self._rgb_to_hex(color)
        elem = (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                f'stroke="{stroke}" stroke-width="{thickness}" />')
        self.elements.append(elem)

    def draw_dashed_line(self, x1, y1, x2, y2, thickness=2, color=(128, 128, 128)):
        stroke = self._rgb_to_hex(color)
        # Define a dash pattern: 8 pixels dash, 4 pixels gap.
        dash_pattern = "8,4"
        elem = (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                f'stroke="{stroke}" stroke-width="{thickness}" '
                f'stroke-dasharray="{dash_pattern}" />')
        self.elements.append(elem)

    def draw_arc(self, cx, cy, rx, ry=None, angle_start_deg=0.0, angle_end_deg=360.0,
                 thickness=2, color=(0, 0, 0)):
        """
        Draw an arc as an SVG path.
        
        The arc is defined by the ellipse centered at (cx,cy) with radii rx and ry.
        The arc starts at angle_start_deg and ends at angle_end_deg.
        """
        if ry is None:
            ry = rx

        stroke = self._rgb_to_hex(color)
        a_start = math.radians(angle_start_deg)
        a_end = math.radians(angle_end_deg)
        x_start = cx + rx * math.cos(a_start)
        y_start = cy + ry * math.sin(a_start)
        x_end = cx + rx * math.cos(a_end)
        y_end = cy + ry * math.sin(a_end)

        # Compute the angle difference (in degrees) to set the large-arc-flag.
        delta_angle = (angle_end_deg - angle_start_deg) % 360
        large_arc_flag = 1 if delta_angle > 180 else 0
        # For simplicity, we always use a sweep_flag of 1.
        sweep_flag = 1

        d = f'M {x_start:.2f},{y_start:.2f} A {rx},{ry} 0 {large_arc_flag} {sweep_flag} {x_end:.2f},{y_end:.2f}'
        elem = (f'<path d="{d}" fill="none" stroke="{stroke}" '
                f'stroke-width="{thickness}" />')
        self.elements.append(elem)

    def draw_filled_arc(self, cx, cy, rx, ry=None, angle_start_deg=0.0, angle_end_deg=360.0,
                          color=(0, 0, 0)):
        """
        Draw a filled arc (wedge) using an SVG path.
        
        The path moves from the center to the start point, draws an arc, and then
        returns to the center.
        """
        if ry is None:
            ry = rx

        fill = self._rgb_to_hex(color)
        a_start = math.radians(angle_start_deg)
        a_end = math.radians(angle_end_deg)
        x_start = cx + rx * math.cos(a_start)
        y_start = cy + ry * math.sin(a_start)
        x_end = cx + rx * math.cos(a_end)
        y_end = cy + ry * math.sin(a_end)

        delta_angle = (angle_end_deg - angle_start_deg) % 360
        large_arc_flag = 1 if delta_angle > 180 else 0
        sweep_flag = 1

        d = (f'M {cx},{cy} L {x_start:.2f},{y_start:.2f} '
             f'A {rx},{ry} 0 {large_arc_flag} {sweep_flag} {x_end:.2f},{y_end:.2f} Z')
        elem = f'<path d="{d}" fill="{fill}" stroke="none" />'
        self.elements.append(elem)


# ----- Test harness -----
if __name__ == "__main__":
    # Create an SVGCanvas and draw some test shapes.
    canvas = SVGCanvas(width=800, height=600)

    # Draw a red filled circle with a black border.
    canvas.draw_filled_circle(400, 300, 50, color=(255, 0, 0))
    canvas.draw_circle_border(400, 300, 50, color=(0, 0, 0), thickness=3)

    # Draw a solid black line.
    canvas.draw_line(100, 100, 700, 100, thickness=4, color=(0, 0, 0))

    # Draw a dashed line.
    canvas.draw_dashed_line(100, 150, 700, 150, thickness=2, color=(0, 0, 255))

    # Draw an arc (non-filled) from 0° to 270°.
    canvas.draw_arc(400, 300, 100, 50, angle_start_deg=0, angle_end_deg=270, thickness=3, color=(0, 128, 0))

    # Draw a filled arc (wedge) from 45° to 135°.
    canvas.draw_filled_arc(400, 300, 80, 80, angle_start_deg=45, angle_end_deg=135, color=(128, 0, 128))

    # Flush the SVG to a file.
    canvas.flush("test_export.svg")
