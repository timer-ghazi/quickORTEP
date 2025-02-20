# x11view/svg_canvas.py

import math

class SVGCanvas:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        # Initialize the SVG document with the header.
        self.svg_data = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">']

    def rgb_to_hex(self, color):
        """Convert an (R, G, B) tuple to a #RRGGBB string."""
        return "#{:02x}{:02x}{:02x}".format(*color)

    def draw_line(self, x1, y1, x2, y2, thickness=1, color=(0,0,0)):
        """
        Draw a line with floating point coordinates.
        The line uses 'round' linecaps and joins to avoid visual gaps.
        """
        col = self.rgb_to_hex(color)
        style = (
            f'stroke="{col}" '
            f'stroke-width="{thickness}" '
            f'stroke-linecap="round" '
            f'stroke-linejoin="round" '
            f'fill="none"'
        )
        line_elem = f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {style} />'
        self.svg_data.append(line_elem)

    def draw_filled_circle(self, cx, cy, radius, color=(0,0,0)):
        """Draw a filled circle at (cx,cy) with the given radius and fill color."""
        col = self.rgb_to_hex(color)
        circle_elem = f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{col}" stroke="none" />'
        self.svg_data.append(circle_elem)

    def draw_circle_border(self, cx, cy, radius, color=(0,0,0), thickness=1):
        """Draw the border of a circle at (cx,cy) with the given radius."""
        col = self.rgb_to_hex(color)
        circle_elem = f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="{col}" stroke-width="{thickness}" />'
        self.svg_data.append(circle_elem)

    def draw_arc(self, cx, cy, rx, ry, angle_start_deg, angle_end_deg, color=(0,0,0), thickness=1):
        """
        Draw an arc (elliptical segment) centered at (cx,cy) with radii rx and ry,
        spanning from angle_start_deg to angle_end_deg.
        """
        col = self.rgb_to_hex(color)
        # Compute start and end points.
        angle_start_rad = math.radians(angle_start_deg)
        angle_end_rad = math.radians(angle_end_deg)
        x_start = cx + rx * math.cos(angle_start_rad)
        y_start = cy + ry * math.sin(angle_start_rad)
        x_end = cx + rx * math.cos(angle_end_rad)
        y_end = cy + ry * math.sin(angle_end_rad)
        # Determine if the arc spans more than 180 degrees.
        large_arc_flag = 1 if ((angle_end_deg - angle_start_deg) % 360) > 180 else 0
        arc_path = f'M {x_start} {y_start} A {rx} {ry} 0 {large_arc_flag} 1 {x_end} {y_end}'
        path_elem = (
            f'<path d="{arc_path}" fill="none" stroke="{col}" stroke-width="{thickness}" '
            f'stroke-linecap="round" stroke-linejoin="round" />'
        )
        self.svg_data.append(path_elem)

    def flush(self, filename=None):
        """Finish the SVG document and either return the SVG string or write it to a file."""
        self.svg_data.append('</svg>')
        svg_content = "\n".join(self.svg_data)
        if filename:
            with open(filename, "w") as f:
                f.write(svg_content)
        else:
            return svg_content

    def clear(self):
        """Clear the canvas and start a new SVG document."""
        self.svg_data = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}">']
