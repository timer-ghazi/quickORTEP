# x11view/svg_canvas.py

import math
from .DejaVuSansMono import EMBEDDED_FONT_WOFF2  

class SVGCanvas:
    def __init__(self, width, height, background_color=(255, 255, 255)):
        """
        Initialize the SVG canvas with the given dimensions and background color.
        
        Parameters:
            width: Canvas width in pixels
            height: Canvas height in pixels
            background_color: RGB tuple (r, g, b) for the background (default: white)
        """
        self.width = width
        self.height = height
        self.background_color = background_color
        
        # Initialize the SVG document with the header
        self.svg_data = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">']
        self.embed_font_once()
        
        # Add background rectangle if not white
        if self.background_color != (255, 255, 255):
            self.add_background_rect()

    def rgb_to_hex(self, color):
        """Convert an (R, G, B) tuple to a #RRGGBB string."""
        return "#{:02x}{:02x}{:02x}".format(*color)

    def add_background_rect(self):
        """Add a background rectangle to the SVG with the current background color."""
        bg_color_hex = self.rgb_to_hex(self.background_color)
        bg_rect = f'<rect x="0" y="0" width="{self.width}" height="{self.height}" fill="{bg_color_hex}" />'
        self.svg_data.append(bg_rect)
    
    def embed_font_once(self):
        """
        Inserts a <style> block in the SVG that defines a @font-face
        using our inlined WOFF2 data. This ensures any <text> that uses
        font-family="DejaVu Sans Mono" will render identically everywhere.
        """
        style_block = f"""
    <style>
    @font-face {{
        font-family: "DejaVu Sans Mono";
        src: url("data:font/woff2;base64,{EMBEDDED_FONT_WOFF2}") format("woff2");
    }}
    </style>
    """
        self.svg_data.append(style_block)

    def set_background_color(self, color):
        """
        Set the background color for the SVG canvas.
        
        Parameters:
            color: RGB tuple (r, g, b) with values in range 0-255
        """
        self.background_color = color
        
        # Clear and recreate the SVG with the new background
        self.clear()

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

    # -------------------------------------------------------
    # New: The draw_text() method using the embedded WOFF2
    # -------------------------------------------------------
    def draw_text(self,
                  x,
                  y,
                  text,
                  color=(0,0,0),
                  font_size=12,
                  font_candidates=None):
        """
        Draw text so that (x, y) is the text baseline (if viewer supports dominant-baseline).
        color is (R,G,B), font_size in px, and we ignore font_candidates
        since we have an embedded WOFF2 font named "DejaVu Sans Mono".
        """
        col_hex = self.rgb_to_hex(color)
        font_family = "DejaVu Sans Mono"

        # By specifying dominant-baseline="alphabetic", we attempt to treat (x, y) as baseline.
        # Not all renderers do precisely the same baseline offset, so test if you need an extra shift.
        text_elem = (
            f'<text x="{x}" y="{y}" '
            f'fill="{col_hex}" '
            f'font-size="{font_size}px" '
            f'font-family="{font_family}" '
            f'dominant-baseline="alphabetic">'
            f'{text}</text>'
        )
        self.svg_data.append(text_elem)

        
    def draw_rect(self,
                  x,
                  y,
                  width,
                  height,
                  color=(0,0,0),
                  thickness=2):
        """
        Draw a rectangle border (unfilled) in the SVG.
        (x, y) is the top-left corner, 'thickness' is the stroke width in pixels.
        """
        col_hex = self.rgb_to_hex(color)
        rect_elem = (
            f'<rect x="{x}" y="{y}" width="{width}" height="{height}" '
            f'stroke="{col_hex}" stroke-width="{thickness}" fill="none" '
            f'stroke-linecap="round" stroke-linejoin="round" />'
        )
        self.svg_data.append(rect_elem)
    
    def draw_filled_rect(self,
                         x,
                         y,
                         width,
                         height,
                         color=(0,0,0)):
        """
        Draw a filled rectangle at (x, y), covering the given width/height.
        """
        col_hex = self.rgb_to_hex(color)
        rect_elem = (
            f'<rect x="{x}" y="{y}" width="{width}" height="{height}" '
            f'fill="{col_hex}" stroke="none" />'
        )
        self.svg_data.append(rect_elem)


    def draw_filled_arc(self,
                        cx: float,
                        cy: float,
                        rx: float,
                        ry: float = None,
                        angle_start_deg: float = 0.0,
                        angle_end_deg: float = 360.0,
                        color=(0, 0, 0)):
        """
        Draw a filled wedge or elliptical sector. (cx, cy) is the center;
        rx, ry are the ellipse radii; angle_start_deg and angle_end_deg are in degrees.
        We assume angles go counterclockwise, matching typical geometry.
        """
    
        import math
        if ry is None:
            ry = rx
    
        col = self.rgb_to_hex(color)
    
        # Compute the extent in degrees, ensuring it's positive for a CCW arc
        extent_deg = angle_end_deg - angle_start_deg
        if extent_deg <= 0:
            extent_deg += 360.0
    
        # Convert angles to radians
        start_rad = math.radians(angle_start_deg)
        end_rad   = math.radians(angle_start_deg + extent_deg)
    
        # Calculate start and end points of the arc
        x_start = cx + rx * math.cos(start_rad)
        y_start = cy + ry * math.sin(start_rad)
        x_end = cx + rx * math.cos(end_rad)
        y_end = cy + ry * math.sin(end_rad)
    
        # Determine if the arc is > 180 degrees for the large_arc_flag
        large_arc_flag = 1 if extent_deg > 180 else 0
    
        # Build the path data:
        # - Move to the arc start point
        # - Draw the elliptical arc
        # - Line back to center
        # - Close the path
        path_data = (
            f"M {x_start} {y_start} "
            f"A {rx} {ry} 0 {large_arc_flag} 1 {x_end} {y_end} "
            f"L {cx} {cy} Z"
        )
    
        path_elem = f'<path d="{path_data}" fill="{col}" stroke="none" />'
        self.svg_data.append(path_elem)


    def draw_dashed_line(self,
                         x1, y1,
                         x2, y2,
                         thickness=2,
                         color=(128, 128, 128)):
        """
        Draw a dashed line from (x1, y1) to (x2, y2) with the given thickness and color.
        We mimic the '8,4' dash pattern from the X11 code.
        """
        col = self.rgb_to_hex(color)
        dash_pattern = "8,4"  # 8px dash, 4px gap, just like our X11 code
    
        style = (
            f'stroke="{col}" '
            f'stroke-width="{thickness}" '
            f'stroke-dasharray="{dash_pattern}" '
            f'stroke-linecap="round" '
            f'stroke-linejoin="round" '
            f'fill="none"'
        )
        line_elem = (
            f'<line x1="{x1}" y1="{y1}" '
            f'x2="{x2}" y2="{y2}" '
            f'{style} />'
        )
        self.svg_data.append(line_elem)

    def draw_triangle(self,
                     x1, y1,
                     x2, y2,
                     x3, y3,
                     color=(0, 0, 0),
                     thickness=2):
        """
        Draw an unfilled triangle in SVG.
        """
        col = self.rgb_to_hex(color)
        points = f"{x1},{y1} {x2},{y2} {x3},{y3}"
        
        polygon_elem = (
            f'<polygon points="{points}" '
            f'fill="none" '
            f'stroke="{col}" '
            f'stroke-width="{thickness}" '
            f'stroke-linejoin="round" />'
        )
        self.svg_data.append(polygon_elem)

    def draw_filled_triangle(self,
                            x1, y1,
                            x2, y2,
                            x3, y3,
                            color=(0, 0, 0)):
        """
        Draw a filled triangle in SVG.
        """
        col = self.rgb_to_hex(color)
        points = f"{x1},{y1} {x2},{y2} {x3},{y3}"
        
        polygon_elem = f'<polygon points="{points}" fill="{col}" stroke="none" />'
        self.svg_data.append(polygon_elem)
    
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
        """Clear the canvas and start a new SVG document with the current background color."""
        self.svg_data = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}">']
        self.embed_font_once()
        
        # Add background rectangle if not white
        if self.background_color != (255, 255, 255):
            self.add_background_rect()