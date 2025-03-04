# zobjects.py

from config import ARC_STYLE, HIGHLIGHT, ATOM_STYLE, VECTOR_STYLE

class ZObject:
    """
    Abstract base for anything that has a 3D depth (z_value) and can draw itself on the 2D canvas.
    """
    def __init__(self, z_value):
        self.z_value = z_value
        self.selected = False  # Selection state

    def draw(self, canvas):
        raise NotImplementedError("Subclasses must implement draw().")


class ZAtom(ZObject):
    """
    Represents one atom's 2D circle (and arcs) along with its 3D depth.
    
    An optional reference to the underlying atom data (e.g. ORTEP_Atom) is stored in self.atom.
    """
    def __init__(self, x2d, y2d, radius, color, z_value):
        super().__init__(z_value)
        self.x2d = x2d
        self.y2d = y2d
        self.radius = radius
        self.color = color
        self.atom = None  # Underlying ORTEP_Atom instance

    def draw(self, canvas):
        # Draw the filled circle.
        canvas.draw_filled_circle(self.x2d, self.y2d, self.radius, color=self.color)
        # Draw the atom border using the configured border color.
        canvas.draw_circle_border(self.x2d, self.y2d, self.radius, color=ATOM_STYLE["border_color"], thickness=2)
        
        # Draw ORTEP arcs.
        meridian_thickness = ARC_STYLE["meridian_thickness"]
        flatten = ARC_STYLE["flatten"]
        canvas.draw_arc(
            self.x2d, self.y2d, 
            self.radius, int(self.radius * flatten),
            angle_start_deg=180, angle_end_deg=360,
            color=ATOM_STYLE["border_color"], thickness=meridian_thickness
        )
        canvas.draw_arc(
            self.x2d, self.y2d, 
            int(self.radius * flatten), self.radius,
            angle_start_deg=270, angle_end_deg=450,
            color=ATOM_STYLE["border_color"], thickness=meridian_thickness
        )

        # Draw highlight if the atom is selected.
        if self.selected:
            # The extra border offset (+4) remains hard-coded; adjust if needed later.
            canvas.draw_circle_border(self.x2d, self.y2d, self.radius + 4, color=HIGHLIGHT["color"], thickness=HIGHLIGHT["thickness_delta"])

    def contains(self, x, y):
        """
        Determine whether the point (x, y) lies within the atom's circle.
        
        Uses the standard circle equation:
        
        $$
        (x - x_{center})^2 + (y - y_{center})^2 \\leq \\text{radius}^2
        $$
        """
        return (x - self.x2d)**2 + (y - self.y2d)**2 <= self.radius**2


class ZSegment(ZObject):
    """
    Represents a short bond segment in 2D with a single z_value (its 3D midpoint).
    
    An optional reference to the underlying bond data is stored in self.bond.
    """
    def __init__(self, x1, y1, x2, y2, z_value, thickness=16, color=(0, 0, 0)):
        super().__init__(z_value)
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.thickness = thickness
        self.color = color
        self.bond = None  # Underlying bond instance
        self.vector = None  # Underlying vector instance (for vector segments)

    def draw(self, canvas):
        # Draw the bond segment.
        canvas.draw_line(self.x1, self.y1, self.x2, self.y2,
                         thickness=self.thickness, color=self.color)
        # If selected, overlay a highlight using the configured settings.
        if self.selected:
            canvas.draw_line(self.x1, self.y1, self.x2, self.y2,
                             thickness=self.thickness + HIGHLIGHT["thickness_delta"], color=HIGHLIGHT["color"])

    def contains(self, x, y, tolerance=3):
        """
        Determine whether the point (x, y) is close to the bond segment.
        
        The method projects (x, y) onto the line defined by the segment endpoints
        and checks if the distance from the point to the projection is within a threshold.
        
        The projection factor is computed as:
        
        $$
        t = \\frac{(x-x_1)(x_2-x_1) + (y-y_1)(y_2-y_1)}{(x_2-x_1)^2 + (y_2-y_1)^2}
        $$
        
        where $t$ is clamped to the interval $[0,1]$. The point is considered within
        tolerance if:
        
        $$
        (x - c_x)^2 + (y - c_y)^2 \\leq \\Bigl(tolerance + \\frac{thickness}{2}\\Bigr)^2
        $$
        """
        dx = self.x2 - self.x1
        dy = self.y2 - self.y1
        if dx == 0 and dy == 0:
            return False
        t = ((x - self.x1) * dx + (y - self.y1) * dy) / (dx * dx + dy * dy)
        t = max(0, min(1, t))
        cx = self.x1 + t * dx
        cy = self.y1 + t * dy
        dist_sq = (x - cx)**2 + (y - cy)**2
        effective_tolerance = tolerance + self.thickness / 2
        return dist_sq <= effective_tolerance**2


class ZArrowHead(ZObject):
    """
    Represents a triangular arrowhead in 2D with a single z_value.
    
    The arrowhead is defined by three points: the tip and two base corners.
    An optional reference to the underlying vector data is stored in self.vector.
    """
    def __init__(self, tip_x, tip_y, corner1_x, corner1_y, corner2_x, corner2_y, 
                 z_value, color=(0, 0, 0), filled=True):
        """
        Initialize an arrowhead with the given tip and corner points.
        
        Parameters:
            tip_x, tip_y: The coordinates of the arrowhead tip
            corner1_x, corner1_y: The coordinates of the first base corner
            corner2_x, corner2_y: The coordinates of the second base corner
            z_value: The depth for z-ordering
            color: The color of the arrowhead
            filled: Whether to draw a filled triangle (True) or just an outline (False)
        """
        super().__init__(z_value)
        self.tip_x = tip_x
        self.tip_y = tip_y
        self.corner1_x = corner1_x
        self.corner1_y = corner1_y
        self.corner2_x = corner2_x
        self.corner2_y = corner2_y
        self.color = color
        self.filled = filled
        self.thickness = max(VECTOR_STYLE["min_thickness_px"], 
                            int(VECTOR_STYLE["thickness"] * 1.5))  # Slightly thicker than shaft
        self.vector = None  # Underlying vector instance

    def draw(self, canvas):
        """
        Draw the arrowhead on the canvas.
        """
        # Use the new triangle drawing methods
        if self.filled:
            # Draw as a filled triangle
            canvas.draw_filled_triangle(
                self.tip_x, self.tip_y,
                self.corner1_x, self.corner1_y,
                self.corner2_x, self.corner2_y,
                color=self.color
            )
        else:
            # Draw as an outline triangle
            canvas.draw_triangle(
                self.tip_x, self.tip_y,
                self.corner1_x, self.corner1_y,
                self.corner2_x, self.corner2_y,
                color=self.color,
                thickness=self.thickness
            )
        
        # Draw highlight if selected
        if self.selected:
            highlight_thickness = self.thickness + HIGHLIGHT["thickness_delta"]
            canvas.draw_triangle(
                self.tip_x, self.tip_y,
                self.corner1_x, self.corner1_y,
                self.corner2_x, self.corner2_y,
                color=HIGHLIGHT["color"],
                thickness=highlight_thickness
            )

    def contains(self, x, y, tolerance=3):
        """
        Determine whether the point (x, y) is inside the arrowhead triangle.
        
        This uses a barycentric coordinate calculation to check if the point
        is inside the triangle. A tolerance value is added to make selection easier.
        
        Parameters:
            x, y: The point to test
            tolerance: Additional radius for hit detection
            
        Returns:
            bool: True if point is within the triangle (plus tolerance)
        """
        # First check if the point is near any of the three edges
        edges = [
            (self.tip_x, self.tip_y, self.corner1_x, self.corner1_y),
            (self.tip_x, self.tip_y, self.corner2_x, self.corner2_y),
            (self.corner1_x, self.corner1_y, self.corner2_x, self.corner2_y)
        ]
        
        for x1, y1, x2, y2 in edges:
            dx = x2 - x1
            dy = y2 - y1
            if dx == 0 and dy == 0:
                continue
                
            # Project onto the line
            t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)
            t = max(0, min(1, t))
            
            # Calculate closest point on the line
            cx = x1 + t * dx
            cy = y1 + t * dy
            
            # Check distance
            dist_sq = (x - cx)**2 + (y - cy)**2
            if dist_sq <= tolerance**2:
                return True
        
        # If not near any edge, check if inside the triangle using barycentric coordinates
        # Area of the triangle
        area_full = 0.5 * abs(
            (self.corner1_x - self.tip_x) * (self.corner2_y - self.tip_y) -
            (self.corner2_x - self.tip_x) * (self.corner1_y - self.tip_y)
        )
        
        if area_full < 1e-6:  # Degenerate triangle
            return False
            
        # Barycentric coordinates
        alpha = 0.5 * abs(
            (self.corner1_x - x) * (self.corner2_y - y) -
            (self.corner2_x - x) * (self.corner1_y - y)
        ) / area_full
        
        beta = 0.5 * abs(
            (self.tip_x - x) * (self.corner2_y - y) -
            (self.corner2_x - x) * (self.tip_y - y)
        ) / area_full
        
        gamma = 1 - alpha - beta
        
        # If all barycentric coordinates are between 0 and 1, the point is inside
        return 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1
