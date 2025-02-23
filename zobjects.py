# zobjects.py

from config import ARC_STYLE, HIGHLIGHT, ATOM_STYLE

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
