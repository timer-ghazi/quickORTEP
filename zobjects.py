# zobjects.py

class ZObject:
    """
    Abstract base for anything that has a 3D depth (z_value) and can draw itself on the 2D canvas.
    """
    def __init__(self, z_value):
        self.z_value = z_value

    def draw(self, canvas):
        raise NotImplementedError("Subclasses must implement draw().")


class ZAtom(ZObject):
    """
    Represents one atom's 2D circle (and arcs) along with its 3D depth.
    """
    def __init__(self, x2d, y2d, radius, color, z_value):
        super().__init__(z_value)
        self.x2d = x2d
        self.y2d = y2d
        self.radius = radius
        self.color = color
        # ... any other style fields ...

    def draw(self, canvas):
        # same code you had in _draw_atoms for drawing a single atom
        canvas.draw_filled_circle(self.x2d, self.y2d, self.radius, color=self.color)
        canvas.draw_circle_border(self.x2d, self.y2d, self.radius, color=(0,0,0), thickness=2)
        # plus your ORTEP arcs, etc.

        # ORTEP arcs
        MERIDIAN_THICK = 2
        H_FLATTEN = 0.4

        # horizontal ellipse (front half: 180..360)
        canvas.draw_arc(
            self.x2d, self.y2d, 
            self.radius, int(self.radius * H_FLATTEN),
            angle_start_deg=180, angle_end_deg=360,
            color=(0,0,0), thickness=MERIDIAN_THICK
        )
        # vertical ellipse (front half: 270..450)
        canvas.draw_arc(
            self.x2d, self.y2d, 
            int(self.radius * H_FLATTEN), self.radius,
            angle_start_deg=270, angle_end_deg=450,
            color=(0,0,0), thickness=MERIDIAN_THICK
        )


class ZSegment(ZObject):
    """
    Represents a short bond segment in 2D with a single z_value (its 3D midpoint).
    """
    def __init__(self, x1, y1, x2, y2, z_value, thickness=16, color=(0,0,0)):
        super().__init__(z_value)
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.thickness = thickness
        self.color = color

    def draw(self, canvas):
        canvas.draw_line(self.x1, self.y1, self.x2, self.y2,
                         thickness=self.thickness, color=self.color)
