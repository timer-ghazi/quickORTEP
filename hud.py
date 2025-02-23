# hud.py

class HUDPanel:
    """
    HUDPanel manages the Heads-Up Display (HUD) overlay for the ORTEP viewer.

    Responsibilities:
      - Store and update a list of HUD text lines.
      - Render these lines on a canvas with configurable position, spacing, font size, and color.

    Attributes:
      - x (int): X-coordinate for the HUD text starting point.
      - y_offset (int): Vertical offset from the bottom of the canvas.
      - line_spacing (int): Vertical spacing between consecutive HUD lines.
      - font_size (int): Font size for the HUD text.
      - color (tuple): Color (R, G, B) for the text.
      - lines (list of str): The current lines to display in the HUD.
    """
    def __init__(self, x=10, y_offset=20, line_spacing=20, font_size=14, color=(0, 0, 0)):
        self.x = x
        self.y_offset = y_offset
        self.line_spacing = line_spacing
        self.font_size = font_size
        self.color = color
        self.lines = []  # Initialize with an empty list of HUD lines.

    def update_lines(self, lines):
        """
        Update the HUD content.

        Parameters:
          - lines (list of str): New list of strings to display in the HUD.
        """
        self.lines = lines

    def draw(self, canvas):
        """
        Draw the HUD on the provided canvas. The HUD is rendered starting at a fixed offset
        from the bottom left corner, with each line spaced by 'line_spacing'.

        Parameters:
          - canvas: A canvas object that must support a 'draw_text' method with parameters:
                    x, y, text, color, and font_size.
        """
        for i, line in enumerate(self.lines):
            # Calculate the y-coordinate for each line.
            y = canvas.height - self.y_offset - i * self.line_spacing
            canvas.draw_text(self.x, y, line, color=self.color, font_size=self.font_size)
