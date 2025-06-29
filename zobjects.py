# zobjects.py

from config import ARC_STYLE, HIGHLIGHT, ATOM_STYLE, VECTOR_STYLE

class ZObject:
    """
    Abstract base for anything that has a 3D depth (z_value) and can draw itself on the 2D canvas.
    """
    def __init__(self, z_value):
        self.z_value = z_value
        self.selected = False  # Selection state

    def draw(self, canvas, color_override=None):
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

    def draw(self, canvas, color_override=None):
        # Use the color override if provided (for fog effect)
        atom_color = color_override if color_override is not None else self.color
        
        # Add shadow for 3D effect (DRAW THIS FIRST)
        if ATOM_STYLE.get("shadow", {}).get("enabled", True):
            # Get shadow parameters from config or use defaults
            shadow_offset_x = ATOM_STYLE.get("shadow", {}).get("offset_x", 2)
            shadow_offset_y = ATOM_STYLE.get("shadow", {}).get("offset_y", 3)
            shadow_radius_factor = ATOM_STYLE.get("shadow", {}).get("radius_factor", 1.1)
            shadow_darkness = ATOM_STYLE.get("shadow", {}).get("darkness", 40)
            shadow_color = (shadow_darkness, shadow_darkness, shadow_darkness)
            
            # Get the current zoom level from VIEWER_INTERACTION if available
            from config import VIEWER_INTERACTION
            # Scale offsets based on zoom ratio (100 is the default zoom level)
            base_zoom = 100.0
            current_zoom = VIEWER_INTERACTION.get("current_zoom", base_zoom)
            zoom_ratio = current_zoom / base_zoom
            
            # Scale the shadow offset with zoom level (minimum 1 pixel)
            scaled_offset_x = max(1, int(shadow_offset_x * zoom_ratio))
            scaled_offset_y = max(1, int(shadow_offset_y * zoom_ratio))
            
            # Calculate shadow position (slightly offset, typically down-right)
            shadow_radius = int(self.radius * shadow_radius_factor)
            shadow_x = self.x2d + scaled_offset_x
            shadow_y = self.y2d + scaled_offset_y
            
            # Draw the shadow
            canvas.draw_filled_circle(shadow_x, shadow_y, shadow_radius, color=shadow_color)
            
        # Draw the filled circle.
        canvas.draw_filled_circle(self.x2d, self.y2d, self.radius, color=atom_color)
        
        # Get the current zoom level for scaled thickness calculations
        from config import VIEWER_INTERACTION
        current_zoom = VIEWER_INTERACTION.get("current_zoom", 100.0)
        
        # Calculate scaled border thickness based on zoom level
        border_thickness = max(
            ATOM_STYLE.get("min_border_thickness_px", 1),
            int(ATOM_STYLE.get("border_thickness", 0.02) * current_zoom)
        )
        
        # Use border color with fog effect if color override is provided
        # For fog effects, don't fully override the border color, but blend it with the theme color
        if color_override is not None:
            # Use a less aggressive fog effect for borders - blend 70% override and 30% original border color 
            r1, g1, b1 = color_override
            r2, g2, b2 = ATOM_STYLE["border_color"]
            border_color = (
                int(0.7 * r1 + 0.3 * r2),
                int(0.7 * g1 + 0.3 * g2),
                int(0.7 * b1 + 0.3 * b2)
            )
        else:
            border_color = ATOM_STYLE["border_color"]
        
        # Draw the atom border with scaled thickness
        canvas.draw_circle_border(
            self.x2d, self.y2d, 
            self.radius,
            color=border_color, 
            thickness=border_thickness
        )
        
        # Calculate scaled meridian arc thickness
        meridian_thickness = max(
            ARC_STYLE.get("min_meridian_thickness_px", 1),
            int(ARC_STYLE.get("meridian_thickness", 0.03) * current_zoom)
        )
        
        # Get arc flattening factor
        flatten = ARC_STYLE["flatten"]
        
        # Use the same border_color that has been adjusted for fog
        canvas.draw_arc(
            self.x2d, self.y2d, 
            self.radius, int(self.radius * flatten),
            angle_start_deg=180, angle_end_deg=360,
            color=border_color, 
            thickness=meridian_thickness
        )
        
        canvas.draw_arc(
            self.x2d, self.y2d, 
            int(self.radius * flatten), self.radius,
            angle_start_deg=270, angle_end_deg=450,
            color=border_color, 
            thickness=meridian_thickness
        )
        
        # Add 3D highlight effect
        if ATOM_STYLE.get("highlight", {}).get("enabled", True):
            # Get highlight parameters from config or use defaults
            size_ratio = ATOM_STYLE.get("highlight", {}).get("size_ratio", 0.3)
            offset_ratio = ATOM_STYLE.get("highlight", {}).get("offset_ratio", 0.4)
            brightness_factor = ATOM_STYLE.get("highlight", {}).get("brightness_factor", 1.5)
            use_white_for_dark = ATOM_STYLE.get("highlight", {}).get("use_white_for_dark", True)
            
            # Calculate highlight position (upper-left quadrant by default)
            highlight_radius = max(2, int(self.radius * size_ratio))
            highlight_offset_x = -int(self.radius * offset_ratio)
            highlight_offset_y = -int(self.radius * offset_ratio)
            highlight_x = self.x2d + highlight_offset_x
            highlight_y = self.y2d + highlight_offset_y
            
            # Calculate highlight color based on atom color (use the overridden color if provided)
            base_r, base_g, base_b = atom_color
            # Use perceptual luminance instead of simple RGB average (ITU-R BT.709 coefficients)
            luminance = 0.299 * base_r + 0.587 * base_g + 0.114 * base_b
            if use_white_for_dark and luminance < 70:
                # For dark atoms, use white highlight
                highlight_color = (255, 255, 255)
            else:
                # For lighter atoms, create a light tinted highlight by blending with white
                # This works better for saturated colors like pure red oxygen
                white_blend_factor = 0.6  # Blend 60% white with 40% original color
                highlight_r = min(255, int(base_r * (1 - white_blend_factor) + 255 * white_blend_factor))
                highlight_g = min(255, int(base_g * (1 - white_blend_factor) + 255 * white_blend_factor))
                highlight_b = min(255, int(base_b * (1 - white_blend_factor) + 255 * white_blend_factor))
                highlight_color = (highlight_r, highlight_g, highlight_b)
            
            # Draw the highlight
            canvas.draw_filled_circle(highlight_x, highlight_y, highlight_radius, color=highlight_color)

        # Draw atom label if enabled
        if ATOM_STYLE.get("label", {}).get("enabled", False) and self.atom is not None:
            # Get label settings from config or use defaults
            font_size = ATOM_STYLE.get("label", {}).get("font_size", 12)
            offset_x = ATOM_STYLE.get("label", {}).get("offset_x", -8)
            offset_y = ATOM_STYLE.get("label", {}).get("offset_y", -8)
            background = ATOM_STYLE.get("label", {}).get("background", True)
            background_color = ATOM_STYLE.get("label", {}).get("background_color", (0, 0, 0, 180))
            show_symbols = ATOM_STYLE.get("label", {}).get("show_symbols", False)
            padding = ATOM_STYLE.get("label", {}).get("padding", 2)
            
            # Get current zoom level for scaling offsets
            zoom_ratio = current_zoom / 100.0
            scaled_offset_x = int(offset_x * zoom_ratio)
            scaled_offset_y = int(offset_y * zoom_ratio)
            
            # Calculate label position (top-left quadrant by default)
            label_x = self.x2d + scaled_offset_x
            label_y = self.y2d + scaled_offset_y
            
            # Create label text
            if show_symbols and hasattr(self.atom, 'symbol'):
                label_text = f"{self.atom.symbol}{self.atom.index}"
            else:
                label_text = f"{self.atom.index}"
            
            # Get label color - use provided color or atom color adjusted for visibility
            label_color = ATOM_STYLE.get("label", {}).get("color", None)
            if label_color is None:
                from config import CURRENT_THEME
                label_color = CURRENT_THEME.get("default_text_color", (255, 255, 255))
            
            # Apply fog to label color if color_override is provided
            if color_override is not None:
                # Blend original label color (70%) with fog color (30%)
                r1, g1, b1 = label_color
                r2, g2, b2 = color_override
                label_color = (
                    int(0.7 * r1 + 0.3 * r2),
                    int(0.7 * g1 + 0.3 * g2),
                    int(0.7 * b1 + 0.3 * b2)
                )
            
            # Draw label background if enabled
            if background:
                # Get text dimensions to calculate background size
                text_width, text_height = canvas.get_text_dimensions(label_text, font_size)
                bg_x = label_x - padding
                bg_y = label_y - text_height - padding
                bg_width = text_width + 2 * padding
                bg_height = text_height + 2 * padding
                
                # Handle RGBA colors (for SVG) vs RGB colors (for X11)
                # Convert (r, g, b, a) to (r, g, b) if needed
                if len(background_color) == 4:
                    r, g, b, a = background_color
                    # For X11, just use the RGB components
                    bg_color = (r, g, b)
                else:
                    bg_color = background_color
                
                # Draw background rectangle
                canvas.draw_filled_rect(bg_x, bg_y, bg_width, bg_height, color=bg_color)
            
            # Draw label text
            canvas.draw_text(label_x, label_y, label_text, color=label_color, font_size=font_size)

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

    def draw(self, canvas, color_override=None):
        # Use the color override if provided (for fog effect)
        segment_color = color_override if color_override is not None else self.color
        
        # Draw the bond segment.
        canvas.draw_line(self.x1, self.y1, self.x2, self.y2,
                         thickness=self.thickness, color=segment_color)
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

    def draw(self, canvas, color_override=None):
        """
        Draw the arrowhead on the canvas.
        
        Parameters:
            canvas: The canvas to draw on
            color_override: Optional color override for fog effect
        """
        # Use the color override if provided (for fog effect)
        arrow_color = color_override if color_override is not None else self.color
        
        # Use the canvas drawing methods
        if self.filled:
            # Draw as a filled triangle
            # Convert points to integers to avoid potential type issues
            x1, y1 = int(self.tip_x), int(self.tip_y)
            x2, y2 = int(self.corner1_x), int(self.corner1_y)
            x3, y3 = int(self.corner2_x), int(self.corner2_y)
            
            canvas.draw_filled_triangle(
                x1, y1,
                x2, y2,
                x3, y3,
                color=arrow_color
            )
        else:
            # Draw as an outline triangle
            canvas.draw_triangle(
                self.tip_x, self.tip_y,
                self.corner1_x, self.corner1_y,
                self.corner2_x, self.corner2_y,
                color=arrow_color,
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
