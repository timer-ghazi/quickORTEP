# vectors.py

import math
import numpy as np
from config import VECTOR_STYLE, AXES_STYLE
from geometry_utils import project_point, project_points
from zobjects import ZSegment, ZArrowHead

class Vector:
    """
    Represents a 3D vector with a starting point and direction.
    
    A vector has:
    - start_point: (x, y, z) coordinates of the vector's start
    - direction: (dx, dy, dz) normalized direction vector
    - length: magnitude of the vector
    - color: (r, g, b) tuple for the vector's color
    """
    def __init__(self, start_point, end_point=None, direction=None, length=None, color=None):
        """
        Initialize a vector with either:
        1. start_point, end_point
        2. start_point, direction, length
        
        Parameters:
            start_point: (x, y, z) coordinates of the vector's start
            end_point: (x, y, z) coordinates of the vector's end (optional)
            direction: (dx, dy, dz) normalized direction vector (optional)
            length: magnitude of the vector (optional)
            color: (r, g, b) tuple for the vector's color (None = use default)
        """
        self.start_point = np.array(start_point)
        self.selected = False  # For selection highlighting
        
        # Set the color (use default if not specified)
        self.color = color if color is not None else VECTOR_STYLE["color"]
        
        # Calculate direction and length from end point if provided
        if end_point is not None:
            end_point = np.array(end_point)
            displacement = end_point - self.start_point
            self.length = np.linalg.norm(displacement)
            
            # Avoid division by zero for zero-length vectors
            if self.length > 1e-6:
                self.direction = displacement / self.length
            else:
                self.direction = np.array([0.0, 0.0, 0.0])
                self.length = 0.0
        # Or use provided direction and length
        elif direction is not None and length is not None:
            self.direction = np.array(direction)
            # Normalize direction if it's not already normalized
            norm = np.linalg.norm(self.direction)
            if norm > 1e-6:
                self.direction = self.direction / norm
            self.length = length
        else:
            raise ValueError("Must provide either end_point or both direction and length")
    
    def get_end_point(self):
        """
        Calculate the endpoint of the vector based on start point, direction, and length.
        
        Returns:
            np.array: (x, y, z) coordinates of the vector's end
        """
        return self.start_point + self.direction * self.length
    
    def get_segments(self, rotated_coords, view_params):
        """
        Generate the segments representing the vector's shaft and arrowhead.
        
        Parameters:
            rotated_coords: tuple containing ((x1, y1, z1, r1), (x2, y2, z2, r2))
                           for start and end points after rotation
            view_params: view parameters for projection
            
        Returns:
            list: ZSegment objects for the shaft and arrowhead
        """
        # Extract rotated coordinates
        (x1, y1, z1, _), (x2, y2, z2, _) = rotated_coords
        
        # Create vectors in rotated space
        p1 = np.array([x1, y1, z1])
        p2 = np.array([x2, y2, z2])
        
        # Calculate vector properties in rotated space
        v = p2 - p1
        dist = np.linalg.norm(v)
        
        if dist < 1e-6:
            return []  # Can't draw zero-length vector
        
        # Leave space for arrowhead
        arrowhead_length = VECTOR_STYLE["arrowhead"]["length"]
        shaft_end_factor = 1.0 - (arrowhead_length / dist) if dist > arrowhead_length else 0.5
        
        # Create shaft segments
        segments = self._get_shaft_segments(p1, p2, shaft_end_factor, view_params)
        
        # Create arrowhead segments
        arrowhead_segments = self._get_arrowhead_segments(p1, p2, shaft_end_factor, view_params)
        
        # Combine all segments
        segments.extend(arrowhead_segments)
        
        return segments
    
    def _get_shaft_segments(self, p1, p2, shaft_end_factor, view_params):
        """
        Generate segments for the vector's shaft.
        
        Parameters:
            p1: np.array with (x, y, z) of start point in rotated space
            p2: np.array with (x, y, z) of end point in rotated space
            shaft_end_factor: fraction of the total length where shaft ends
            view_params: view parameters for projection
            
        Returns:
            list: ZSegment objects for the shaft
        """
        # Calculate shaft end point
        v = p2 - p1
        shaft_end = p1 + v * shaft_end_factor
        
        # Create segments
        segment_length = VECTOR_STYLE["segment_length"]
        shaft_dist = np.linalg.norm(shaft_end - p1)
        
        # If shaft is too short, just return a single segment
        if shaft_dist <= segment_length:
            # Project points to screen space
            start_2d = project_point(p1[0], p1[1], p1[2], view_params)
            end_2d = project_point(shaft_end[0], shaft_end[1], shaft_end[2], view_params)
            
            # Create segment
            z_mid = (p1[2] + shaft_end[2]) / 2
            segment = ZSegment(
                x1=start_2d[0],
                y1=start_2d[1],
                x2=end_2d[0],
                y2=end_2d[1],
                z_value=z_mid,
                thickness=max(VECTOR_STYLE["min_thickness_px"], 
                             int(VECTOR_STYLE["thickness"] * view_params.scale)),
                color=self.color
            )
            return [segment]
        
        # Multiple segments needed
        N = max(1, int(math.ceil(shaft_dist / segment_length)))
        dt = 1.0 / N
        
        # Generate points along the shaft
        t_values = np.linspace(0, 1.0, N + 1)
        p1_array = np.array([p1]*len(t_values[:-1]))
        dir_array = np.array([shaft_end - p1]*len(t_values[:-1]))
        
        segment_starts = p1_array + dir_array * t_values[:-1, np.newaxis]
        segment_ends = p1_array + dir_array * t_values[1:, np.newaxis]
        
        # Project all points to screen space
        starts_2d = project_points(segment_starts, view_params)
        ends_2d = project_points(segment_ends, view_params)
        
        # Create segments
        shaft_thickness = max(VECTOR_STYLE["min_thickness_px"], 
                             int(VECTOR_STYLE["thickness"] * view_params.scale))
        
        segments = []
        for i in range(N):
            z_mid = (segment_starts[i, 2] + segment_ends[i, 2]) / 2
            segment = ZSegment(
                x1=int(starts_2d[i, 0]),
                y1=int(starts_2d[i, 1]),
                x2=int(ends_2d[i, 0]),
                y2=int(ends_2d[i, 1]),
                z_value=z_mid,
                thickness=shaft_thickness,
                color=self.color
            )
            segments.append(segment)
        
        return segments
    
    def _get_arrowhead_segments(self, p1, p2, shaft_end_factor, view_params):
        """
        Generate the arrowhead for the vector.
        
        Parameters:
            p1: np.array with (x, y, z) of start point in rotated space
            p2: np.array with (x, y, z) of end point in rotated space
            shaft_end_factor: fraction of the total length where shaft ends
            view_params: view parameters for projection
            
        Returns:
            list: Objects for the arrowhead (ZArrowHead or ZSegment)
        """
        # Calculate arrow properties
        arrowhead_length = VECTOR_STYLE["arrowhead"]["length"]
        arrowhead_width = arrowhead_length * VECTOR_STYLE["arrowhead"]["width_ratio"]
        
        # Direction of the vector in rotated space
        v = p2 - p1
        v_length = np.linalg.norm(v)
        direction = v / v_length if v_length > 1e-6 else np.array([0, 0, 0])
        
        # Calculate shaft endpoint (where arrowhead begins)
        shaft_end = p1 + direction * (v_length * shaft_end_factor)
        
        # Find perpendicular vectors for arrowhead width
        if abs(direction[2]) < 0.9:  # If not pointing mostly along z
            perp = np.array([direction[1], -direction[0], 0])
        else:  # If pointing mostly along z
            perp = np.array([1, 0, -direction[0]/direction[2] if abs(direction[2]) > 1e-6 else 0])
        
        # Normalize perpendicular vector
        perp_length = np.linalg.norm(perp)
        if perp_length > 1e-6:
            perp = perp / perp_length
        
        # Calculate arrowhead points
        arrowhead_width_half = arrowhead_width / 2
        corner1 = shaft_end + perp * arrowhead_width_half
        corner2 = shaft_end - perp * arrowhead_width_half
        
        # Project to screen space
        tip_2d = project_point(p2[0], p2[1], p2[2], view_params)
        corner1_2d = project_point(corner1[0], corner1[1], corner1[2], view_params)
        corner2_2d = project_point(corner2[0], corner2[1], corner2[2], view_params)
        
        # Determine thickness
        thickness = max(VECTOR_STYLE["min_thickness_px"], 
                       int(VECTOR_STYLE["thickness"] * view_params.scale))
        
        # Get arrowhead color (use the same color as shaft if not specified)
        arrowhead_color = VECTOR_STYLE["arrowhead"]["color"] if VECTOR_STYLE["arrowhead"]["color"] else self.color
        
        segments = []
        
        # Use ZArrowHead if we want a filled arrowhead
        if VECTOR_STYLE["arrowhead"]["filled"]:
            arrow_head = ZArrowHead(
                tip_x=tip_2d[0],
                tip_y=tip_2d[1],
                corner1_x=corner1_2d[0],
                corner1_y=corner1_2d[1],
                corner2_x=corner2_2d[0],
                corner2_y=corner2_2d[1],
                z_value=(p2[2] + shaft_end[2]) / 2,  # Average z-value for z-ordering
                color=arrowhead_color,
                filled=True
            )
            # Add a reference to the source vector
            arrow_head.vector = self
            segments.append(arrow_head)
        else:
            # Create line segments for an unfilled arrowhead
            # First line from base corner 1 to tip
            segments.append(ZSegment(
                x1=corner1_2d[0],
                y1=corner1_2d[1],
                x2=tip_2d[0],
                y2=tip_2d[1],
                z_value=(corner1[2] + p2[2]) / 2,
                thickness=thickness,
                color=arrowhead_color
            ))
            
            # Second line from base corner 2 to tip
            segments.append(ZSegment(
                x1=corner2_2d[0],
                y1=corner2_2d[1],
                x2=tip_2d[0],
                y2=tip_2d[1],
                z_value=(corner2[2] + p2[2]) / 2,
                thickness=thickness,
                color=arrowhead_color
            ))
            
            # Base line between the base corners
            segments.append(ZSegment(
                x1=corner1_2d[0],
                y1=corner1_2d[1],
                x2=corner2_2d[0],
                y2=corner2_2d[1],
                z_value=(corner1[2] + corner2[2]) / 2,
                thickness=thickness,
                color=arrowhead_color
            ))
        
        return segments


class AxisSystem:
    """
    Represents a 3D coordinate axis system (x, y, z axes).
    """
    def __init__(self, origin=(0, 0, 0)):
        """
        Initialize a 3D axis system at the given origin.
        
        Parameters:
            origin: (x, y, z) coordinates of the axis system origin
        """
        self.origin = np.array(origin)
        self.vectors = []
        self.axis_labels = ["X", "Y", "Z"]
        self._create_axes()
    
    def _create_axes(self):
        """
        Create the three axis vectors (x, y, z).
        """
        length = AXES_STYLE["length"]
        thickness_multiplier = AXES_STYLE["thickness_multiplier"]
        
        # Create X axis (red)
        x_axis = Vector(
            start_point=self.origin,
            direction=(1, 0, 0),
            length=length,
            color=AXES_STYLE["x_color"]
        )
        
        # Create Y axis (green)
        y_axis = Vector(
            start_point=self.origin,
            direction=(0, 1, 0),
            length=length,
            color=AXES_STYLE["y_color"]
        )
        
        # Create Z axis (blue)
        z_axis = Vector(
            start_point=self.origin,
            direction=(0, 0, 1),
            length=length,
            color=AXES_STYLE["z_color"]
        )
        
        self.vectors = [x_axis, y_axis, z_axis]
    
    def get_vectors(self):
        """
        Return the vectors representing the axes.
        
        Returns:
            list: The three axis vectors
        """
        return self.vectors
    
    def get_labels(self):
        """
        Return information about axis labels (positions and text).
        
        Returns:
            list: Tuples of (position, text, color) for each axis label
        """
        if not AXES_STYLE["show_labels"]:
            return []
        
        labels = []
        label_offset = AXES_STYLE["label_offset"]
        
        for i, (vector, label_text) in enumerate(zip(self.vectors, self.axis_labels)):
            # Calculate position (a bit past the end of the vector)
            end_point = vector.get_end_point()
            label_pos = end_point + vector.direction * label_offset
            
            # Store label information
            labels.append((label_pos, label_text, vector.color))
        
        return labels