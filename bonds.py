# bonds.py

from abc import ABC, abstractmethod
import math
import numpy as np
from config import COVALENT_BOND, NCI_BOND, TS_BOND, DISTANCE_BOND
from geometry_utils import project_point, project_points
from zobjects import ZSegment

class Bond(ABC):
    """
    Abstract base class for bonds.
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.selected = False  # Persist selection state

    @abstractmethod
    def get_segments(self, rotated_coords, view_params):
        """
        Given the rotated coordinates of both atoms and the view parameters,
        returns a list of ZSegment objects to be drawn.
        """
        pass

    def compute_visible_region(self, r_i, r_j, dist):
        """
        Compute the fraction along the bond where drawing should start and end.
        
        $$ t_{start} = \\frac{r_i}{\\text{dist}}, \\quad t_{end} = 1 - \\frac{r_j}{\\text{dist}} $$
        """
        t_start = r_i / dist
        t_end = 1.0 - (r_j / dist)
        return t_start, t_end

    def _get_vectorized_segments(self, rotated_coords, view_params, bond_config, dash=False):
        """
        Compute segments using NumPy for vectorized operations.
        
        Parameters:
          rotated_coords: tuple containing ((x1, y1, z1, r_i), (x2, y2, z2, r_j))
          view_params: view parameters (should include a 'scale' attribute)
          bond_config: a dict (e.g., COVALENT_BOND, NCI_BOND, etc.) with keys:
            'segment_length', 'thickness', 'color', and optionally 'thickness_factor'
          dash: if True, only every other segment is used (for dashed bonds).
        """
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        # Create NumPy arrays for the start and end points.
        p1 = np.array([x1, y1, z1])
        p2 = np.array([x2, y2, z2])
        v = p2 - p1
        dist = np.linalg.norm(v)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []

        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        seg_length = bond_config["segment_length"]
        N = max(1, int(math.ceil(visible_length / seg_length)))
        dt = (t_end - t_start) / N

        # Generate an array of t values.
        t_values = np.linspace(t_start, t_end, N + 1)
        # Compute the segment endpoints in a vectorized fashion.
        p1_segments = p1 + (t_values[:-1])[:, None] * v  # shape (N, 3)
        p2_segments = p1 + (t_values[1:])[:, None] * v     # shape (N, 3)

        # Compute midpoints for z-ordering.
        z_mid = 0.5 * (p1_segments[:, 2] + p2_segments[:, 2])

        # Use vectorized projection:
        # Rather than using a list comprehension, use the vectorized project_points function
        proj1 = project_points(p1_segments, view_params, as_float=view_params.as_float)
        proj2 = project_points(p2_segments, view_params, as_float=view_params.as_float)
        X1, Y1 = proj1[:, 0], proj1[:, 1]
        X2, Y2 = proj2[:, 0], proj2[:, 1]

        # For dashed bonds, select only even-indexed segments.
        if dash:
            indices = np.arange(N)[::2]
            X1, Y1 = X1[indices], Y1[indices]
            X2, Y2 = X2[indices], Y2[indices]
            z_mid = z_mid[indices]
            N = len(indices)

        # Calculate the thickness in pixels.
        thickness_factor = bond_config.get("thickness_factor", 1.0)
        bond_thickness_px = max(1, int(bond_config["thickness"] * thickness_factor * view_params.scale))

        segments = []
        for i in range(N):
            seg_obj = ZSegment(
                x1=X1[i],
                y1=Y1[i],
                x2=X2[i],
                y2=Y2[i],
                z_value=z_mid[i],
                thickness=bond_thickness_px,
                color=bond_config["color"]
            )
            seg_obj.bond = self  # Attach underlying bond data.
            segments.append(seg_obj)

        self.length = dist
        return segments

class CovalentBond(Bond):
    """
    Represents a covalent bond as a solid line.
    """
    def get_segments(self, rotated_coords, view_params):
        return self._get_vectorized_segments(rotated_coords, view_params, COVALENT_BOND)

class NCIBond(Bond):
    """
    Represents a non-covalent interaction (NCI) bond as a dashed line.
    """
    def get_segments(self, rotated_coords, view_params):
        return self._get_vectorized_segments(rotated_coords, view_params, NCI_BOND, dash=True)

class TS_Bond(Bond):
    """
    Represents a TS bond as a solid line using the TS_BOND settings.
    """
    def get_segments(self, rotated_coords, view_params):
        return self._get_vectorized_segments(rotated_coords, view_params, TS_BOND)

class Distance_Bond(Bond):
    """
    Represents a distance bond, used for tracking distances between atoms.
    Drawn as a solid line with minimal thickness.
    """
    def get_segments(self, rotated_coords, view_params):
        return self._get_vectorized_segments(rotated_coords, view_params, DISTANCE_BOND)
