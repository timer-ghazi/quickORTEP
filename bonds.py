# bonds.py

from abc import ABC, abstractmethod
import math
from config import BOND_SEGMENT_LENGTH_ANG, BOND_THICKNESS_ANG
from geometry_utils import project_point
from zobjects import ZSegment

class Bond(ABC):
    """
    Abstract base class for bonds.
    
    This class encapsulates common bond logic, such as computing the visible
    region (based on effective atomic radii) and splitting a bond into segments.
    """

    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    @abstractmethod
    def get_segments(self, rotated_coords, view_params):
        """
        Given the rotated coordinates of both atoms and the view parameters,
        returns a list of segment descriptors (e.g. ZSegment objects) to be drawn.
        """
        pass

    def compute_visible_region(self, r_i, r_j, dist):
        """
        Compute the fraction along the bond where drawing should start and end.
        
        Uses the effective atomic radii to trim the bond:
        
        $$
        t_{start} = \\frac{r_i}{\\text{dist}}, \\quad t_{end} = 1 - \\frac{r_j}{\\text{dist}}
        $$
        
        :param r_i: Effective radius for atom1.
        :param r_j: Effective radius for atom2.
        :param dist: 3D distance between the two atoms.
        :return: Tuple (t_start, t_end) representing the normalized positions.
        """
        t_start = r_i / dist
        t_end = 1.0 - (r_j / dist)
        return t_start, t_end

    def split_into_segments(self, start_point, end_point, num_segments):
        """
        Helper function to split the bond (from start_point to end_point) into
        `num_segments` parts.
        
        :param start_point: Starting 3D coordinate (x, y, z).
        :param end_point: Ending 3D coordinate (x, y, z).
        :param num_segments: Number of segments.
        :return: List of tuples, each containing (seg_start, seg_end).
        """
        segments = []
        dt = 1.0 / num_segments
        for seg_index in range(num_segments):
            t1 = seg_index * dt
            t2 = (seg_index + 1) * dt
            seg_start = tuple(s + (e - s) * t1 for s, e in zip(start_point, end_point))
            seg_end = tuple(s + (e - s) * t2 for s, e in zip(start_point, end_point))
            segments.append((seg_start, seg_end))
        return segments


class CovalentBond(Bond):
    """
    Represents a covalent bond as a solid line.
    
    This bond type will be rendered as a series of segments (solid, black),
    where the endpoints are trimmed based on the effective atomic radii.
    """

    def get_segments(self, rotated_coords, view_params):
        """
        Compute the drawable segments for a covalent bond.
        
        :param rotated_coords: Tuple containing two tuples:
            ((x1, y1, z1, r_eff1), (x2, y2, z2, r_eff2)).
        :param view_params: View parameters (including scale) used for projection.
        :return: List of ZSegment objects.
        """
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        # Compute the vector between the atoms and its 3D distance.
        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx*vx + vy*vy + vz*vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            # No visible bond segment if atoms overlap or nearly so.
            return []

        # Compute where the bond should start and end based on effective radii.
        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        # Determine the number of segments.
        N = max(1, int(math.ceil(visible_length / BOND_SEGMENT_LENGTH_ANG)))
        dt = (t_end - t_start) / N

        segments = []
        for seg_index in range(N):
            t1 = t_start + seg_index * dt
            t2 = t_start + (seg_index + 1) * dt

            # Compute 3D coordinates for the endpoints of the segment.
            p1 = (x1 + vx * t1, y1 + vy * t1, z1 + vz * t1)
            p2 = (x1 + vx * t2, y1 + vy * t2, z1 + vz * t2)
            
            # Project the 3D points to 2D screen coordinates.
            X1, Y1 = project_point(*p1, view_params)
            X2, Y2 = project_point(*p2, view_params)
            
            # Use the midpoint's z for painter's algorithm sorting.
            zm = 0.5 * (p1[2] + p2[2])
            
            # Compute bond thickness in pixels.
            bond_thickness_px = max(1, int(BOND_THICKNESS_ANG * view_params.scale))
            
            # Create a segment object (solid, black bond).
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=(0, 0, 0)
            )
            segments.append(seg_obj)
        return segments
