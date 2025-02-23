# bonds.py

from abc import ABC, abstractmethod
import math
from config import BOND_SEGMENT_LENGTH_ANG, BOND_THICKNESS_ANG
from geometry_utils import project_point
from zobjects import ZSegment

class Bond(ABC):
    """
    Abstract base class for bonds.
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.selected = False  # New attribute to persist selection state

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
        
        $$ t_{start} = \\frac{r_i}{dist}, \\quad t_{end} = 1 - \\frac{r_j}{dist} $$
        """
        t_start = r_i / dist
        t_end = 1.0 - (r_j / dist)
        return t_start, t_end

    def split_into_segments(self, start_point, end_point, num_segments):
        """
        Helper function to split the bond into num_segments segments.
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
    """
    def get_segments(self, rotated_coords, view_params):
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        # Compute the vector between the atoms.
        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx*vx + vy*vy + vz*vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []

        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        N = max(1, int(math.ceil(visible_length / BOND_SEGMENT_LENGTH_ANG)))
        dt = (t_end - t_start) / N

        segments = []
        for seg_index in range(N):
            t1 = t_start + seg_index * dt
            t2 = t_start + (seg_index + 1) * dt

            p1 = (x1 + vx * t1, y1 + vy * t1, z1 + vz * t1)
            p2 = (x1 + vx * t2, y1 + vy * t2, z1 + vz * t2)
            
            X1, Y1 = project_point(*p1, view_params)
            X2, Y2 = project_point(*p2, view_params)
            
            zm = 0.5 * (p1[2] + p2[2])
            
            bond_thickness_px = max(1, int(BOND_THICKNESS_ANG * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=(0, 0, 0)
            )
            # Attach underlying bond data.
            seg_obj.bond = self
            segments.append(seg_obj)
        # Store the bond length for later use.
        self.length = dist
        return segments

class NCIBond(Bond):
    """
    Represents a non-covalent interaction (NCI) bond as a dashed grey line.
    """
    def get_segments(self, rotated_coords, view_params):
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx*vx + vy*vy + vz*vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []

        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        N = max(1, int(math.ceil(visible_length / BOND_SEGMENT_LENGTH_ANG)))
        dt = (t_end - t_start) / N

        segments = []
        for seg_index in range(N):
            if seg_index % 2 == 1:
                continue  # Skip every other segment to create dashed effect

            t1 = t_start + seg_index * dt
            t2 = t_start + (seg_index + 1) * dt

            p1 = (x1 + vx * t1, y1 + vy * t1, z1 + vz * t1)
            p2 = (x1 + vx * t2, y1 + vy * t2, z1 + vz * t2)
            
            X1, Y1 = project_point(*p1, view_params)
            X2, Y2 = project_point(*p2, view_params)
            
            zm = 0.5 * (p1[2] + p2[2])
            
            bond_thickness_px = max(1, int((BOND_THICKNESS_ANG * 0.5) * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=(128, 128, 128)
            )
            # Attach underlying bond data.
            seg_obj.bond = self
            segments.append(seg_obj)
        self.length = dist
        return segments
