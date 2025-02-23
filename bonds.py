# bonds.py

from abc import ABC, abstractmethod
import math
from config import COVALENT_BOND, NCI_BOND, TS_BOND, DISTANCE_BOND
from geometry_utils import project_point
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
        dist = math.sqrt(vx * vx + vy * vy + vz * vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []

        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        seg_length = COVALENT_BOND["segment_length"]
        N = max(1, int(math.ceil(visible_length / seg_length)))
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
            
            bond_thickness_px = max(1, int(COVALENT_BOND["thickness"] * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=COVALENT_BOND["color"]
            )
            # Attach underlying bond data.
            seg_obj.bond = self
            segments.append(seg_obj)
        
        self.length = dist
        return segments

class NCIBond(Bond):
    """
    Represents a non-covalent interaction (NCI) bond as a dashed line.
    """
    def get_segments(self, rotated_coords, view_params):
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx * vx + vy * vy + vz * vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []

        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        seg_length = NCI_BOND["segment_length"]
        N = max(1, int(math.ceil(visible_length / seg_length)))
        dt = (t_end - t_start) / N

        segments = []
        for seg_index in range(N):
            # Skip every other segment to create a dashed effect.
            if seg_index % 2 == 1:
                continue

            t1 = t_start + seg_index * dt
            t2 = t_start + (seg_index + 1) * dt

            p1 = (x1 + vx * t1, y1 + vy * t1, z1 + vz * t1)
            p2 = (x1 + vx * t2, y1 + vy * t2, z1 + vz * t2)
            
            X1, Y1 = project_point(*p1, view_params)
            X2, Y2 = project_point(*p2, view_params)
            
            zm = 0.5 * (p1[2] + p2[2])
            
            # Apply the thickness factor for NCI bonds.
            thickness = NCI_BOND["thickness"] * NCI_BOND["thickness_factor"]
            bond_thickness_px = max(1, int(thickness * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=NCI_BOND["color"]
            )
            seg_obj.bond = self
            segments.append(seg_obj)
        
        self.length = dist
        return segments

class TS_Bond(Bond):
    """
    Represents a TS bond as a solid line using the TS_BOND settings.
    """
    def get_segments(self, rotated_coords, view_params):
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx * vx + vy * vy + vz * vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []
        
        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        seg_length = TS_BOND["segment_length"]
        N = max(1, int(math.ceil(visible_length / seg_length)))
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
            
            bond_thickness_px = max(1, int(TS_BOND["thickness"] * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=TS_BOND["color"]
            )
            seg_obj.bond = self
            segments.append(seg_obj)
        
        self.length = dist
        return segments

class Distance_Bond(Bond):
    """
    Represents a distance bond, used for tracking distances between atoms.
    Drawn as a solid line with minimal thickness.
    """
    def get_segments(self, rotated_coords, view_params):
        (x1, y1, z1, r_i), (x2, y2, z2, r_j) = rotated_coords

        vx = x2 - x1
        vy = y2 - y1
        vz = z2 - z1
        dist = math.sqrt(vx * vx + vy * vy + vz * vz)
        if dist < 1e-6 or (r_i + r_j >= dist):
            return []
        
        t_start, t_end = self.compute_visible_region(r_i, r_j, dist)
        visible_length = (t_end - t_start) * dist

        seg_length = DISTANCE_BOND["segment_length"]
        N = max(1, int(math.ceil(visible_length / seg_length)))
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
            
            bond_thickness_px = max(1, int(DISTANCE_BOND["thickness"] * view_params.scale))
            
            seg_obj = ZSegment(
                x1=X1, y1=Y1, x2=X2, y2=Y2,
                z_value=zm,
                thickness=bond_thickness_px,
                color=DISTANCE_BOND["color"]
            )
            seg_obj.bond = self
            segments.append(seg_obj)
        
        self.length = dist
        return segments
