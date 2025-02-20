import math
from geometry_utils import rotate_point, project_point
from elements_table import Elements
from zobjects import ZAtom, ZSegment  # Our new drawable objects

from config import SCALE_ATOM_SPHERE, ANGSTROM_TO_PIXEL, BOND_THICKNESS_ANG

def hex_to_rgb(hex_color):
    """
    Helper to convert a #RRGGBB string into an (R,G,B) tuple.
    Fallback to (200,200,200) if invalid.
    """
    if not (hex_color.startswith('#') and len(hex_color) == 7):
        return (200, 200, 200)
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    return (r, g, b)

class ORTEP_MoleculeRenderer:
    """
    Responsible for drawing atoms and bonds of an ORTEP_Molecule using
    a painter's algorithm (z-sorting). Bonds are split into segments that
    are trimmed to the edges of the atomic spheres (using the actual drawn radii).
    """

    def build_render_list(self, ortep_molecule, vp):
        """
        Build a list of ZObject instances (ZAtom, ZSegment, etc.)
        for all atoms and bonds.
        
        For each atom:
        - Rotate its 3D coordinates.
        - Project to 2D.
        - Compute the drawn radius (px_r) using the configurable scale.
        - Compute an effective 3D radius for bond clipping:
        
          $r_{\text{eff}} = r_{\text{covalent}} \times \text{SCALE_ATOM_SPHERE}$
        
        For each bond:
        - Trim its ends so the bond doesn't intrude into the atomic sphere.
        - Subdivide the remaining line into several segments.
        - Unify endpoints across segments to avoid gaps.
        """
        render_list = []
        rotated_info = []  # List of (x_rot, y_rot, z_rot, r_eff)

        # --- 1) Process atoms ---
        for atom in ortep_molecule.atoms:
            # Rotate the atom's position
            x_rot, y_rot, z_rot = rotate_point(atom.x, atom.y, atom.z,
                                               vp.rx, vp.ry, vp.rz)
            # Project to 2D screen coordinates
            x2d, y2d = project_point(x_rot, y_rot, z_rot, vp)

            # Get the covalent radius (in Angstroms)
            try:
                r_covalent = Elements.covalent_radius(atom.symbol, order="single",
                                                       source="cordero", unit="Ang")
            except KeyError:
                r_covalent = 1.0  # Fallback if unknown
            
            # Compute the drawn radius in pixels
            px_r = max(2, int(r_covalent * SCALE_ATOM_SPHERE * vp.scale))
            # Compute effective 3D radius for bond clipping (in Å)
            r_eff = r_covalent * SCALE_ATOM_SPHERE

            # Create a ZAtom for drawing
            z_atom = ZAtom(
                x2d=x2d,
                y2d=y2d,
                radius=px_r,
                color=self._get_atom_color(atom),
                z_value=z_rot
            )
            render_list.append(z_atom)

            # Save rotated position and effective radius for bond clipping
            rotated_info.append((x_rot, y_rot, z_rot, r_eff))

        # --- 2) Process bonds ---
        for bond in ortep_molecule.bonds:
            ai = bond.atom1
            aj = bond.atom2

            # Find indices to retrieve precomputed rotated info.
            try:
                i = ortep_molecule.atoms.index(ai)
                j = ortep_molecule.atoms.index(aj)
            except ValueError:
                continue

            x_i, y_i, z_i, r_i = rotated_info[i]
            x_j, y_j, z_j, r_j = rotated_info[j]

            # Vector from atom i to atom j in rotated 3D space
            vx = x_j - x_i
            vy = y_j - y_i
            vz = z_j - z_i
            dist = math.sqrt(vx*vx + vy*vy + vz*vz)
            if dist < 1e-6:
                continue  # Avoid division by zero

            # Skip bond if effective spheres overlap or touch
            if r_i + r_j >= dist:
                continue

            # Compute t values for bond visibility:
            # Bond starts at a distance r_i from atom i and ends r_j before atom j.
            t_start = r_i / dist
            t_end = 1.0 - (r_j / dist)

            # Subdivide the visible portion into N segments
            N = 6
            dt = (t_end - t_start) / N
            if dt <= 0:
                continue

            # Generate N+1 points along the bond in 3D
            points_3d = []
            for seg_index in range(N+1):
                t = t_start + seg_index * dt
                x_3d = x_i + vx * t
                y_3d = y_i + vy * t
                z_3d = z_i + vz * t
                points_3d.append((x_3d, y_3d, z_3d))

            # Determine coordinate type based on view parameters.
            # For SVG export, we want float coordinates.
            as_float = getattr(vp, "as_float", False)

            # Project all 3D points to 2D once
            projected_points = []
            for (px, py, pz) in points_3d:
                X, Y = project_point(px, py, pz, vp, as_float=as_float)
                projected_points.append((X, Y, pz))

            # Build bond segments using consecutive projected points
            for seg_index in range(N):
                (X1, Y1, Z1) = projected_points[seg_index]
                (X2, Y2, Z2) = projected_points[seg_index+1]

                # Use the midpoint z for painter's sorting
                zm = 0.5 * (Z1 + Z2)

                # Compute bond thickness in pixels
                bond_thickness_px = max(1, int(BOND_THICKNESS_ANG * vp.scale))
                
                seg_obj = ZSegment(
                    x1=X1,
                    y1=Y1,
                    x2=X2,
                    y2=Y2,
                    z_value=zm,
                    thickness=bond_thickness_px,  
                    color=(0, 0, 0)
                )
                render_list.append(seg_obj)

        return render_list

    def _get_atom_color(self, atom):
        """
        Returns the RGB color tuple for an atom using the Elements module.
        """
        try:
            hexcol = Elements.color(atom.symbol, palette="Rasmol")
            rgbcol = hex_to_rgb(hexcol)
        except KeyError:
            rgbcol = (200, 200, 200)
        return rgbcol

    def draw_molecule(self, canvas, ortep_molecule, view_params):
        """
        Draw the molecule by building a render list (ZObjects), sorting them by z-value,
        and drawing them in order.
        """
        # 1) Build the render list (atoms + bond segments)
        z_list = self.build_render_list(ortep_molecule, view_params)

        # 2) Sort objects by z-value (farthest first)
        z_list.sort(key=lambda obj: obj.z_value, reverse=True)

        # 3) Draw each object
        for obj in z_list:
            obj.draw(canvas)
