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
        for all atoms and bonds. For each atom, we:
          - Rotate its 3D coordinates.
          - Project to 2D.
          - Compute the drawn radius (px_r) using a fixed multiplier (40).
          - Derive an effective 3D radius for bond clipping:
            
              r_eff = (r_ang * 40) / vp.scale

        For each bond, we trim its ends so that the bond does not intrude
        into the atomic sphere (using the effective radii) and subdivide
        the remaining line into several segments.
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
            
            ## Compute the drawn radius in pixels using the configurable constants:
            ##   px_r = max(2, int(r_covalent * SCALE_ATOM_SPHERE * ANGSTROM_TO_PIXEL))
            #px_r = max(2, int(r_covalent * SCALE_ATOM_SPHERE * ANGSTROM_TO_PIXEL))
            #
            ## Compute the effective 3D radius for bond clipping.
            ##   r_eff = (r_covalent * SCALE_ATOM_SPHERE * ANGSTROM_TO_PIXEL) / vp.scale
            #r_eff = (r_covalent * SCALE_ATOM_SPHERE * ANGSTROM_TO_PIXEL) / vp.scale

            px_r = max(2, int(r_covalent * SCALE_ATOM_SPHERE * vp.scale))
            r_eff = r_covalent * SCALE_ATOM_SPHERE  # since (r_covalent * SCALE_ATOM_SPHERE * vp.scale) / vp.scale simplifies

            # Create a ZAtom to be drawn
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

            # Find indices so we can retrieve the precomputed rotated info.
            try:
                i = ortep_molecule.atoms.index(ai)
                j = ortep_molecule.atoms.index(aj)
            except ValueError:
                # If for some reason an atom isn't in the list, skip this bond.
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

            # If the effective spheres overlap or just touch, skip drawing the bond.
            if r_i + r_j >= dist:
                continue

            # Compute the parameters (t values) where the bond becomes visible.
            # We clip t from 0 to 1 so that:
            #   - The bond starts at a distance r_i from atom i,
            #   - and ends r_j before atom j.
            t_start = r_i / dist
            t_end   = 1.0 - (r_j / dist)

            # Subdivide the visible portion into N segments for smooth occlusion.
            N = 6
            dt = (t_end - t_start) / N
            if dt <= 0:
                continue  # Should not happen, but safety first

            for seg_index in range(N):
                # Determine parameters for the segment endpoints
                t1 = t_start + seg_index * dt
                t2 = t_start + (seg_index + 1) * dt

                # Compute 3D endpoints for this bond segment
                x1_3d = x_i + vx * t1
                y1_3d = y_i + vy * t1
                z1_3d = z_i + vz * t1

                x2_3d = x_i + vx * t2
                y2_3d = y_i + vy * t2
                z2_3d = z_i + vz * t2

                # Use the midpoint for z-sorting
                xm = 0.5 * (x1_3d + x2_3d)
                ym = 0.5 * (y1_3d + y2_3d)
                zm = 0.5 * (z1_3d + z2_3d)

                # Project the endpoints to 2D
                X1, Y1 = project_point(x1_3d, y1_3d, z1_3d, vp)
                X2, Y2 = project_point(x2_3d, y2_3d, z2_3d, vp)

                # Compute bond thickness in pixels by multiplying your chosen Å thickness by vp.scale.
                # Also ensure a minimum of 1 pixel so bonds never disappear at small zoom.
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

                # # Create a ZSegment for this portion of the bond.
                # seg_obj = ZSegment(
                #     x1=X1,
                #     y1=Y1,
                #     x2=X2,
                #     y2=Y2,
                #     z_value=zm,
                #     thickness=16,
                #     color=(0, 0, 0)
                # )
                render_list.append(seg_obj)

        return render_list

    def _get_atom_color(self, atom):
        """
        Helper: returns the RGB color tuple for an atom using the Elements module.
        """
        try:
            hexcol = Elements.color(atom.symbol, palette="Rasmol")
            rgbcol = hex_to_rgb(hexcol)
        except KeyError:
            rgbcol = (200, 200, 200)
        return rgbcol

    def draw_molecule(self, canvas, ortep_molecule, view_params):
        """
        High-level method to draw the molecule: we build a list of
        ZObjects, sort them by z_value (farthest first), and then
        draw them in order.
        """
        # 1) Build the render list (atoms + bond segments)
        z_list = self.build_render_list(ortep_molecule, view_params)

        # 2) Sort by z_value descending (farthest first)
        z_list.sort(key=lambda obj: obj.z_value, reverse=True)

        # 3) Draw all objects
        for obj in z_list:
            obj.draw(canvas)
