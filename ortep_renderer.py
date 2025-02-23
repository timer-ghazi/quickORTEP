# ortep_renderer.py

import math
from geometry_utils import rotate_point, project_point
from elements_table import Elements
from zobjects import ZAtom, ZSegment  # Our drawable objects
from config import SCALE_ATOM_SPHERE, ARC_FLATTEN, BOND_THICKNESS_ANG, BOND_SEGMENT_LENGTH_ANG

def hex_to_rgb(hex_color):
    """
    Convert a hex color string (e.g. "#FF0000") to an (R,G,B) tuple.
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
    a painter's algorithm (z-sorting). Bond segmentation is delegated to the bond objects.
    """
    def build_render_list(self, ortep_molecule, vp):
        render_list = []
        rotated_info = []  # List of tuples: (x_rot, y_rot, z_rot, r_eff)

        # --- Process atoms ---
        for atom in ortep_molecule.atoms:
            # Rotate the atom's position.
            x_rot, y_rot, z_rot = rotate_point(atom.x, atom.y, atom.z, vp.rx, vp.ry, vp.rz)
            # Project to 2D screen coordinates.
            x2d, y2d = project_point(x_rot, y_rot, z_rot, vp)
            # Get the covalent radius (in Ångströms).
            try:
                r_covalent = Elements.covalent_radius(atom.symbol, order="single",
                                                       source="cordero", unit="Ang")
            except KeyError:
                r_covalent = 1.0  # Fallback if unknown
            # Compute the drawn radius in pixels.
            px_r = max(2, int(r_covalent * SCALE_ATOM_SPHERE * vp.scale))
            # Compute effective 3D radius for bond clipping (in Å).
            r_eff = r_covalent * SCALE_ATOM_SPHERE

            # Create a ZAtom for drawing.
            z_atom = ZAtom(
                x2d=x2d,
                y2d=y2d,
                radius=px_r,
                color=self._get_atom_color(atom),
                z_value=z_rot
            )
            # Attach underlying atom data and propagate its persistent selection state.
            z_atom.atom = atom
            z_atom.selected = atom.selected
            render_list.append(z_atom)
            rotated_info.append((x_rot, y_rot, z_rot, r_eff))

        # --- Process bonds ---
        for bond in ortep_molecule.bonds:
            try:
                i = ortep_molecule.atoms.index(bond.atom1)
                j = ortep_molecule.atoms.index(bond.atom2)
            except ValueError:
                continue
            rotated_coords = (rotated_info[i], rotated_info[j])
            segments = bond.get_segments(rotated_coords, vp)
            for seg in segments:
                seg.bond = bond
                # Propagate persistent selection state from the underlying bond.
                seg.selected = getattr(bond, 'selected', False)
                render_list.append(seg)
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
        Draw the molecule by building a render list, sorting by z-value,
        and drawing each object in order.
        """
        z_list = self.build_render_list(ortep_molecule, view_params)
        # Sort objects so that those with higher z_value (closer) are drawn last.
        z_list.sort(key=lambda obj: obj.z_value, reverse=True)
        for obj in z_list:
            obj.draw(canvas)
