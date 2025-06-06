# ortep_renderer.py

import numpy as np
from geometry_utils import rotate_points, project_points
from elements_table import Elements
from zobjects import ZAtom, ZSegment, ZArrowHead
from config import ATOM_STYLE, ARC_STYLE, AXES_STYLE

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
    Also handles rendering of vectors and coordinate axes.
    """
    def build_render_list(self, ortep_molecule, vp):
        render_list = []
        rotated_info = []  # List of tuples: (x_rot, y_rot, z_rot, r_eff)
        n_atoms = len(ortep_molecule.atoms)
        
        if n_atoms == 0:
            # If no atoms, we may still have vectors (like coordinate axes)
            if ortep_molecule.vectors:
                render_list.extend(self._build_vector_render_list(ortep_molecule, vp))
            return render_list

        # --- Vectorized processing of atom positions ---
        # Collect atom coordinates into an (n, 3) NumPy array.
        coords = np.array([[atom.x, atom.y, atom.z] for atom in ortep_molecule.atoms])
        
        # Perform vectorized rotation using the angles specified in view parameters.
        rotated_coords = rotate_points(coords, vp.rx, vp.ry, vp.rz)
        
        # Project all rotated coordinates onto 2D screen in one go.
        screen_coords = project_points(rotated_coords, vp)

        # Process each atom with its rotated and projected coordinates.
        for idx, atom in enumerate(ortep_molecule.atoms):
            x_rot, y_rot, z_rot = rotated_coords[idx]
            # Retrieve the covalent radius in Ångströms.
            try:
                r_covalent = Elements.covalent_radius(atom.symbol, order="single",
                                                       source="cordero", unit="Ang")
            except KeyError:
                r_covalent = 1.0  # Fallback if unknown

            # Compute the drawn radius in pixels.
            px_r = max(ATOM_STYLE["min_radius"], int(r_covalent * ATOM_STYLE["scale"] * vp.scale))
            # Compute effective 3D radius for bond clipping (in Å).
            r_eff = r_covalent * ATOM_STYLE["scale"]

            # Get the 2D screen coordinates from the projected results.
            x2d, y2d = screen_coords[idx]

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
            bond_rotated_coords = (rotated_info[i], rotated_info[j])
            segments = bond.get_segments(bond_rotated_coords, vp)
            for seg in segments:
                seg.bond = bond
                # Propagate persistent selection state from the underlying bond.
                seg.selected = getattr(bond, 'selected', False)
                render_list.append(seg)
        
        # --- Process vectors ---
        vector_objects = self._build_vector_render_list(ortep_molecule, vp)
        render_list.extend(vector_objects)
        
        return render_list

    def _build_vector_render_list(self, ortep_molecule, vp):
        """
        Build the render list for vectors (including coordinate axes).
        
        Parameters:
            ortep_molecule: The ORTEP_Molecule instance containing vectors
            vp: The view parameters
            
        Returns:
            list: A list of ZObject instances for the vectors
        """
        render_list = []
        
        for vector in ortep_molecule.vectors:
            # Get start and end points
            start_point = vector.start_point
            end_point = vector.get_end_point()
            
            # Rotate start and end points
            rotated_points = rotate_points(
                np.array([start_point, end_point]), 
                vp.rx, vp.ry, vp.rz
            )
            
            # Get rotated coordinates
            start_rot = rotated_points[0]
            end_rot = rotated_points[1]
            
            # Radius is 0 for vectors (no clipping)
            rotated_coords = (
                (start_rot[0], start_rot[1], start_rot[2], 0.0),
                (end_rot[0], end_rot[1], end_rot[2], 0.0)
            )
            
            # Get segments for the vector shaft and arrowhead
            segments = vector.get_segments(rotated_coords, vp)
            
            # Add segments to render list
            for seg in segments:
                seg.vector = vector  # Link back to source vector
                seg.selected = getattr(vector, 'selected', False)
                render_list.append(seg)
            
        # If axis system exists and has labels, render them
        if ortep_molecule.axis_system and AXES_STYLE["show_labels"]:
            labels = ortep_molecule.axis_system.get_labels()
            for position, text, color in labels:
                # Rotate the label position
                rotated_pos = rotate_points(np.array([position]), vp.rx, vp.ry, vp.rz)[0]
                
                # Project to screen coordinates
                screen_pos = project_points(np.array([rotated_pos]), vp)[0]
                
                # Create a render object for the label
                # Since we don't have a ZText object, we'll create a custom object
                # that draws text at the specified position
                class ZLabel(ZSegment):
                    def __init__(self, x, y, z, text, color):
                        super().__init__(x, y, x, y, z, 1, color)  # Zero-length segment
                        self.text = text
                    
                    def draw(self, canvas):
                        # Draw text at position
                        canvas.draw_text(
                            self.x1, self.y1, 
                            self.text, 
                            color=self.color, 
                            font_size=AXES_STYLE["label_font_size"]
                        )
                
                # Create and add the label
                z_label = ZLabel(
                    int(screen_pos[0]), int(screen_pos[1]), 
                    rotated_pos[2], 
                    text, color
                )
                render_list.append(z_label)
        
        return render_list

    def _get_atom_color(self, atom):
        """
        Returns the RGB color tuple for an atom using the Elements module.
        """
        try:
            # Use the color palette specified in ATOM_STYLE
            hexcol = Elements.color(atom.symbol, palette=ATOM_STYLE.get("color_palette", "Rasmol"))
            rgbcol = hex_to_rgb(hexcol)
        except KeyError:
            # Use the fallback color from ATOM_STYLE
            rgbcol = ATOM_STYLE.get("fallback_color", (200, 200, 200))
        return rgbcol

    def draw_molecule(self, canvas, ortep_molecule, view_params):
        """
        Draw the molecule by building a render list, sorting by z-value,
        and drawing each object in order with optional fog effect.
        """
        import math
        from config import FOG_STYLE
        
        z_list = self.build_render_list(ortep_molecule, view_params)
        # Sort objects so that those with higher z_value (closer) are drawn last.
        z_list.sort(key=lambda obj: obj.z_value, reverse=True)
        
        # Get fog mode from view parameters or default to FOG_STYLE
        fog_mode = view_params.fog_mode if view_params.fog_mode is not None else FOG_STYLE["mode"]
        
        # If fog is enabled, calculate z range for the fog effect
        if fog_mode > 0 and z_list:
            # Determine min/max z values in the render list
            z_min = min(obj.z_value for obj in z_list)  # Farthest objects (smallest z)
            z_max = max(obj.z_value for obj in z_list)  # Closest objects (largest z)
            z_range = z_max - z_min
            
            # Get fog density parameter from view_params if available, otherwise use default
            fog_density = view_params.fog_current_density if view_params.fog_current_density is not None else FOG_STYLE["density"]
            
            # Determine fog color - use background color if not specified
            fog_color = FOG_STYLE["color"]
            if fog_color is None:
                fog_color = canvas.background_color
                
        # Draw all objects with fog effect applied if enabled
        for obj in z_list:
            # Apply fog if enabled
            if fog_mode > 0 and z_range > 0:
                # Calculate fog factor based on the object's z_value using exponential fog
                # Normalize Z based on distance from the *farthest* point (z_min)
                # normalized_z = 1 for farthest (z=z_min), decreases towards 0 for closest (z=z_max)
                normalized_z = max(0, (obj.z_value - z_min) / z_range) if z_range > 0 else 0
                
                # Use exponential falloff: 1 - e^(-density * normalized_distance_from_viewer^2)
                # Invert normalized_z so farther objects get more fog
                # normalized_z = 1.0 - normalized_z
                factor = 1.0 - math.exp(-(fog_density * normalized_z)**2)

                # Clamp between 0 and 1 (should already be in this range, but good practice)
                factor = max(0, min(1, factor))

                # --- Interpolate between object's color and fog color ---
                if hasattr(obj, 'color'):
                    r1, g1, b1 = obj.color
                    r2, g2, b2 = fog_color
                    # Linear interpolation: final = original + factor * (target - original)
                    final_color = (
                        int(r1 + factor * (r2 - r1)),
                        int(g1 + factor * (g2 - g1)),
                        int(b1 + factor * (b2 - b1))
                    )
                    # Draw with fogged color
                    obj.draw(canvas, color_override=final_color)
                else:
                    # Object doesn't have a color attribute, draw normally
                    obj.draw(canvas)
            else:
                # No fog, draw normally
                obj.draw(canvas)        
