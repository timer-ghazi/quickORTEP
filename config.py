# config.py

# Default conversion factor: how many pixels represent 1 Ångström.
ANGSTROM_TO_PIXEL = 100

# Atom sphere scaling factor: scales the atom’s covalent radius (in Å) to its visual size.
# For example, a covalent radius of $r_{covalent}$ will be scaled as:
#   r_scaled = r_covalent * SCALE_ATOM_SPHERE
SCALE_ATOM_SPHERE = 0.4

# Arc flattening factor: used for drawing ORTEP arcs to simulate the sphere’s curvature.
# It controls the ratio between the minor and major axes of the arc ellipses.
ARC_FLATTEN = 0.4

# Bond thickness in ang
BOND_THICKNESS_ANG = 0.18
