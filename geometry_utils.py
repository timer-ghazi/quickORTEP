import math

def rotate_point(x, y, z, rx_deg, ry_deg, rz_deg):
    """
    Rotate the 3D point (x, y, z) by rx_deg, ry_deg, rz_deg about X, Y, Z axes, respectively.
    Returns the rotated coordinates (x_rot, y_rot, z_rot).
    """
    rx = math.radians(rx_deg)
    ry = math.radians(ry_deg)
    rz = math.radians(rz_deg)

    cosx = math.cos(rx); sinx = math.sin(rx)
    cosy = math.cos(ry); siny = math.sin(ry)
    cosz = math.cos(rz); sinz = math.sin(rz)

    # Rotate around X
    y2 = y*cosx - z*sinx
    z2 = y*sinx + z*cosx

    # Rotate around Y
    x3 = x*cosy + z2*siny
    z3 = -x*siny + z2*cosy

    # Rotate around Z
    x4 = x3*cosz - y2*sinz
    y4 = x3*sinz + y2*cosz

    return (x4, y4, z3)


def project_point(x, y, z, view_params, as_float=False):
    """
    Takes the rotated 3D point (x, y, z) and projects it onto the 2D screen
    using the scale and offsets in view_params. Returns (X_screen, Y_screen).
    If as_float=True, returns floating-point coords. Otherwise, integer coords.
    """
    X_screen = x * view_params.scale + view_params.x_offset
    Y_screen = y * view_params.scale + view_params.y_offset

    if not as_float:
        X_screen = int(X_screen)
        Y_screen = int(Y_screen)

    return (X_screen, Y_screen)
