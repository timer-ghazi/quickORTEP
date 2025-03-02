import numpy as np
import math  # Keep for backward compatibility

def rotate_points(points, rx_deg, ry_deg, rz_deg):
    """
    Rotate an array of 3D points by rx_deg, ry_deg, rz_deg about X, Y, Z axes, respectively.
    
    Parameters:
    -----------
    points : numpy.ndarray
        Array of shape (n, 3) containing n points with (x, y, z) coordinates
    rx_deg, ry_deg, rz_deg : float
        Rotation angles in degrees around X, Y, and Z axes
        
    Returns:
    --------
    numpy.ndarray
        Array of shape (n, 3) containing the rotated points
    """
    # Convert angles to radians
    rx = np.radians(rx_deg)
    ry = np.radians(ry_deg)
    rz = np.radians(rz_deg)
    
    # Create rotation matrices
    # Rotation around X axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rx), -np.sin(rx)],
        [0, np.sin(rx), np.cos(rx)]
    ])
    
    # Rotation around Y axis
    Ry = np.array([
        [np.cos(ry), 0, np.sin(ry)],
        [0, 1, 0],
        [-np.sin(ry), 0, np.cos(ry)]
    ])
    
    # Rotation around Z axis
    Rz = np.array([
        [np.cos(rz), -np.sin(rz), 0],
        [np.sin(rz), np.cos(rz), 0],
        [0, 0, 1]
    ])
    
    # Combine rotation matrices - order matters: first X, then Y, then Z
    R = Rz @ Ry @ Rx
    
    # Apply rotation to all points at once
    rotated_points = points @ R.T
    
    return rotated_points

def project_points(points, view_params, as_float=False):
    """
    Project an array of 3D points onto a 2D screen using the scale and offsets.
    
    Parameters:
    -----------
    points : numpy.ndarray
        Array of shape (n, 3) containing n points with (x, y, z) coordinates
    view_params : object
        Object with attributes: scale, x_offset, y_offset
    as_float : bool, optional
        If True, return floating-point coordinates, otherwise integer coordinates
        
    Returns:
    --------
    numpy.ndarray
        Array of shape (n, 2) containing the screen coordinates
    """
    # Apply scale and offset (only to x and y coordinates)
    screen_coords = np.zeros((len(points), 2))
    screen_coords[:, 0] = points[:, 0] * view_params.scale + view_params.x_offset
    screen_coords[:, 1] = points[:, 1] * view_params.scale + view_params.y_offset
    
    # Convert to integers if needed
    if not as_float:
        screen_coords = screen_coords.astype(int)
    
    return screen_coords

def rotate_and_project_points(points, rx_deg, ry_deg, rz_deg, view_params, as_float=False):
    """
    Rotate and project an array of 3D points in one operation.
    
    Parameters:
    -----------
    points : numpy.ndarray
        Array of shape (n, 3) containing n points with (x, y, z) coordinates
    rx_deg, ry_deg, rz_deg : float
        Rotation angles in degrees around X, Y, and Z axes
    view_params : object
        Object with attributes: scale, x_offset, y_offset
    as_float : bool, optional
        If True, return floating-point coordinates, otherwise integer coordinates
        
    Returns:
    --------
    numpy.ndarray
        Array of shape (n, 2) containing the screen coordinates
    """
    rotated_points = rotate_points(points, rx_deg, ry_deg, rz_deg)
    return project_points(rotated_points, view_params, as_float)


# Backward compatibility wrappers for the original API
def rotate_point(x, y, z, rx_deg, ry_deg, rz_deg):
    """
    Rotate the 3D point (x, y, z) by rx_deg, ry_deg, rz_deg about X, Y, Z axes, respectively.
    Returns the rotated coordinates (x_rot, y_rot, z_rot).
    
    This is a wrapper for the vectorized rotate_points function.
    """
    point = np.array([[x, y, z]])
    rotated = rotate_points(point, rx_deg, ry_deg, rz_deg)
    return tuple(rotated[0])

def project_point(x, y, z, view_params, as_float=False):
    """
    Takes the rotated 3D point (x, y, z) and projects it onto the 2D screen
    using the scale and offsets in view_params. Returns (X_screen, Y_screen).
    If as_float=True, returns floating-point coords. Otherwise, integer coords.
    
    This is a wrapper for the vectorized project_points function.
    """
    point = np.array([[x, y, z]])
    projected = project_points(point, view_params, as_float)
    return tuple(projected[0])