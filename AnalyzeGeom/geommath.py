import numpy as np
from typing import List
from numpy.typing import NDArray

def magnitude(u: NDArray[np.float64]) -> float:
    """
    Calculate the magnitude (Euclidean norm) of vector u.
    
    Parameters
    ----------
    u : List[float]
        Vector of interest.
        
    Returns
    -------
    float
        Magnitude of the vector.
    """
    magnitude = np.sqrt(u.dot(u))
    return magnitude
    
def normalize(u):
    """ Normalize a vector to have a unit length. """
    return u / np.linalg.norm(u)

def vector(u):
    return np.array(u)

def distance(r1, r2):
    """
    Calculate the Euclidean distance between two 3D points.
    
    Parameters
    ----------
    r1 : List[float]
        First 3D point.
    r2 : List[float]
        Second 3D point.
        
    Returns
    -------
    float
        Euclidean distance between the points.
    """
    return np.sqrt(sum([(a - b) ** 2 for a, b in zip(r1, r2)]))


def vector_sum(vectors):
    """
    Calculate the element-wise sum of multiple vectors.
    
    Parameters
    ----------
    vectors : List[List[float]]
        List of vectors to sum.
        
    Returns
    -------
    List[float]
        Summed vector.
    """
    return [sum(x) for x in zip(*vectors)]

def scaled_distance(point_fixed, point_scaled, desired_distance):
    """
    Calculate a new point at a specific distance from the fixed point along the direction of another point.

    Parameters
    ----------
    point_fixed : np.ndarray
        The starting (fixed) point.
    point_scaled : np.ndarray
        The point used to define the direction.
    desired_distance : float
        The distance to move along the direction from point_fixed.
        
    Returns
    -------
    np.ndarray
        The new point at the specified distance from point_fixed.
    """
    direction_vector = point_scaled - point_fixed
    unit_vector = normalize(direction_vector)
    scaled_vector = unit_vector * desired_distance
    new_point = point_fixed + scaled_vector
    return new_point

def angle_between_vectors(v1, v2):
    """Calculate the angle in radians between two vectors."""
    v1_u = normalize(v1)
    v2_u = normalize(v2)
    dot_product = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)  # Ensure within [-1, 1]
    return np.arccos(dot_product)


def rotation_matrix(axis, angle):
    """
    Create a rotation matrix for rotating around a given axis by a specified angle.
    
    Parameters:
    - axis: A 3-element array specifying the axis of rotation.
    - angle: The angle in radians by which to rotate.
    
    Returns:
    - A 3x3 rotation matrix.
    """
    axis = normalize(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c + a*d), 2*(b*d - a*c)],
        [2*(b*c - a*d), a*a + c*c - b*b - d*d, 2*(c*d + a*b)],
        [2*(b*d + a*c), 2*(c*d - a*b), a*a + d*d - b*b - c*c]
    ])

def rotate_coordinates(coords, center, rotation_matrix):
    """
    Rotate a set of coordinates around a specified center using a rotation matrix.
    
    Parameters:
    - coords: An Nx3 array of coordinates to rotate.
    - center: The center of rotation as a 3-element array.
    - rotation_matrix: The 3x3 rotation matrix.
    
    Returns:
    - The rotated Nx3 array of coordinates.
    """
    # Translate coordinates to origin
    translated_coords = coords - center
    # Apply the rotation
    rotated_coords = np.dot(translated_coords, rotation_matrix.T)
    # Translate back to original position
    final_coords = rotated_coords + center
    return final_coords

def translate_bysettingorigin(coords, idx_toorigin):
    translated_coords = coords - coords[idx_toorigin]
    return translated_coords
