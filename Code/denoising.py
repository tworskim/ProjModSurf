import numpy as np
import math


def calculate_dihedral_angle(edge):
    pt1 = edge[0]
    pt2 = edge[1]

    dot = np.dot(pt1, pt2)

    angle = dot / (np.linalg.norm(pt1) * np.linalg.norm(pt2))

    return math.acos(angle)


def calculate_average_dihedral_angles(edges):
    edge_angles = [calculate_dihedral_angle(edge) for edge in edges]

    return np.average(edge_angles)


