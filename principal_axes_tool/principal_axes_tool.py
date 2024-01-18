import numpy as np
from pymol import cmd
from pymol.cgo import *

"""
--- Principal Axes Tool ---
Author  : Michele Bonus
Program : Principal Axes Tool aka zjob
Date    : May 2021
To Do   : Add functionality to align_axes to calculate the principal axes
          on a different atom set than the one that is translated
"""

__author__ = "Michele Bonus"
__copyright__ = "Copyright 2020, Michele Bonus"
__credits__ = ["Michele Bonus"]
__license__ = "GPL"
__version__ = "0.93"
__maintainer__ = "Michele Bonus"
__email__ = "Michele.Bonus@hhu.de"
__status__ = "Development"


def rotation_matrix_from_vectors(vec1, vec2):
    """Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def calculate_principal_axes_and_cog(selection, state):
    coords = cmd.get_coords(selection, state)
    cog = np.mean(coords, 0)
    coords = coords - cog

    inertia = np.dot(coords.transpose(), coords)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()

    axes = np.array([axis1, axis2, axis3])
    evals = np.array([eval1, eval2, eval3])
    return axes, cog


def draw_principal_axes(axes_and_cog, scale_factor):

    point1 = 3 * scale_factor * axes_and_cog[0][0] + axes_and_cog[1]
    point2 = 2 * scale_factor * axes_and_cog[0][1] + axes_and_cog[1]
    point3 = 1 * scale_factor * axes_and_cog[0][2] + axes_and_cog[1]

    axis1_cgo = [BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, VERTEX, axes_and_cog[1][0], axes_and_cog[1][1], axes_and_cog[1][2], VERTEX, point1[0], point1[1], point1[2], END]
    axis2_cgo = [BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, VERTEX, axes_and_cog[1][0], axes_and_cog[1][1], axes_and_cog[1][2], VERTEX, point2[0], point2[1], point2[2], END]
    axis3_cgo = [BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, VERTEX, axes_and_cog[1][0], axes_and_cog[1][1], axes_and_cog[1][2], VERTEX, point3[0], point3[1], point3[2], END]

    cmd.load_cgo(axis1_cgo, "axis1")
    cmd.load_cgo(axis2_cgo, "axis2")
    cmd.load_cgo(axis3_cgo, "axis3")
    cmd.set("cgo_line_width", 4)

    print("\n")
    print("First principal axis (in red)")
    print("coordinates: ", axes_and_cog[0][0])
    print("\n")
    print("Second principal axis (in green)")
    print("coordinates:", axes_and_cog[0][1])
    print("\n")
    print("Third principal axis (in blue)")
    print("coordinates:", axes_and_cog[0][2])


def draw_coordinate_axes(scale_factor=20.0):
    """
    DESCRIPTION

        "draw_coordinate_axes" draws the coordinate axes x, y, z at
        the origin.

    USAGE

        align_axes scale_factor

    ARGUMENTS

        scale_factor = float: the scaling factor for the cgo axes
                              representation {default: 20.0}

    PYMOL API

        not implemented

    EXAMPLES

        draw_coordinate_axes 20.0

    NOTES

        To be made.
    """
    scale_factor = float(scale_factor)

    x_cgo = [BEGIN, LINES, COLOR, 0.5, 0.5, 0.5, VERTEX, 0.0, 0.0, 0.0, VERTEX, 1.0 * scale_factor, 0.0, 0.0, END]
    y_cgo = [BEGIN, LINES, COLOR, 0.5, 0.5, 0.5, VERTEX, 0.0, 0.0, 0.0, VERTEX, 0.0, 1.0 * scale_factor, 0.0, END]
    z_cgo = [BEGIN, LINES, COLOR, 0.5, 0.5, 0.5, VERTEX, 0.0, 0.0, 0.0, VERTEX, 0.0, 0.0, 1.0 * scale_factor, END]

    cmd.load_cgo(x_cgo, "x_axis")
    cmd.load_cgo(y_cgo, "y_axis")
    cmd.load_cgo(z_cgo, "z_axis")
    cmd.set("cgo_line_width", 4)


def align_axes(selection1="all", selection2="all", state=1, principal_axis=1, coordinate_axis="z", vector=None, translate_to_origin=1):
    """
    DESCRIPTION

        "align_axes" aligns one of the principal axes of a selection
        to a either a coordinate axis (x, y, or z) or an axis given
        by a float vector.

    USAGE

        align_axes selection1 [, selection2 [, state [, principal_axis [, coordinate_axis [, vector [, translate_to_origin ]]]]]]

    ARGUMENTS

        selection1 = string: atoms to calculate principal axes on
                             {default: all}

        selection2 = string: atoms for which transformation should be
                             performed
                             {default: all}

        state = integer: state of the (selection in the) object for
                         which transformation should be performed
                         {default: 1}

        principal_axis = integer: the principal axis that should be
                                  aligned to the coordinate axis or
                                  vector
                                  {default: 1}

        coordinate_axis = string: the coordinate axis to be aligned
                                  to; if used, do not specify vector
                                  {default: z}

        vector = float vector: the vector of the axis to be aligned
                               to; if used, do not specify coordinate_axis
                               {default: None}

        translate_to_origin: whether or not a translation to the origin
                             should be performed after rotation.
                             {default: 1}

    PYMOL API

        not implemented

    EXAMPLES

        align_axes 5H3O and name CA, 5H3O, state=1, principal_axis = 1, coordinate_axis = z, translate_to_origin = 1
        align_axes 5H3O and name CA, 5H3O, state=1, principal_axis = 1, vector = [0.70,0.30,1.35], translate_to_origin = 0

    NOTES

        To be made.
    """
    selection1 = str(selection1)
    selection2 = str(selection2)
    state = int(state)
    principal_axis = int(principal_axis)
    coordinate_axis = str(coordinate_axis)
    translate_to_origin = int(translate_to_origin)

    if vector is not None:
        axis_vector = [float(string) for string in vector.strip("][").split(",")]
    elif coordinate_axis == "x":
        axis_vector = [1.0, 0.0, 0.0]
    elif coordinate_axis == "y":
        axis_vector = [0.0, 1.0, 0.0]
    elif coordinate_axis == "z":
        axis_vector = [0.0, 0.0, 1.0]

    axes_and_cog = calculate_principal_axes_and_cog(selection1, state)

    if principal_axis == 1:
        mat = rotation_matrix_from_vectors(axes_and_cog[0][0], axis_vector)
    elif principal_axis == 2:
        mat = rotation_matrix_from_vectors(axes_and_cog[0][1], axis_vector)
    elif principal_axis == 3:
        mat = rotation_matrix_from_vectors(axes_and_cog[0][2], axis_vector)

    mat = np.hstack((mat, np.zeros((mat.shape[0], 1))))
    print(translate_to_origin)
    if translate_to_origin == 1:
        mat = np.vstack((mat, np.append(-axes_and_cog[1], 0)))
    else:
        mat = np.vstack((mat, np.zeros((mat.shape[1]))))
    print(mat)
    mat = list(np.concatenate(mat))
    print(mat)

    cmd.transform_selection(selection2, mat, state)


def get_principal_axes(selection="all", state=1, scale_factor=20.0):
    """
    DESCRIPTION

        "get_principal_axes" calculates and displays the three principal
        axes in the selection, a particular object or object-state.

    USAGE

        get_principal_axes selection [, state [, scale_factor ]]

    ARGUMENTS

        selection = string: atoms for which principal axes should be
                            calculated
                            {default: all}

        state = integer: state of the (selection in the) object for
                         which principal axes should be calculated
                         {default: 1}

        scale_factor = float: the scaling factor for the cgo axes
                              representation
                              {default: 20.0}

    PYMOL API

        not implemented

    EXAMPLES

        get_principal_axes 5H3O, state=1, scale_factor=20.0

    NOTES

        To be made.
    """
    selection = str(selection)
    state = int(state)
    scale_factor = float(scale_factor)

    axes_and_cog = calculate_principal_axes_and_cog(selection, state)
    draw_principal_axes(axes_and_cog, scale_factor)


cmd.extend("get_principal_axes", get_principal_axes)
cmd.extend("align_axes", align_axes)
cmd.extend("draw_coordinate_axes", draw_coordinate_axes)