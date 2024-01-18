from pymol import cmd
from pymol import stored

"""
--- Get Colors ---
Author  : Michele Bonus
Program : Get Colors
Date    : May 2022
To Do   : Make it easier to use
"""

__author__ = "Michele Bonus"
__copyright__ = "Copyright 2022, Michele Bonus"
__credits__ = ["Michele Bonus"]
__license__ = "GPL"
__version__ = "0.5"
__maintainer__ = "Michele Bonus"
__email__ = "Michele.Bonus@hhu.de"
__status__ = "Development"

def get_rgb_color(color):
    """Get the rgb values of a named color in PyMOL
    :param color: A PyMOL color name
    :return colorlst: A list of type [R, G, B] (printed to stdout)
    """
    colorstr = str(color)
    colorlst = [
        round(color * 255) for color in cmd.get_color_tuple(cmd.get_color_index("red"))
    ]
    print(colorlst)


def get_rgb_colors_sele(selection):
    """Get the rgb values of all atoms in a selection in PyMOL
    :param selection: A PyMOL selection
    :return color_rgbs: A list of type [[R1,G1,B1],[R2,G2,B2],...,[Rn,Gn,Bn]] (printed to stdout)
    """
    stored.colors = []
    cmd.iterate(selection, "stored.colors.append(color)")
    color_tuples = [cmd.get_color_tuple(stored_color) for stored_color in stored.colors]
    color_rgbs = [
        [
            round(color_tuple[0] * 255),
            round(color_tuple[1] * 255),
            round(color_tuple[2] * 255),
        ]
        for color_tuple in color_tuples
    ]
    print(color_rgbs)


cmd.extend("get_rgb_color", get_rgb_color)
cmd.extend("get_rgb_colors_sele", get_rgb_colors_sele)
