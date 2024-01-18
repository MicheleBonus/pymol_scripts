"""
--- Disulfide Bond Information Tool ---
Author  : Michele Bonus
Program : Disulfide Bond Information Tool
Date    : 17.01.2024
Purpose : Extract and print information about disulfide bonds in protein structures
          in LEaP bond command format using PyMOL.
"""

__author__ = "Michele Bonus"
__copyright__ = "Copyright 2024, Michele Bonus"
__credits__ = ["Michele Bonus"]
__license__ = "MIT"
__version__ = "0.3"
__maintainer__ = "Michele Bonus"
__email__ = "Michele.Bonus@hhu.de"
__status__ = "Development"

from pymol import cmd, stored

# Define constants
SELECTION_EVERYTHING = "everything"
RESIDUE_CYS_CYX = "CYS+CYX"
ELEMENT_SULFUR = "S"
SELECTION_NONE = "none"
SELECTION_SKIP = "skip"
SELECTION_S1 = "s1"
SELECTION_S2 = "s2"

def get_disulfide_info_for_leap(selection="everything"):
    """
    This function identifies pairs of cysteine SG atoms that form disulfide bonds in a protein
    structure, and prints the bond information in the format used by the LEaP tool from the AMBER
    molecular dynamics package.

    The function uses PyMOL's selection and iteration capabilities to find pairs of sulfur (S) atoms
    from cysteine (CYS/CYX) residues that are bound to each other, indicating a disulfide bond.
    For each pair found, it prints a line in the format

    "bond X.<residue number>.SG X.<residue number>.SG",

    which can be modified and used in a LEaP input script to specify the disulfide bond.

    Note that the script cannot take into account the internal renumbering of residues by LEaP.
    Although LEaP should now use the residue numbering of the input file, it will still renumber
    residues internally if they have an insertion code, for example, or if the numbering begins anew
    after a chain break. Ideally, you should therefore use a PDB file without a chain identifier
    and with consecutive residue numbering as input.

    Args:
        selection (str, optional): A PyMOL selection string specifying the part of the molecule to
                                   search for disulfide bonds.
                                   The default is "everything", which means the entire molecule
                                   will be searched.

    Raises:
        PyMOLCmdException: If an error occurs during the PyMOL commands
                           (e.g., the selection string is invalid).

    Example:
        To find and print disulfide bonds in the entire molecule:
        >>> get_disulfide_info_for_leap all

        To find and print disulfide bonds in chain A only:
        >>> get_disulfide_info_for_leap chain A
    """
    cmd.select(SELECTION_SKIP, SELECTION_NONE)
    disulfide_cys_sg_indices = cmd.index(
        f"{selection} and ((resn {RESIDUE_CYS_CYX} and elem {ELEMENT_SULFUR}) and"
        f" bound_to (resn {RESIDUE_CYS_CYX} and elem {ELEMENT_SULFUR}))"
    )

    for _, sg_atom in enumerate(disulfide_cys_sg_indices):
        selection_string_sg = f"({sg_atom[0]} and index {sg_atom[1]})"

        if cmd.select(SELECTION_S1, selection_string_sg) and cmd.select(
            SELECTION_S2,
            f"(elem {ELEMENT_SULFUR}) and (bound_to {selection_string_sg}) and not"
            f" ?{SELECTION_SKIP}",
        ):
            stored.info = []
            cmd.iterate(
                f"{SELECTION_S1} or {SELECTION_S2}",
                "stored.info.append((chain, resn, resi, name))",
            )
            atom1, atom2 = stored.info
            print(f"bond X.{atom1[2]}.SG X.{atom2[2]}.SG")
            cmd.select(SELECTION_SKIP, f"{SELECTION_S1} or {SELECTION_S2} or ?{SELECTION_SKIP}")

    cmd.delete(SELECTION_SKIP)
    cmd.delete(SELECTION_S1)
    cmd.delete(SELECTION_S2)

cmd.extend("get_disulfide_info_for_leap", get_disulfide_info_for_leap)
