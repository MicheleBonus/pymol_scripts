from itertools import groupby
from operator import itemgetter

import numpy as np
import pymol
from pymol import cmd, selector, stored


def get_range(selection = "all", state = 1):
    selection = str(selection)
    state = int(state)
    
    stored.resi_of_atoms = []
    cmd.iterate_state(state, selection, "stored.resi_of_atoms.append(resi)")
    resi_of_atoms_arr = np.array(stored.resi_of_atoms).astype(int)
    resi_unique = np.unique(resi_of_atoms_arr)
    consecutive_groups = []
    for k, g in groupby(enumerate(resi_unique), lambda ix: ix[0] - ix[1]):
        consecutive_groups.append(list(map(itemgetter(1), g)))
    for consecutive_group in consecutive_groups:
        if len(consecutive_group) > 1:
            print(str(consecutive_group[0]) + "-" + str(consecutive_group[-1]))
        else:
            print(str(consecutive_group[0]))

cmd.extend("get_range", get_range)
