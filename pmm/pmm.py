'''Functions needed for MD-PMM calculations'''

import MDAnalysis as mda
import numpy as np
import conversions as conv


def convert2Universe(geometry: np.ndarray) -> mda.Universe:
    '''Converts a geometry expressed as a numpy.ndarry to a mda.Universe
    object.

    Parameters:
        geometry (np.ndarray): (n_atoms, 4) matrix, containing atomic masses
            and xyz coordinates.
    Returns:
        univ_geom (mda.Universe): xyz coordinates of the systems and atom
            types.
    '''
    univ_geom = mda.Universe.empty(geometry.shape[0], trajectory=True)
    #print(univ_geom)
    # take the first element of each row (corresponding to the atomic masses)
    masses = [atom[0] for atom in geometry]
    univ_geom.add_TopologyAttr('mass', masses)
    # convert masses to atom type
    atom_types = [conv.mass2symbol[atom[0]] for atom in geometry]
    univ_geom.add_TopologyAttr('type', atom_types)
    # add xyz coordinates
    univ_geom.atoms.positions = geometry[:,1:]
    #print(univ_geom.atoms.positions, univ_geom.atoms.masses)
    return univ_geom


def pmm_matrix(energies, rot_dip_matrix, el_field):
    pass

if __name__ == '__main__':
    h2o = np.array([[16.0, 0,        0,       0      ],  # oxygen
                [1.0, 0.95908, -0.02691, 0.03231],  # hydrogen
                [1.0, -0.28004, -0.58767, 0.70556]])
    convert2Universe(h2o)