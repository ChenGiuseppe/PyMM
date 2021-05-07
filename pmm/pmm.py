'''Functions needed for MD-PMM calculations'''

from conversions import Bohr2Ang
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
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
    # print(univ_geom)
    # take the first element of each row (corresponding to the atomic masses)
    masses = [atom[0] for atom in geometry]
    univ_geom.add_TopologyAttr('mass', masses)
    # convert masses to atom type
    #print([atom[0] for atom in geometry])
    atom_types = [conv.mass2symbol[atom[0]] for atom in geometry]
    univ_geom.add_TopologyAttr('type', atom_types)
    # add xyz coordinates
    univ_geom.atoms.positions = geometry[:, 1:]
    # print(univ_geom.atoms.positions, univ_geom.atoms.masses)
    return univ_geom


def split_qc_solv(traj: mda.Universe,
                  qc_indexes: str) -> mda.core.groups.AtomGroup:
    '''Select AtomGroups for the QC (quantum center) and the "solvent"
    (i.e. the perturbing field).

    Parameters:
        traj (mda.Universe): trajectory sampled by the MD simulation.
        qc_indexes (str): string expressing the selection of the QC atoms in
            the formalism used in MDAnalysis (bynum, so starting from 1).

    Returns:
        qc (mda.core.groups.AtomGroup): AtomGroup containing the QC.
        solv (mda.core.groups.AtomGroup): AtomGroup containing the solvent
           (complementary to the QC).

    Example:
        split_qc_solv(traj, "1:10 or 12")  # NOTE: it includes the extremes.
    '''
    qc = traj.select_atoms(f'bynum {qc_indexes}')
    solv = traj.select_atoms(f'not bynum {qc_indexes}')
    return qc, solv


def cut_qc(qc: mda.Universe, qc_indexes: str) -> mda.core.groups.AtomGroup:
    '''Select only a portion of the system used for the QM calculation.

    Parameters:
        qc (mda.Universe): QC as it was considered for the QM calculation.
        qc_indexes (str): string expressing the selection of the portion of
            the system used in the QM calculation to be used in the MD-PMM
            calculation, in the formalism used in MDAnalysis (bynum, so
            starting from 1).

    Returns:
        qc_pmm (mda.core.groups.AtomGroup): AtomGroup of the portion of the
            system used in the QM calculation to be used in the MD-PMM
            calculation.
    '''
    qc_pmm = qc.select_atoms(f'bynum {qc_indexes}')
    return qc_pmm


def rotate_dip_matrix(dip_matrix: np.ndarray,
                      traj_geom: mda.core.groups.AtomGroup,
                      ref_geom: mda.Universe) -> np.ndarray:
    '''Rotate the electric dipole moment matrix in order to align the
    geometry in the simulation trajectory frame to the reference geometry
    (that is the one used in the QM calculation).

    Parameters:
        dip_matrix (np.ndarray): matrix of the electric dipole moments.
            In a.u..
        frame_geom (mda.core.groups.AtomGroup): geometry of the QC in the
            considered frame of the trajectory. In nm if the simulation is
            done in Gromacs.
        ref_geom (mda.Universe): geometry of the QC as used in the QM
            calculation. In a.u..

    Returns:
        rot_dip_matrix (np.ndarray): rotated electric dipole moment matrix. In a.u..
    '''
    # shift origin to the centers of mass of frame_geom to (0, 0, 0).
    cdm = traj_geom.atoms.center_of_mass()
    # NOTE: it changes the positions in traj_geom permanently.
    traj_geom.atoms.positions = traj_geom.atoms.positions - cdm
    rot_matrix, rmsd = align.rotation_matrix(traj_geom.atoms.positions,
                                             ref_geom.atoms.positions,
                                             weights=ref_geom.atoms.masses)
    # rot_matrix is transposed to obtain the inverse. This way dip_matrix is
    # rotated into the Gromacs reference system.
    rot_dip_matrix = np.einsum('ij,klj->kli', rot_matrix.T, dip_matrix)
    return rot_dip_matrix


def calc_el_field_pot(solv_coor: np.ndarray, charges: np.ndarray,
                      ref_origin: np.ndarray) -> tuple[np.ndarray, float]:
    '''Calculate the electric field and potential due to the solvent
    considered as a distribution of point charges sampled by the simulation
    on a specified point (in the PMM: center of mass of the Quantum Center).

    Parameters:
        solv_coor (np.ndarray): (n_atoms, 3) array containing the xyz
            coordinates of the solvent. Units in nm that the function converts
            in Bohr.
        charges (np.ndarray): (n_atoms) array containing the force field
            charges of the solvent.
        ref_origin (np.ndarray): point on which the electric field is applied
            (in the PMM: center of mass the Quantum Center). Units in Bohr.

    Returns:
        el_field (np.ndarray): electric field expressed in its xyz components.
            In a.u..
        potential (float): electric potential. In a.u..
    '''
    # converts the coordinates from the trajectory in a.u.
    xyz_distances = ref_origin - solv_coor * 10 / Bohr2Ang
    distances = np.sqrt(np.einsum('ij,ij->i', xyz_distances, xyz_distances))
    el_field = (((charges * xyz_distances.T) / distances ** 3).T).sum(axis=0)
    potential = (charges / distances).sum()
    return el_field, potential


def pmm_matrix(energies: np.ndarray, rot_dip_matrix: np.ndarray,
               el_field: np.ndarray, potential: float, qc_qtot=0) -> np.ndarray:
    '''Construct PMM matrix.

    Parameters:
        energies (np.ndarray): energies of the electronic states (eigenvalues
            of the unperturbed Hamiltonian) expressed in Hartree (a.u.).
        rot_dip_matrix (np.ndarray): electric dipole matrix rotated in order
            to align the reference geometry to the geometry in the simulation
            trajectory frame. Values expressed in a.u..
        el_field (np.ndarray): electric field generated by the solvent
            considered as a distribution of point charges sampled by the
            simulation. Expressed in a.u..
        potential (float): electric potential generated by the solvent
            considered as a distribution of point charges sampled by the
            simulation. Expressed in a.u..

    Returns:
        pmm_matrix (np.ndarray): Hamiltonian matrix of the perturbed system
            as calculated (... cite article).
        '''
    # TODO #1 Add reference to PMM article.
    # diagonal elements
    pmm_matrix = np.diag(energies + qc_qtot*potential) + -1*np.einsum('i,jki->jk', el_field, rot_dip_matrix)
    return pmm_matrix



if __name__ == '__main__':
    h2o = np.array([[16.0, 0,        0,       0],  # oxygen
                   [1.0, 0.95908, -0.02691, 0.03231],  # hydrogen
                   [1.0, -0.28004, -0.58767, 0.70556]])
    convert2Universe(h2o)
