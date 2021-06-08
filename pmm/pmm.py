'''Functions needed for MD-PMM calculations'''

from pmm.inputs import read_pmm_inputs, read_tot_input_gauss
from pmm.spectra import calc_pert_matrix, calc_uv
import sys
from timeit import default_timer as timer
from argparse import ArgumentParser, FileType
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np
from numba import njit, prange
from scipy import linalg
import pmm.conversions as conv
from pmm.conversions import Bohr2Ang


def convert2Universe(geometry: np.ndarray) -> mda.Universe:
    '''Converts a geometry expressed as a numpy.ndarry to a mda.Universe
    object.

    Parameters:
        geometry (np.ndarray): (n_atoms, 4) matrix, containing atomic masses
            and xyz coordinates. Units assumed to be in Angstrom.
    Returns:
        univ_geom (mda.Universe): xyz coordinates of the systems and atom
            types. Units in Angstrom.
    '''
    univ_geom = mda.Universe.empty(geometry.shape[0], trajectory=True)
    # print(univ_geom)
    # take the first element of each row (corresponding to the atomic masses)
    masses = [atom[0] for atom in geometry]
    univ_geom.add_TopologyAttr('mass', masses)
    # convert masses to atom type
    # print([atom[0] for atom in geometry])
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
        rot_dip_matrix (np.ndarray): rotated electric dipole moment matrix.
            Expressed in a.u..
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
    return rot_dip_matrix, rot_matrix


@njit
def calc_el_field_pot_qc(solv_coor: np.ndarray, solv_charges: np.ndarray,
                         qc_coor: np.ndarray, cdm: np.ndarray,
                         q_tot: int):
    '''Calculate the electric field on the center of mass and the energy
    contribution given by the interaction between the QC charge and the
    electric potential produced by the solvent in the framework of the QC-based
    expansion (the perturbing field).

    Parameters:
        solv_coor (np.ndarray): (n_solv_atoms, 3) array containing the xyz
            coordinates of the solvent. Units Angstrom that the function
            converts in Bohr.
        solv_charges (np.ndarray): (n_solv_atoms) array containing the force
            field charges of the solvent.
        qc_coor (np.ndarray): (n_qc_atoms, 3) array containing the xyz
            coordinates of the QC in the MD trajectory frame.
        cdm (np.ndarray): QC center of mass in the MD trajectory frame.
        q_tot (int): QC total charge.

    Returns:
        el_field (np.ndarray): electric field expressed in its xyz components.
            In a.u..
        potential (float or np.ndarray): energy contributions to the matrix
            elements of the perturbed Hamiltonian given by the interactions
            between the charge/s of the QC and the electric potential
            produced by the solvent (perturbing field). It is obtained by
            considering the QC a point-charge in its center of mass. Expressed
            in a.u..
    '''
    # converts the coordinates from the trajectory in a.u.
    xyz_distances = (cdm - solv_coor) / Bohr2Ang
    # distances = np.sqrt(np.einsum('ij,ij->i', xyz_distances, xyz_distances))
    distances = np.sqrt((xyz_distances**2).sum(axis=1))
    el_field = (((solv_charges * xyz_distances.T) /
                 (distances ** 3)).T).sum(axis=0)
    potential = q_tot * (solv_charges / distances).sum()
    # potential = np.full(qc_charges.shape[0],
    #                    q_tot * (solv_charges / distances).sum())
    return el_field, potential


@njit(parallel=True)
def calc_el_field_pot_atom(solv_coor: np.ndarray, solv_charges: np.ndarray,
                           qc_coor: np.ndarray, cdm: np.ndarray,
                           qc_charges: np.ndarray):
    '''Calculate the electric field on the center of mass and the energy
    contribution given by the interaction between the QC charge and the
    electric potential produced by the solvent (the perturbing field) in the
    framework of the atom-based expansion.

    Parameters:
        solv_coor (np.ndarray): (n_solv_atoms, 3) array containing the xyz
            coordinates of the solvent. Units Angstrom that the function
            converts in Bohr.
        solv_charges (np.ndarray): (n_solv_atoms) array containing the force
            field charges of the solvent.
        qc_coor (np.ndarray): (n_qc_atoms, 3) array containing the xyz
            coordinates of the QC in the MD trajectory frame.
        cdm (np.ndarray): QC center of mass in the MD trajectory frame.
        qc_charges (np.ndarray): arrays providing the atomic charge
            distributions of the QC.

    Returns:
        el_field (np.ndarray): electric field expressed in its xyz components.
            In a.u..
        potential (float or np.ndarray): energy contributions to the matrix
            elements of the perturbed Hamiltonian given by the interactions
            between the charge/s of the QC and the electric potential
            produced by the solvent (perturbing field). It is obtained by
            considering the charges on each atom of the QC. Expressed in a.u..
    '''
    # converts the coordinates from the trajectory in a.u.
    xyz_distances = (cdm - solv_coor) / Bohr2Ang
    # distances = np.sqrt(np.einsum('ij,ij->i', xyz_distances, xyz_distances))
    distances = np.sqrt((xyz_distances**2).sum(axis=1))
    el_field = (((solv_charges * xyz_distances.T) /
                 (distances ** 3)).T).sum(axis=0)
    potential = np.zeros(qc_charges.shape[0])
    qc_distances = np.zeros((qc_coor.shape[0], solv_coor.shape[0]))
    for i in prange(qc_coor.shape[0]):
        qc_xyz_distances = (qc_coor[i] - solv_coor) / Bohr2Ang
        # qc_distances[i, :] = np.sqrt(np.einsum('ij,ij->i', qc_xyz_distances,
        #                                 qc_xyz_distances))
        qc_distances[i, :] = np.sqrt((qc_xyz_distances**2).sum(axis=1))
        for j in prange(qc_charges.shape[0]):
            potential[j] += qc_charges[j, i] *\
                (solv_charges / qc_distances[i, :]).sum()
    return el_field, potential


def calc_pmm_matrix(energies: np.ndarray, rot_dip_matrix: np.ndarray,
                    el_field: np.ndarray,
                    potential: float, qc_ch_swith: bool) -> np.ndarray:
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
        potential (float or np.ndarray): energy contributions to the matrix
            elements of the perturbed Hamiltonian given by the interactions
            between the charge/s of the QC and the electric potential
            produced by the solvent (perturbing field). If can be obtained by
            considering the QC a point-charge in its center of mass or by
            considering the charges on each atom of the QC. Expressed in a.u..
        qc_ch_switch (bool): choose whether to use the QC-based expansion or
            the atom-based expansion.

    Returns:
        pmm_matrix (np.ndarray): Hamiltonian matrix of the perturbed system
            as calculated (... cite article).
        '''
    # TODO #1 Add reference to PMM article.
    # print('E0 + mu', np.diag(energies) - 1 * np.einsum('i,jki->jk',
    #       el_field, rot_dip_matrix))
    if qc_ch_swith:
        pmm_matrix = - 1 * np.einsum('i,jki->jk', el_field, rot_dip_matrix)
        np.fill_diagonal(pmm_matrix, energies + potential)
    else:
        pmm_matrix = np.diag(energies + potential) \
            - 1 * np.einsum('i,jki->jk', el_field, rot_dip_matrix)
    return pmm_matrix


def main():
    '''Program with CLI to perform MD-PMM calculations.
    TODO: #3 add documentation.
    '''
    parser = ArgumentParser(prog='pmm',
                            description='Perform MD-PMM calculation.')
    # Get the QC QM properties from text files.
    parser.add_argument('-g', '--ref-geom', action='store', type=str,
                        help='QC reference (QM) geometry filename')
    parser.add_argument('-gu', '--geom-units', choices=['angstrom', 'bohr',
                        'nm'], help='specify units used in the reference ' +
                        'geometry', default='angstrom')
    parser.add_argument('-dm', '--dip-matrix', action='store', type=str,
                        help='QC unperturbed electric dipole moment matrix ' +
                        'filename')
    parser.add_argument('-e', '--energies', action='store', type=str,
                        help='QC unperturbed electronic states energies ' +
                        'filename')
    parser.add_argument('-ch', '--charges', action='store', type=str,
                        default=False, help='file with the QC atomic ' +
                        'charges for each unperturbed electronic state ' +
                        '(default=False)')
    parser.add_argument('-traj', '--trajectory-path', action='store',
                        type=str, help='XTC file of the MD simulation ' +
                        'trajectory')
    parser.add_argument('-top', '--topology-path', action='store', type=str,
                        help='TPR file of the MD simulation')
    parser.add_argument('-q', '--qc-charge', action='store', default=0,
                        type=int, help='total QC charge (default=0)')
    parser.add_argument('-nm', '--mm-indexes', action='store', type=str,
                        help='indexes of the QC in the MM trajectory')
    parser.add_argument('-nq', '--qm-indexes', action='store', type=str,
                        default=False, help='indexes of the portion of the ' +
                        'reference (QM) geometry to be considered in the ' +
                        'MD-PMM calculation')
    parser.add_argument('-o', '--output', action='store', type=str,
                        default='eigenval.txt', help='perturbed QC energies ' +
                        '(default: eigenval.txt)')
    parser.add_argument('-uv', action='store_true', help='calculate UV spectra')
    cmdline = parser.parse_args()

    start = timer()
    # determine if the QC total charge is to be used of if the charge
    # distributions are provided.
    qc_ch_switch = isinstance(cmdline.charges, str)
    # print(qc_ch_switch)
    # gather the electronic properties of the QC and load the MM trajectory
    qm_inputs, mm_traj = read_pmm_inputs(cmdline)
    # print(qm_inputs['energies'])
    # geometry units: converts to MDAnalysis defaults (Angstrom).
    if cmdline.geom_units.lower() == 'bohr':
        qm_inputs['geometry'][:, 1:] *= Bohr2Ang
    elif cmdline.geom_units.lower() == 'nm':
        qm_inputs['geometry'][:, 1:] *= 10
    # print(qm_inputs['geometry'])
    # print(qm_inputs)
    qc_ref = convert2Universe(qm_inputs['geometry'])
    # cut only the portion of interest of the QC.
    if cmdline.qm_indexes:
        qc_ref = cut_qc(qc_ref, cmdline.qm_indexes)
    # shift the origin to the center of mass.
    qc_ref.atoms.positions -= qc_ref.atoms.center_of_mass()
    # MD-PMM calculation.
    eigvals = np.zeros((mm_traj.trajectory.n_frames,
                       qm_inputs['energies'].shape[0]))
    eigvecs = np.zeros((mm_traj.trajectory.n_frames,
                        qm_inputs['energies'].shape[0],
                        qm_inputs['energies'].shape[0]))
    # Divide the coordinates in the frame between QC and solvent.
    # NOTE: the indexes are inclusive of the extremes.
    qc_traj, solv_traj = split_qc_solv(mm_traj, cmdline.mm_indexes)
    # print(mm_traj.atoms.charges)
    # print(qc_traj.atoms.center_of_mass().dtype,
    #       solv_traj.atoms.positions.dtype)
    for i, frame in enumerate(mm_traj.trajectory):
        # print(solv_traj.atoms.positions)
        # cdm_qc_traj = qc_traj.atoms.center_of_mass().astype('float32')
        # print(cdm_qc_traj.dtype, solv_traj.atoms.positions.dtype)
        if qc_ch_switch:
            el_field, potential = calc_el_field_pot_atom(solv_traj.atoms.positions,
                                                         solv_traj.atoms.charges,
                                                         qc_traj.atoms.positions,
                                                         qc_traj.atoms.center_of_mass(),
                                                         qm_inputs['charges'])
        else:
            el_field, potential = calc_el_field_pot_qc(solv_traj.atoms.positions,
                                                       solv_traj.atoms.charges,
                                                       qc_traj.atoms.positions,
                                                       qc_traj.atoms.center_of_mass(),
                                                       q_tot=cmdline.qc_charge)
        rot_dip_matrix, rot_matrix = rotate_dip_matrix(qm_inputs['dip_matrix'],
                                                       qc_traj, qc_ref)
        pmm_matrix = calc_pmm_matrix(qm_inputs['energies'], rot_dip_matrix,
                                     el_field, potential, qc_ch_switch)
        # print('potential', potential, '\nel field', el_field)
        # print('H tot', pmm_matrix.diagonal())
        eigval, eigvec = linalg.eigh(pmm_matrix)
        # eigvals.append(eigval)
        # eigvecs.append(eigvec)
        eigvals[i, :] = eigval
        eigvecs[i, :, :] = eigvec
        # print('eigenvec', eigvec)
    np.savetxt(cmdline.output, eigvals,
               header='Perturbed QC energies:')
    np.save('eigenvecs', eigvecs)
    end = timer()
    print('The calculation took: ', end - start)
    # print(eigvecs)
    # print(solv_traj.atoms.positions.shape)
    if cmdline.uv:
        pert_matrix = calc_pert_matrix(qm_inputs['dip_matrix'], eigvecs)
        calc_uv(eigvals, pert_matrix)
    return 0


if __name__ == '__main__':
    main()
