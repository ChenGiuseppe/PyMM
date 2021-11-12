'''Functions needed for MD-PMM calculations'''

from pymm.inputs import read_pmm_inputs, write_geom
import itertools
import logging
import sys
from timeit import default_timer as timer
from argparse import ArgumentParser, FileType
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np
from numba import njit, prange
from scipy import linalg
import pymm.conversions as conv
from pymm.conversions import Bohr2Ang


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
    atom_types = [conv.mass2symbol[round(atom[0])] for atom in geometry]
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
        split_qc_solv(traj, "1:10 12")  # NOTE: it includes the extremes.
    '''

    new_qc_indexes = ' or bynum '.join(qc_indexes.split())
    #print(new_qc_indexes)
    qc = traj.select_atoms(f'bynum {new_qc_indexes}')
    #print(qc)
    solv = traj.select_atoms(f'not bynum {new_qc_indexes}')
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


def order_mass(traj_geom: mda.core.groups.AtomGroup,
               ref_geom: mda.Universe) -> list:
    '''Order atoms according to their mass.
        Parameters:
        traj_geom (mda.core.groups.AtomGroup): geometry of the QC in the
            considered frame of the trajectory.
        ref_geom (mda.Universe): geometry of the QC as used in the QM
            calculation. In a.u..

    Returns:
        new_indices (list): ordered indices.
    '''

    ordered = False
    indices = [i for i in range(ref_geom.atoms.n_atoms)]
    # print(ref_geom.atoms.masses)
    # print(traj_geom.atoms.masses)
    while not ordered:
        for i, j in enumerate(indices):
            # print(j, i, ref_geom.atoms.masses[j], traj_geom.atoms.masses[i])
            if round(ref_geom.atoms.masses[j]) == round(traj_geom.atoms.masses[i]):
                if i == ref_geom.atoms.n_atoms - 1:
                    ordered = True
                continue
            else:
                indices[i:] = indices[i+1:] + [indices[i]]
                break
        # print(ordered)
    return indices


def match_qc(indices: list, traj_geom: mda.core.groups.AtomGroup,
             ref_geom: mda.Universe) -> list:
    '''Find the correct order of the QC atoms to match the MD simulation.
    A brute force approach was employed (obtain all the permutations and
    select the one with the smallest value of RMSD).

    Parameters:
        indices (list): list of indices ordered in order to already match
            the masses of the QC in the MD simulation.
        traj_geom (mda.core.groups.AtomGroup): geometry of the QC in the
            considered frame of the trajectory.
        ref_geom (mda.Universe): geometry of the QC as used in the QM
            calculation. In a.u..

    Returns:
        new_indices: list of indices that match the QC atoms order in the
        MD simulation.
    '''

    cdm = traj_geom.atoms.center_of_mass()
    coor = traj_geom.atoms.positions.copy()
    coor -= cdm
    # cdm2 = np.zeros(3)
    # for i in range(coor.shape[0]):
    #     cdm2 += traj_geom.atoms.masses[i]*(coor[i,:])
    # print('cdm', cdm2)

    elements = set([round(i) for i in ref_geom.atoms.masses])

    # Not ordered masses
    masses = ref_geom.atoms.masses.copy()
    # print('indices', indices)
    tot_ndx = []
    for element in elements:
        element_ndx = []
        for i, j in enumerate(indices):
            if round(masses[j]) == element:
                element_ndx.append([i,j])
        # print(element, element_ndx)
        perms = list(itertools.permutations([i[1] for i in element_ndx]))
        rmsd_min = 100000.
        new_element_ndx = []
        # print(element, [i[0] for i in element_ndx])

        if not len([i[0] for i in element_ndx]) == 1:
            for perm in perms:
                #print('perm', list(perm))
                perm = list(perm)
                geom_tmp = ref_geom.atoms.positions[perm,:]
                masses_tmp = ref_geom.atoms.masses[perm]
                #print(masses_tmp)
                #print([i[0] for i in perm])
                matrix, rmsd_tmp = align.rotation_matrix(coor[[i[0] for i in element_ndx],:],
                                                         geom_tmp, weights=masses_tmp)
                if rmsd_tmp < rmsd_min:
                    rmsd_min = rmsd_tmp
                    new_element_ndx = perm
        else:
            new_element_ndx = [i[1] for i in element_ndx]
        #print(new_element_ndx)
        tot_ndx += [[element_ndx[i][0], new_element_ndx[i]] for i in range(len(element_ndx))]

    new_indices = indices[:]
    for i in tot_ndx:
        #print(i)
        new_indices[i[0]] = i[1]

    logging.info(' ** Old indices:\n' +' '.join([str(i) for i in indices]) + '\n' +
                 ' ** Ordered indices:\n' + ' '.join([str(i) for i  in new_indices]) +
                 '\n')

    return new_indices



def rotate_dip_matrix(dip_matrix: np.ndarray,
                      traj_geom: mda.core.groups.AtomGroup,
                      ref_geom: mda.Universe) -> np.ndarray:
    '''Rotate the electric dipole moment matrix in order to align the
    geometry in the simulation trajectory frame to the reference geometry
    (that is the one used in the QM calculation).

    Parameters:
        dip_matrix (np.ndarray): matrix of the electric dipole moments.
            In a.u..
        traj_geom (mda.core.groups.AtomGroup): geometry of the QC in the
            considered frame of the trajectory.
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


def pmm(cmdline):
    '''Program to perform MD-PMM calculations.
    TODO: #3 add documentation.
    '''

    # determine if the QC total charge is to be used or if the charge
    # distributions are provided.
    qc_ch_switch = isinstance(cmdline.charges, str)
    if qc_ch_switch:
        logging.info(' * PMM-MD calculation using the charges approximation ' +
                     'has been selected.\n')
    else:
        logging.info(' * PMM-MD calculation using the dipole approximation ' +
                     'has been selected.\n')
    # print(qc_ch_switch)
    # gather the electronic properties of the QC and load the MM trajectory
    qc, mm_traj = read_pmm_inputs(cmdline)
    # print(qm_inputs['energies'])
    # geometry units: converts to MDAnalysis defaults (Angstrom).
    if cmdline.geom_units.lower() == 'bohr':
        qc.geom[:, 1:] *= Bohr2Ang
    elif cmdline.geom_units.lower() == 'nm':
        qc.geom[:, 1:] *= 10

    # print(qm_inputs['geometry'])
    # print(qm_inputs)
    qc_ref = convert2Universe(qc.geom)
    # cut only the portion of interest of the QC.
    if cmdline.qm_indexes:
        qc_ref = cut_qc(qc_ref, cmdline.qm_indexes)
    # shift the origin to the center of mass.
    qc_ref.atoms.positions -= qc_ref.atoms.center_of_mass()

    # Divide the coordinates in the frame between QC and solvent.
    # NOTE: the indexes are inclusive of the extremes.
    qc_traj, solv_traj = split_qc_solv(mm_traj, cmdline.mm_indexes)

    n_qc_atoms = qc.geom.shape[0]
    logging.info('\n\n'
                 '==========================================================================\n'
                 '|                                                                        |\n'
                 '|                Properties of the unperturbed QC:                       |\n'
                 '|                                                                        |\n'
                 '==========================================================================\n'
                 '\n')
    logging.info(' * Initial QC geometry (a.u.):\nNumber of atoms: {}'.format(n_qc_atoms))
    for i in range(n_qc_atoms):
        logging.info('{:7.4} {:12.7f} {:12.7f} {:12.7f}'.format(qc.geom[i, 0],
                                          *qc.geom[i, 1:] / Bohr2Ang))
    logging.info('--------------------------------------------------------------------------\n')

    if cmdline.match:
        logging.info(' * The reference QC geometry will be reordered to match '
                     'the atoms order in the MD simulation.\n')
        # match QC atoms order of the reference geometry with the MD simulation. 
        new_indices = match_qc(order_mass(qc_traj, qc_ref), qc_traj, qc_ref)
        qc_ref.atoms.masses = qc_ref.atoms.masses[new_indices]
        qc_ref.atoms.positions = qc_ref.atoms.positions[new_indices,:]
        qc.geom = qc.geom[new_indices,:]
        if qc_ch_switch:
            qc.charges = qc.charges[:,new_indices]

        logging.info('--------------------------------------------------------------------------\n')
        logging.info(' * Final QC geometry (a.u.):\nNumber of atoms: {}'.format(n_qc_atoms))
        for i in range(n_qc_atoms):
            logging.info('{:7.4} {:12.7f} {:12.7f} {:12.7f}'.format(qc.geom[i, 0],
                                                                    *qc.geom[i, 1:] / Bohr2Ang))
        logging.info('\n--------------------------------------------------------------------------\n')

    n_qc_states = qc.energies.shape[0]
    logging.info(' * QC electronic states energies (a.u.):\n' +
                 'Number of electronic states: {}'.format(n_qc_states))
    for i in range(n_qc_states):
        logging.info('{:5d} {:12.7f}'.format(i, qc.energies[i]))
    logging.info('\n--------------------------------------------------------------------------\n')

    logging.info(' * Electric dipole matrix (a.u.):')
    for i in range(n_qc_states):
        for j in range(n_qc_states):
            logging.info('{:5d} {:5d} {:12.7f} {:12.7f} {:12.7f}'.format(i, j, *qc.dip_matrix[i,j,:]))
    logging.info('\n--------------------------------------------------------------------------\n')

    if qc_ch_switch:
        logging.info(' * QC atomic charges:\n\n')
        for i in range(n_qc_states):
            logging.info(' ** Electronic state {}:'.format(i))
            for j in range(n_qc_atoms):
                logging.info('{:7.4f} {:12.7f}'.format(qc_ref.atoms.masses[j], qc.charges[i,j]))
            logging.info('\n')
        logging.info('\n--------------------------------------------------------------------------\n')

    logging.info('\n\n'
                 '==========================================================================\n'
                 '|                                                                        |\n'
                 '|                       MD simulation details                            |\n'
                 '|                                                                        |\n'
                 '==========================================================================\n')
    logging.info(' * Total number of atoms: {}\n'.format(mm_traj.atoms.n_atoms) +
                 ' * Number of frames in the trajectory: {}\n'.format(mm_traj.trajectory.n_frames) +
                 ' * QC atoms indeces in the total system (starting from 1):\n' +
                 ' '.join([str(i + 1) for i in qc_ref.atoms.indices]) +
                 '\n\n--------------------------------------------------------------------------\n')

    # MD-PMM calculation.
    eigvals = np.zeros((mm_traj.trajectory.n_frames,
                       qc.energies.shape[0]))
    eigvecs = np.zeros((mm_traj.trajectory.n_frames,
                        qc.energies.shape[0],
                        qc.energies.shape[0]))

    # print(mm_traj.atoms.charges)
    # print(qc_traj.atoms.center_of_mass().dtype,
    #       solv_traj.atoms.positions.dtype)
    for i, frame in enumerate(mm_traj.trajectory):
        # print(i)
        # print(solv_traj.atoms.positions)
        # cdm_qc_traj = qc_traj.atoms.center_of_mass().astype('float32')
        # print(cdm_qc_traj.dtype, solv_traj.atoms.positions.dtype)
        if qc_ch_switch:
            el_field, potential = calc_el_field_pot_atom(solv_traj.atoms.positions,
                                                         solv_traj.atoms.charges,
                                                         qc_traj.atoms.positions,
                                                         qc_traj.atoms.center_of_mass(),
                                                         qc.charges)
        else:
            el_field, potential = calc_el_field_pot_qc(solv_traj.atoms.positions,
                                                       solv_traj.atoms.charges,
                                                       qc_traj.atoms.positions,
                                                       qc_traj.atoms.center_of_mass(),
                                                       q_tot=cmdline.qc_charge)
        rot_dip_matrix, rot_matrix = rotate_dip_matrix(qc.dip_matrix,
                                                       qc_traj, qc_ref)
        pmm_matrix = calc_pmm_matrix(qc.energies, rot_dip_matrix,
                                     el_field, potential, qc_ch_switch)
        # print('potential', potential, '\nel field', el_field)
        # print('H tot', pmm_matrix.diagonal())
        eigval, eigvec = linalg.eigh(pmm_matrix)
        # eigvals.append(eigval)
        # eigvecs.append(eigvec)
        eigvals[i, :] = eigval
        eigvecs[i, :, :] = eigvec
        # print('eigenvec', eigvec)

    np.savetxt('{}_eigvals.dat'.format(cmdline.output), eigvals,
               header='Perturbed QC energies:')
    np.save('{}_eigvecs'.format(cmdline.output), eigvecs)

    xyz_fn = '{}_qc_geom.xyz'.format(cmdline.output)
    write_geom(xyz_fn, qc.geom)
    logging.info('==========================================================================\n'
                 '|                                                                        |\n'
                 '|                               OUTPUT                                   |\n'
                 '|                                                                        |\n'
                 '==========================================================================\n'
                 '\n'
                 ' * Results have been saved in:\n' +
                 ' ** Perturbed energies: {}_eigvals.dat\n'.format(cmdline.output) +
                 ' ** Eigenvectors: {}_eigvecs.npy\n'.format(cmdline.output) +
                 ' ** QC reference geometry in .xyz format: {}'.format(xyz_fn))
    # print(eigvecs)
    # print(solv_traj.atoms.positions.shape)

    return 0


if __name__ == '__main__':
    pass
