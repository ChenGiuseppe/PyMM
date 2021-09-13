import MDAnalysis as mda
import numpy as np
import pmm.conversions as conv


class Qc:
    def __init__(self, geom, dip_matrix, energies, charges=None):
        self.geom = geom
        self.dip_matrix = dip_matrix
        self.energies = energies
        self.charges = charges

def read_raw_geom(geom_fn: str) -> np.ndarray:
    '''Read the geometry directly from a text file
    (for each X atom: "X x y z\n").

    Parameters:
        geom_fn (str): geometry filename.

    Returns:
        geom_tmp (np.ndarray): geometry as (n_atoms, 4) numpy array.
    '''
    with open(geom_fn, 'r') as geom_file:
        geom = []
        for i, line in enumerate(geom_file):
            content = line.split()
            # Check if the first line is used for specifying the number of
            # atoms or something else.
            if i == 0 and not len(content) == 4:
                continue
            # NOTE: skip empty line. It was added having in mind the
            # possibility that the last line could be empty, as it is now it
            # accounts for all empty lines.
            elif content == []:
                continue
            elif not len(content) == 4:
                raise ValueError("There's an error in the geometry file " +
                                 f'at line {i}: there are ' +
                                 f'{len(content)} elements.\n' +
                                 f'{line}')
            try:
                geom.append([float(j) for j in content])
                continue
            except ValueError:
                pass
            try:
                # Try to convert the first element of the list (identifying
                # the element) to its corresponding atomic mass (float).
                # It takes only the first character, expecting something like:
                # Examples: 'H', 'H1', 'HA', etc. and not '1H', 'XH', etc..
                geom.append([conv.symbol2mass[content[0][0]]] +
                            [float(j) for j in content[1:]])
            except ValueError:
                print('Expects the fist element of each line',
                      'to begin with the atomic mass (float)',
                      'or a string that begins with the letter',
                      'corresponding to its symbol.')
    geom_tmp = np.zeros_like(geom)
    geom_tmp[:] = geom
    return geom_tmp


def read_raw_matrix(matrix_fn: str) -> np.ndarray:
    '''Read a (N, N, 3) matrix directly from a text file
    (for each matrix element i j: "i j x y z\n").

    Parameters:
        matrix_fn (str): matrix filename.

    Returns:
        matrix (np.ndarray): matrix as a numpy array.
    '''
    with open(matrix_fn, 'r') as matrix_file:
        matrix_x = []
        matrix_y = []
        matrix_z = []
        for i, line in enumerate(matrix_file):
            content = line.split()
            # print(i, content)
            # Check if the first line is to be read.
            if i == 0 and not (len(content) == 3 or len(content) == 5):
                continue
            # NOTE: skip empty line. It was added having in mind the
            # possibility that the last line could be empty, as it is now it
            # accounts for all empty lines.
            elif content == []:
                continue
            elif not (len(content) == 3 or len(content) == 5):
                raise ValueError("There's an error in the matrix file " +
                                 f'at line {i}: there are ' +
                                 f'{len(content)} elements.\n' +
                                 f'{line}')
            if len(content) == 3:
                matrix_x.append(float(content[0]))
                matrix_y.append(float(content[1]))
                matrix_z.append(float(content[2]))
            elif len(content) == 5:
                matrix_x.append(float(content[2]))
                matrix_y.append(float(content[3]))
                matrix_z.append(float(content[4]))
    n_rows = round(np.sqrt(len(matrix_x)))
    # Workaround to the maximum dimension of npdarray (32).
    matrix_tmp = np.zeros_like(matrix_x)
    # print(len(matrix_x), np.sqrt(len(matrix_x)), n_rows)
    matrix_tmp[:] = matrix_x
    matrix_x = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
    matrix_tmp[:] = matrix_y
    matrix_y = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
    matrix_tmp[:] = matrix_z
    matrix_z = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
    # print(matrix_x, matrix_y, matrix_z)
    matrix = np.stack((matrix_x, matrix_y, matrix_z), axis=2)
    return matrix


def read_raw_energies(energies_fn: str) -> np.ndarray:
    '''Read the electronic state energies from a text file
    (for each ith electronic state: "Ei\n").

    Parameters:
        energies_fn (str): energies filename.

    Returns:
        energies_tmp (np.ndarray): energies as a numpy array.
    '''
    with open(energies_fn, 'r') as energies_file:
        energies = []
        for line in energies_file:
            if line.split() == []:
                continue
            try:
                energies.append(float(line.strip()))
            except ValueError:
                print('There are values that are not float',
                      'in the energies file.')
    energies_tmp = np.zeros_like(energies)
    energies_tmp[:] = energies
    return energies_tmp


def read_raw_charges(charges_fn: str) -> np.ndarray:
    '''Read the RESP charges from a text file
    (a line for each electronic state, containing the charges for each atom).

    Parameters:
        charges_fn (str): charges filename.

    Returns:
        charges_tmp (np.ndarray): charges as a numpy matrix file
        (n_el_states, n_atoms).
    '''
    charges = []
    with open(charges_fn, 'r') as charges_file:
        for line in charges_file:
            content = line.split()
            if content == []:
                continue
            try:
                charges.append([float(i) for i in content])
            except ValueError:
                print('There are values that are not float in the',
                      'charges file.')
    charges_tmp = np.zeros_like(charges)
    charges_tmp[:, :] = charges
    return charges_tmp


def read_raw_inputs(geom_fn: str, dip_mat_fn: str,
                   energies_fn: str, charges_fn=False) -> None:
    if charges_fn:
        qm_charges = read_raw_charges(charges_fn)
    else:
        qm_charges = None
    qc = Qc(read_raw_geom(geom_fn),
            read_raw_matrix(dip_mat_fn),
            read_raw_energies(energies_fn),
            qm_charges)
    return qc

def read_pmm_inputs(cmdline):
    '''Obtain the pmm inputs according to the different sources.

    Parameters:
        cmdline (argparse.Namespace): input given through the command line.
            If the legacy format is used it consists of the input filename,
            otherwise the direct parsing of the output files from other
            software packages is requested using the arguments given in the
            command line.

    Returns:
        qm_inputs (dict): "geometry": geometry in Angstrom (numpy.darray,
                shape=(n_atoms, 4)).
            "energies": electronic states energies in a.u. (numpy.darray,
                shape=n_el_states).
            "dip_matrix": electric dipole moment matrix in a.u. (numpy.darray,
                shape=(n_el_states, n_el_states, 3)). INCOMPLETE: values
                only for the (0,i), (i,0) and the (el_state, el_state)
                elements.
            "charges": RESP charges for the el_state electronic state
                (numpy.darray, shape=(n_el_states, n_atoms)).
        mm_traj (mda.Universe): MM simulation trajectory.
    '''
    qc = read_raw_inputs(cmdline.ref_geom, cmdline.dip_matrix,
                               cmdline.energies, cmdline.charges)
    if '.tpr' in cmdline.topology_path or '.xtc' in cmdline.trajectory_path:
            mm_traj = mda.Universe(cmdline.topology_path,
                                   cmdline.trajectory_path)
    return qc, mm_traj

if __name__ == '__main__':
    pass