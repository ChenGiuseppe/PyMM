
import MDAnalysis as mda
import numpy as np
import pmm.conversions as conv


def get_raw_geom(geom_fn: str) -> np.ndarray:
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


def get_raw_matrix(matrix_fn: str) -> np.ndarray:
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


def get_raw_energies(energies_fn: str) -> np.ndarray:
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


def get_raw_charges(charges_fn: str) -> np.ndarray:
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


def get_raw_inputs(geom_fn: str, dip_mat_fn: str,
                   energies_fn: str, charges_fn=False) -> None:
    qm_inputs = {}
    qm_inputs['geometry'] = get_raw_geom(geom_fn)
    qm_inputs['dip_matrix'] = get_raw_matrix(dip_mat_fn)
    qm_inputs['energies'] = get_raw_energies(energies_fn)
    if charges_fn:
        qm_inputs['charges'] = get_raw_charges(charges_fn)
    else:
        qm_inputs['charges'] = None
    return qm_inputs


def legacy_input(filename: str, cmdline) -> dict:
    '''Provide the inputs necessary to the MD-PMM calculation using the legacy
    format:
    1.  geometry filename
    2.  dipole matrix filename
    3.  roots
    4.  solvent indexes filename
    5.  QC indexes filename
    6.  topology with atom charges filename
    7.  number of total atoms
    8.  energies filename
    9.  gromacs .xtc file
    10. first frame
    11. last frame
    12. output roots
    13. starting state ?
    14. arriving state ?
    15. starting state ?
    16. arriving state ?
    17. output filename
    18. total QC charge
    19. RESP charges
    20. diagok
    ---------------------
    21. Gromacs .tpr file NOTE: this is added simply for ease of use.
    '''
    # NOTE: still to be finished and tested.
    qm_inputs = {}
    with open(filename, 'rt') as input_file:
        input_file.readlines()
        geom_fn = input_file[0].strip()
        dip_matrix_fn = input_file[1].strip()
        cmdline.roots = int(input_file[2].strip)
        qc_ndx_fn = input_file[4].strip()
        # top_fn = input_file[5].strip()
        energies_fn = input_file[7].strip()
        traj_fn = input_file[8].strip()
        cmdline.qc_charge = input_file[17].strip()
        charges_fn = input_file[18].strip()
        tpr_fn = input_file[20].strip()
    # QC geometry
    with open(geom_fn, 'r') as geom_file:
        n_qc_atoms = int(next(geom_file).strip())
        geom = np.zeros((n_qc_atoms, 4))
        for i, line in enumerate(geom_file):
            geom[i, :] = [float(j) for j in line.split()]
        qm_inputs['geometry'] = geom
    # QC dip_matrix
    with open(dip_matrix_fn, 'r') as dip_matrix_file:
        dip_matrix = np.zeros((cmdline.roots, cmdline.roots, 3))
        for i in range(cmdline.roots):
            for j in range(cmdline.roots):
                dip_matrix[i, j, :] =\
                    [float(k) for k in next(dip_matrix_file).split()[2:]]
        qm_inputs['dip_matrix'] = dip_matrix
    # QC energies
    with open(energies_fn, 'r') as energies_file:
        energies = np.zeros(cmdline.roots)
        for i, line in enumerate(energies_file):
            energies[i] = float(line)
        qm_inputs['energies'] = energies
    # QC charges
    with open(charges_fn, 'r') as charges_file:
        for state in range(cmdline.roots):
            qm_inputs[f'charges_{state}'] = np.array([float(i) for i in
                                                     next(charges_file).split()
                                                      ])
    # MD system charges
    mm_traj = mda.Universe(tpr_fn, traj_fn)
    return qm_inputs, mm_traj


def get_input_gauss(filename: str, n_el_states: int, el_state: int,
                    get_geom=True, get_energies=True, get_trans_dip=True,
                    get_diag_dip=True, get_charges=False) -> dict:
    '''Extract the geometry, the electronic states energies, the elements
    of the electric dipole moment matrix and/or the RESP charges from the
    Gaussian (specifically 16, support for other versions is to be verified)
    software calculation outputs.
    The geometry extracted is in the "standard orientation" if this can be
    found, while if the Gaussian calculation was performed using the keyword
    NoSymmetry, "Input orientation" is selected instead.
    The permanent dipole moment and the RESP charges extracted refer to
    the only electronic state considered in the Gaussian calculaton
    (ground state or excited state).
    All the returned properties are expressed in atomic units.

    Parameters:
        filename (str): A string indicating the name of the Gaussian
            calculation output file.
        n_el_states (int): A integer referring to the total number of
            electronic states considered (including the ground state).
        el_state (int): A integer indicating the electronic state considered
            in the calculation (starting from 0 for the ground state).
        get_geom (bool): A boolean that specifies if extracting the geometry
            is wanted.
        get_energies (bool): A boolean that specifies if extracting the
            electronic states energies is wanted.
        get_trans_dip (bool): A boolean that specifies if extracting the
            not diagonal elements of the electric dipole moment matrix
            is wanted.
        get_diag_dip (bool): A boolean that specifies if extracting the
            diagonal elements of the electric dipole moment matrix is wanted.
        get_charges (bool): A boolean that specifies if extracting the RESP
            charges is desired.

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
    '''
    qm_inputs = {}
    with open(filename, 'r') as gout:
        geom = []
        # geometry using the input orientation instead of the standard one.
        geom_no_symm = []
        energies = np.zeros(n_el_states)
        dip_matrix = np.zeros((n_el_states, n_el_states, 3))
        charges = []
        for line in gout:
            if 'Standard orientation:' in line and get_geom:
                for i in range(4):
                    next(gout)
                content = next(gout).split()
                while '---' not in content[0]:
                    geom.append([conv.Z2mass[content[1]]] +
                                [float(coor) for coor in content[3:6]])
                    content = next(gout).split()
            elif 'Input orientation:' in line and get_geom:
                for i in range(4):
                    next(gout)
                content = next(gout).split()
                while '---' not in content[0]:
                    geom_no_symm.append([conv.Z2mass[content[1]]] +
                                        [float(coor) for coor in content[3:6]])
                    content = next(gout).split()
            elif 'SCF Done:  E(' in line and get_energies:
                energies[0] = float(line.split()[4])
            elif 'Excited State' in line and get_energies:
                # check if this condition is sufficient to cover all the
                # possible cases we take the excitation energy expressed
                # in eV, which is then converted to a.u.
                content = line.split()
                energies[int(content[2][:-1])] =\
                    energies[0] + float(content[4])/conv.au2eV
            # 'Ground to excited state transition electric
            # dipole moments (Au):' is too long. In order to abide to PEP8
            # guidelines a shorter string is used.
            elif 'transition electric dipole moments (Au):'\
                    in line and get_trans_dip:
                next(gout)
                content = next(gout).split()
                while content[0].isnumeric():
                    state_tmp = int(content[0])
                    dip_matrix[0, state_tmp, :] =\
                        dip_matrix[state_tmp, 0, :] =\
                        [float(dip_mom) for dip_mom in content[1:4]]
                    content = next(gout).split()
            # 'Dipole moment (field-independent basis, Debye):' is too long.
            # In order to abide to PEP8 a shorter string is used.
            elif 'Dipole moment (field-independent' in line and get_diag_dip:
                # The values in the Gaussian output are reported in Debye:
                # a conversion from Debye to a.u. is performed.
                dip_matrix[el_state, el_state, :] =\
                    [float(dip_mom)*conv.Debye2au
                     for dip_mom in next(gout).split()[1:6:2]]
            elif 'ESP charges:' in line and get_charges:
                next(gout)
                content = next(gout).split()
                while content[0].isnumeric():
                    charges.append(float(content[2]))
                    content = next(gout).split()
        if get_geom and geom:
            geom_tmp = np.zeros_like(geom)
            geom_tmp[:, :] = geom
            qm_inputs['geometry'] = geom_tmp
        elif get_geom and geom_no_symm:
            geom_tmp = np.zeros_like(geom_no_symm)
            geom_tmp[:, :] = geom_no_symm
            qm_inputs['geometry'] = geom_tmp
        if get_energies:
            energies_tmp = np.zeros_like(energies)
            energies_tmp[:] = energies
            qm_inputs['energies'] = energies_tmp
        if get_diag_dip or get_trans_dip:
            qm_inputs['dip_matrix'] = dip_matrix
        if get_charges:
            charges_tmp = np.zeros_like(charges)
            charges_tmp[:] = charges
            qm_inputs['charges'] = np.array(charges_tmp)
        return qm_inputs


def get_tot_input_gauss(filename_scheme: str, n_el_states: int,
                        get_charges=False) -> dict:
    '''Obtain the electronic properties need for basic MD-PMM calculations
    from calculations performed with Gaussian(16, compatibility with other
    versions still to be verified).

    Parameters:
        filename_scheme (str): A string indicating the scheme used to name the
            n_el_states Gaussian calculations (for each electronic state).
            Must contain "{}" in order to iterate over all the electronic
            states.
        n_el_states (int): A integer referring to the total number of
            electronic states considered (including the ground state).
        get_charges (bool): A boolean that specifies if extracting the RESP
            charges is desired

    Returns:
        qm_inputs (dict): "geometry": geometry in Angstrom (numpy.darray,
                shape=(n_atoms, 4)).
            "energies": electronic states energies in a.u. (numpy.darray,
                shape=n_el_states).
            "dip_matrix": complete electric dipole moment matrix in a.u.
                (numpy.darray, shape=(n_el_states, n_el_states, 3)).
            "charges_i": RESP charges for the ith electronic state
                (numpy.darray, shape=(n_el_states, n_atoms)).
    '''
    if '{}' not in filename_scheme:
        print('filename_scheme must contain "{}"')
        raise IOError
    qm_inputs = {}
    for state in range(n_el_states):
        if state == 0:
            qm_inputs_tmp = get_input_gauss(filename_scheme.format(state),
                                            n_el_states, state,
                                            get_energies=False,
                                            get_trans_dip=False,
                                            get_charges=get_charges)
            qm_inputs['geometry'] = qm_inputs_tmp['geometry']
            qm_inputs['dip_matrix'] = qm_inputs_tmp['dip_matrix']
            if get_charges:
                qm_inputs[f'charges_{state}'] = qm_inputs_tmp['charges']
        elif state == 1:
            qm_inputs_tmp = get_input_gauss(filename_scheme.format(state),
                                            n_el_states, state,
                                            get_geom=False,
                                            get_charges=get_charges)
            qm_inputs['energies'] = qm_inputs_tmp['energies']
            qm_inputs['dip_matrix'] += qm_inputs_tmp['dip_matrix']
            if get_charges:
                qm_inputs[f'charges_{state}'] = qm_inputs_tmp['charges']
        else:
            qm_inputs_tmp = get_input_gauss(filename_scheme.format(state),
                                            n_el_states, state,
                                            get_geom=False,
                                            get_energies=False,
                                            get_charges=get_charges)
            qm_inputs['dip_matrix'][state, state, :] =\
                qm_inputs_tmp['dip_matrix'][state, state, :]
            if get_charges:
                qm_inputs[f'charges_{state}'] = qm_inputs_tmp['charges']
    return qm_inputs


def get_pmm_inputs(cmdline) -> tuple[dict, mda.Universe]:
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
    qm_inputs = get_raw_inputs(cmdline.ref_geom, cmdline.dip_matrix,
                               cmdline.energies, cmdline.charges)
    if '.tpr' in cmdline.topology_path or '.xtc' in cmdline.trajectory_path:
            mm_traj = mda.Universe(cmdline.topology_path,
                                   cmdline.trajectory_path)
    return qm_inputs, mm_traj


def get_pmm_inputs_old(cmdline) -> tuple[dict, mda.Universe]:
    '''NOTE: old version of get_pmm_input. To be reworked in the future.
    Obtain the pmm inputs according to the different sources.

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
    try:
        qm_inputs, mm_traj = legacy_input(cmdline.legacy_input, cmdline)
    except ValueError:
        if cmdline.qm_source.lower() == 'gaussian':
            qm_inputs = get_tot_input_gauss(cmdline.qm_path, cmdline.roots)
        if cmdline.mm_source.lower() == 'gromacs':
            mm_traj = mda.Universe(cmdline.topology_path,
                                   cmdline.trajectory_path)
    return qm_inputs, mm_traj
