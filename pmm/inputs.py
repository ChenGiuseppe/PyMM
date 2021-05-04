import numpy as np

# Dictionary to help convert atomic numbers to atomic masses. Declared here in case it is needed for different functions; position
# in the code to be still evaluated. Only the most common ones are listed, to be expanded in the future.
Z_to_mass = {'1': 1.0, '6': 12.0, '7': 14.0, '8': 16.0, '16': 32.0}

def get_inputs_gaussian(filename: str, n_el_states: int, el_state: int, get_geom=True, get_energies=True, get_dip_matrix=True, get_charges=False) -> dict:
    '''Extract the geometry, the electronic states energies, the elements of the electric dipole moment matrix and/or
    the RESP charges from the Gaussian (specifically 16, support for other versions is to be verified) software 
    calculation outputs.
    The permanent dipole moment and the RESP charges extracted refer to the only electronic state considered in the Gaussian calculaton
    (ground state or excited state).
    All the returned properties are expressed in atomic units.

    Parameters:
        filename (str): A string indicating the name of the Gaussian calculation output file.
        n_el_states (int): A integer referring to the total number of electronic states considered (including the ground state).
        el_state (int): A integer indicating the electronic state considered in the calculation (starting from 0 for the ground state).
        get_geom (bool): A boolean that specifies if extracting the geometry is wanted.
        get_energies (bool): A boolean that specifies if extracting the electronic states energies is wanted.
        get_dip_matrix (bool): A boolean that specifies if extracting the elements of the electric dipole moment matrix is wanted.
        get_charges (bool): A boolean that specifies if extracting the RESP charges is desired.

    Returns:
        pmm_inputs (dict): "geometry": geometry (numpy.darray, shape=(n_atoms, 4)).
            "energies": electronic states energies (numpy.darray, shape=n_el_states).
            "dip_matrix": electric dipole moment matrix (numpy.darray, shape=(n_el_states, n_el_states, 3)). INCOMPLETE: values only for 
            the (0,i), (i,0) and the (el_state, el_state) elements.
            "charges": RESP charges for the el_state electronic state (numpy.darray, shape=(n_el_states, n_atoms)).
    ''' 
    pmm_inputs = {}
    with open(filename, 'r') as gout:
        geom = []
        energies = np.zeros(n_el_states)
        dip_matrix = np.zeros((n_el_states, n_el_states,3))
        charges = []
        for line in gout:
            if 'Standard orientation:' in line and get_geom:
                for i in range(4):
                    next(gout)
                content = next(gout).split()
                while '---' not in content[0]:
                    geom.append([Z_to_mass[content[1]]]+[float(coor) for coor in content[3:6]])
                    content = next(gout).split()
            elif 'SCF Done:  E(' in line and get_energies:
                energies[0] = float(line.split()[4])
            elif 'Excited State' in line and get_energies:
                # check if this condition is sufficient to cover all the possible cases
                # we take the excitation energy expressed in eV, which is then converted to a.u.
                energies[int(line.split()[2][:-1])] = energies[0] + float(line.split()[4])/27.211324570273
            elif 'Ground to excited state transition electric dipole moments (Au):' in line and get_dip_matrix:
                next(gout)
                content = next(gout).split()
                while content[0].isnumeric():
                    dip_matrix[0,int(content[0]),:] = dip_matrix[int(content[0]),0,:] = [float(dip_mom) for dip_mom in content[1:4]]
                    content = next(gout).split()
            elif 'Dipole moment (field-independent basis, Debye):' in line and get_dip_matrix:
                dip_matrix[el_state,el_state,:] = [float(dip_mom) for dip_mom in next(gout).split()[1:6:2]]
            elif 'ESP charges:' in line and get_charges:
                next(gout)
                content = next(gout).split()
                while content[0].isnumeric():
                    charges.append(float(content[2]))
                    content = next(gout).split()
        if get_geom:
            pmm_inputs['geometry'] = np.array(geom)
            print(pmm_inputs['geometry'].shape)
        if get_energies:
            pmm_inputs['energies'] = energies
            print(pmm_inputs['energies'].shape)
        if get_dip_matrix:
            pmm_inputs['dip_matrix'] = dip_matrix
            print(pmm_inputs['dip_matrix'].shape)
        if get_charges:
            pmm_inputs['charges'] = np.array(charges)
            print(pmm_inputs['charges'].shape)
        return pmm_inputs

def get_tot_inputs_gaussian(filename_scheme: str, n_el_states: int, get_charges=False) -> dict:
    '''Obtain the electronic properties need for basic MD-PMM calculations from calculations performed 
    with Gaussian(16, compatibility with other versions still to be verified).

    Parameters:
        filename_scheme (str): A string indicating the scheme used to name the n_el_states Gaussian calculations (for each
            electronic state). Must contain "{}" in order to iterate over all the electronic states.
        n_el_states (int): A integer referring to the total number of electronic states considered (including the ground state).
        get_charges (bool): A boolean that specifies if extracting the RESP charges is desired.
    
    Returns:
        pmm_inputs (dict): "geometry": geometry (numpy.darray, shape=(n_atoms, 4)).
            "energies": electronic states energies (numpy.darray, shape=n_el_states).
            "dip_matrix": complete electric dipole moment matrix (numpy.darray, shape=(n_el_states, n_el_states, 3)).
            "charges_i": RESP charges for the ith electronic state (numpy.darray, shape=(n_el_states, n_atoms)).
    '''
    if '{}' not in filename_scheme:
        print('filename_scheme must contain "{}"')
        raise IOError
    pmm_inputs = {}
    for state in range(n_el_states):
        if state == 0:
            pmm_inputs_tmp = get_inputs_gaussian(filename_scheme.format(state), n_el_states, state, get_energies=False, get_charges=get_charges)
            pmm_inputs['geometry'] = pmm_inputs_tmp['geometry']
            pmm_inputs['dip_matrix'] = pmm_inputs_tmp['dip_matrix']
            if get_charges:
                pmm_inputs[f'charges_{state}'] = pmm_inputs_tmp['charges']
        elif state == 1:
            pmm_inputs_tmp = get_inputs_gaussian(filename_scheme.format(state), n_el_states, state, get_geom=False, get_charges=get_charges)
            pmm_inputs['energies'] = pmm_inputs_tmp['energies']
            pmm_inputs['dip_matrix'][state,state,:] = pmm_inputs_tmp['dip_matrix'][state,state,:]
            if get_charges:
                pmm_inputs[f'charges_{state}'] = pmm_inputs_tmp['charges']
        else:
            pmm_inputs_tmp = get_inputs_gaussian(filename_scheme.format(state), n_el_states, state, get_geom=False, get_energies=False, get_charges=get_charges)
            pmm_inputs['energies'] = pmm_inputs_tmp['energies']
            pmm_inputs['dip_matrix'][state,state,:] = pmm_inputs_tmp['dip_matrix'][state,state,:]
            if get_charges:
                pmm_inputs[f'charges_{state}'] = pmm_inputs_tmp['charges']
    return pmm_inputs


if __name__ == '__main__':
    get_inputs_gaussian('../../../Photoswitch/Norbornadiene/QC_td5.log',11, 5, get_charges=True)