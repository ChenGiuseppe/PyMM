import numpy as np
import pymm.conversions as conv

def calc_dA(T, en_ini, en_fin, state):
    '''
    Calculate the reaction free energy change considering a single
    ensemble.

    Parameters:
        T (float): temperature of the system.
        en_ini (np.ndarray): perturbed energies trajectory of
            the initial state (expressed in atomic units).
        en_fin (np.ndarray): perturbed energies trajectory of
            the final state (expressed in atomic units).
    Returns:
        dA (float): reaction free energy change (expressed in J/mol).
    '''

    # en_ini = np.loadtxt(file_en_ini)[:,]
    # en_fin = np.loadtxt(file_en_fin)[:,]

    delta = en_fin - en_ini

    if state == 'ini':
        sgn = -1
    if state == 'fin':
        sgn = 1

    q = np.exp((sgn*delta*conv.au2eV/(conv.k_B*T))).mean()

    dA = sgn*conv.R*T*np.log(q)

    return dA


def calc_dA_mean(T, en_ini_ens_ini, en_fin_ens_ini, en_ini_ens_fin, en_fin_ens_fin):
    '''
    Calculate the reaction free energy change considered as the mean calculated over
    two ensembles (initial state and final state).
    '''

    dA_ini = calc_dA(T, en_ini_ens_ini, en_fin_ens_ini, state = 'ini')
    dA_fin = calc_dA(T, en_ini_ens_fin, en_fin_ens_fin, state = 'fin')

    dA_mean = 0.5*(dA_ini + dA_fin)

    return dA_mean

if __name__ == '__main__':
   pass