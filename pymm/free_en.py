import numpy as np

def calc_dA(T, en_ini, en_fin, state):
    '''
    The variable state indicates the state of the
    ensemble over the average is performed.
    '''

    # en_ini = np.loadtxt(file_en_ini)[:,]
    # en_fin = np.loadtxt(file_en_fin)[:,]

    delta = en_fin - en_ini

    if state == 'ini':
        sgn = -1
    if state == 'fin':
        sgn = +1

    q = np.exp((sgn*delta*27.211386*(10**5)/(8.617333*T))).mean()


    # Boltzmann's constant in eV/K
    # k_B = 8.617333*10^(-5) eV/K

    R = 8.314 # J/(K*mol)
    dA = sgn*R*T*np.log(q)

    return dA

def calc_dA_mean(T, en_ini_ens_ini, en_fin_ens_ini, en_ini_ens_fin, en_fin_ens_fin):

    dA_ini = calc_dA(T, en_ini_ens_ini, en_fin_ens_ini, state = 'ini')
    dA_fin = calc_dA(T, en_ini_ens_fin, en_fin_ens_fin, state = 'fin')

    dA_mean = 0.5*(dA_ini + dA_fin)

    return dA_mean

if __name__ == '__main__':
    calc_dA_mean(278, 'en_n_ens_n', 'en_rc_ens_n', 'en_n_ens_rc', 'en_rc_ens_rc')