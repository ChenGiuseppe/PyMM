import numpy as np
import matplotlib.pyplot as plt

py = np.loadtxt('TEST_N_GUA_RC_eigvals.dat')

#fort =  np.loadtxt('/mnt/d/Dottorato/Programs/test_guanosine/ensemble_n/gua_n/en_pert_gua_n')
fort =  np.loadtxt('results/output_n_gua_rc')

x_py = np.arange(py.shape[0])
x_fort = np.arange(fort.shape[0])

plt.plot(x_py, py[:,0], label='PyMM')

plt.plot(x_fort, fort[:,1], label='fortran')

print("Average energy of the perturbed ground state:\n",
      "PyMM:", py[:,0].mean(), "\nFortran:",  fort[:,1].mean())


plt.legend()
plt.show()
