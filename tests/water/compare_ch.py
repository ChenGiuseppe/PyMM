import numpy as np
import matplotlib.pyplot as plt

py = np.loadtxt('TEST_WATER_CHARGES_eigvals.dat')

fort =  np.loadtxt('results/output_ch.dat')

x_py = np.arange(py.shape[0])
x_fort = np.arange(fort.shape[0])

plt.plot(x_py, py[:,1], label='PyMM')

plt.plot(x_fort, fort[:,2], label='fortran')

print("Average energy of the 1st perturbed excited state:\n",
      "PyMM:", py[:,1].mean(), "\nFortran:",  fort[:,2].mean())


plt.legend()
plt.show()
