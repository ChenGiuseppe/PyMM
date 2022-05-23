Example of MD-PMM simulation using PyMM: Mmlecule of water in water box.

1. Run the MD-PMM simulation:
* in the dipole approximation:
```
pymm run_pmm -g data/ref_geom_ang.dat -gu angstrom -dm data/dipmat.dat -e data/energies.dat -traj data/traj.xtc -top data/topol.dat -nm 1:3 -o TEST_WATER_DIPOLE
```

* in the QC atomic charges:
```
pymm run_pmm -g data/ref_geom_ang.dat -gu angstrom -dm data/dipmat.dat -e data/energies.dat -traj data/traj.xtc -top data/topol.dat -nm 1:3 -ch data/charges.dat -o TEST_WATER_CHARGES
```

2. Plot absorption spectrum:
```
pymm calc_abs -dm data/dipmat.dat -el TEST_WATER_CHARGES_eigvals.dat -ev TEST_WATER_CHARGES_eigvecs.npy
```
3. Perturbed eigenvector projection analysis:
```
pymm eig -first 1 -last 2 -state 1 -i TEST_WATER_CHARGES_eigvecs.npy -oc
```
4. Pertur eigenvectors projection analysis
```
pymm eig -state 1 -i TEST_WATER_CHARGES_eigvecs.npy -ot
```
5. Cumulative histogram considering only the first 3 states.
```
pymm eig -i TEST_WATER_CHARGES_eigvecs.npy -state 3 -oh
```

