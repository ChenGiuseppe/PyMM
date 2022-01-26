1. Commands to launch the tests:
* Dipole approximation:
pymm run_pmm -g data/geom -gu bohr -dm data/dipmat -e data/energies -traj data/traj_noPBC_skip1000.xtc -top data/topol -nm 1:16 -o TEST_INDOLE_DIPOLE

* QC atomic charges:
pymm run_pmm -g data/geom -gu bohr -dm data/dipmat -e data/energies -traj data/traj_noPBC_skip1000.xtc -top data/topol -ch indole_charges.dat -nm 1:16 -o TEST_INDOLE_CHARGES

2. Plot absorption spectrum:
pymm calc_abs -dm data/dipmat -el TEST_INDOLE_CHARGES_eigvals.dat -ev TEST_INDOLE_CHARGES_eigvecs.npy

3. Essential analysis
pymm eig -first 1 -last 2 -state 1 -i TEST_INDOLE_CHARGES_eigvecs.npy -oc

4. Total essential analysis
pymm eig -state 1 -i TEST_INDOLE_CHARGES_eigvecs.npy -ot

5. Cumulative histogram
# Consider only the first 4 states.
pymm eig -i TEST_INDOLE_CHARGES_eigvecs.npy -state 4 -oh

--------------------------------------------------------------------------------------
COMPARE TO THE RESULTS OBTAINED IN FORTRAN
* Dipole approximation
Compare the first perturbed excited state energies.

python3 compare_dip.py

* QC atomic charges
Compare the first perturbed excited state energies.

python3 compare_ch.py
