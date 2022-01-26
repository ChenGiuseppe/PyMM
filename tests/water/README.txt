fort.251 frame, autovettori

fort.67  frame, angoli
	 matrice di rotazione

fort.600 struttura rfin ottenuta dal fit di r0_2 su r0_1

fort.60 struttura iniziale r0_1

fort.33 Total CPU time

output_ch energia ground perturbata   

-------------------------------------------------

1. Commands to launch the tests:
* Dipole approximation:
pymm run_pmm -g data/geom_ang -gu angstrom -dm data/dipmat -e data/energies -traj data/trajout_l.xtc -top data/topol -nm 1:3 -o TEST_WATER_DIPOLE

* QC atomic charges:
pymm run_pmm -g data/geom_ang -gu angstrom -dm data/dipmat -e data/energies -traj data/trajout_l.xtc -top data/topol -nm 1:3 -ch data/water_charges.dat -o TEST_WATER_CHARGES

2. Plot absorption spectrum:
pymm calc_abs -dm data/dipmat -el TEST_WATER_CHARGES_eigvals.dat -ev TEST_WATER_CHARGES_eigvecs.npy

3. Essential analysis
pymm eig -first 1 -last 2 -state 1 -i TEST_WATER_CHARGES_eigvecs.npy -oc

4. Total essential analysis
pymm eig -state 1 -i TEST_WATER_CHARGES_eigvecs.npy -ot

5. Cumulative histogram
# Consider only the first 3 states.
pymm eig -i TEST_WATER_CHARGES_eigvecs.npy -state 3 -oh

--------------------------------------------------------------------------------------
COMPARE TO THE RESULTS OBTAINED IN FORTRAN

* Dipole approximation
Compare the first perturbed excited state energies.

python3 compare_dip.py

* QC atomic charges
Compare the first perturbed excited state energies.

python3 compare_ch.py
