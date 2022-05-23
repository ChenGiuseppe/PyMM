# Test: Doxorubicin
Input files to calculate the absorption spectrum of __Doxorubicin__ in water. __1,4-dihydroxy-5-methoxy-9,10-anthraquinone__ has been used as QC.

| WARNING: The length of the MD simulation has been reduced for practical reasons. |
| --- |

1. Run the MD-PMM simulation for the Doxorubicin in water:

> pymm run_pmm -g data/ref_geom.dat -gu bohr -dm data/dipmat.dat -e data/energies.dat -traj data/traj.xtc -top data/topol.tpr -nm 17:44 -o TEST_DX

2. Calculate the absorption spectrum of Doxorubicin:

> pymm calc_abs -dm data/dipmat.dat -el TEST_DX_eigvals.dat -ev TEST_DX_eigvecs.npy -sigma 0.14 -ot TEST_DX

3. Analyse the perturbed wavefunction:

* Total projection analysis of the first perturbed electronic excited state:

> pymm eig -i TEST_DX_eigvecs.npy -state 1 -ot

* Projection analysis of the first perturbed electronic excited state on the first and second unperturbed states:

> pymm eig -i TEST_DX_eigvecs.npy -state 1 -first 1 -last 2 -oc

* Cumulative histogram (considering only the first 5 states):

> pymm eig -i TEST_DX_eigvecs.npy -state 5 -oh
