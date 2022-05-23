# Test: Guanosine
Input files to calculate __guanosine__ (radical cation/neutral) redox potential in aqueous solution using MD-PMM with PyMM. The QC used is __guanine__ in the radical cationic and neutral form. The atom-based expansion approximation has been used.

| WARNING: The length of the MD simulation has been reduced for practical reasons. The electric dipole matrix lacks the elements associated with the i --> j transitions (where i, j are the electronic states and i, j =/= 0) due to the limitations of the software used to calculate the unperturbed properties. Since for the modelling of the redox process only the ground state is needed and since no relevant mixing is observed between the unperturbed ground state and the excited states, such approximation doesn't significantly influence the results of the calculation. |
| --- |

1. Run the MD-PMM simulations for guanosine in both oxidized and reduced form in both ensembles.

* Neutral ensemble, neutral guanine QC:

> pymm run_pmm -g gua_n/geom_gua_n.dat -gu bohr -dm gua_n/dipmat_gua_n.dat -e gua_n/energies_gua_n.dat -traj ensemble_n/traj_n.xtc -top ensemble_n/topol_n.dat -ch gua_n/guanine_n_charges.dat -nm 9,11:25 -o TEST_N_GUA_N

* Neutral ensemble, radical cation guanine QC:

> pymm run_pmm -g gua_rc/geom_gua_rc.dat -gu bohr -dm gua_rc/dipmat_gua_rc.dat -e gua_rc/energies_gua_rc.dat -traj ensemble_n/traj_n.xtc -top  ensemble_n/topol_n.dat -ch gua_rc/guanine_rc_charges.dat -nm 9,11:25 -o TEST_N_GUA_RC

* Radical cation ensemble, neutral guanine QC:

> pymm run_pmm -g gua_n/geom_gua_n.dat -gu bohr -dm gua_n/dipmat_gua_n.dat -e gua_n/energies_gua_n.dat -traj ensemble_rc/traj_rc.xtc -top  ensemble_rc/topol_rc.dat -ch gua_n/guanine_n_charges.dat -nm 9,11:25 -o TEST_RC_GUA_N

* Radical cation ensemble, radical cation guanine QC:

> pymm run_pmm -g gua_rc/geom_gua_rc.dat -gu bohr -dm gua_rc/dipmat_gua_rc.dat -e gua_rc/energies_gua_rc.dat -traj ensemble_rc/traj_rc.xtc -top  ensemble_rc/topol_rc.dat -ch gua_rc/guanine_rc_charges.dat -nm 9,11:25 -o TEST_RC_GUA_RC

2. Calculate the reaction free energy change:

> pymm free_en -T 300 -eii TEST_N_GUA_N_eigvals.dat -efi TEST_N_GUA_RC_eigvals.dat -eif TEST_RC_GUA_N_eigvals.dat -eff TEST_RC_GUA_RC_eigvals.dat