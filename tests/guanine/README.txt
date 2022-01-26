#PMM neutral ensemble, neutral guanine
pymm run_pmm -g gua_n/geom_gua_n -gu bohr -dm gua_n/dipmat_gua_n -e gua_n/energies_gua_n -traj ensemble_n/traj_n.xtc -top  ensemble_n/topol_n -ch gua_n/guanine_n_charges.dat -nm 9,11:25 -o TEST_N_GUA_N

#PMM neutral ensemble, radical cation guanine
pymm run_pmm -g gua_rc/geom_gua_rc -gu bohr -dm gua_rc/dipmat_gua_rc -e gua_rc/energies_gua_rc -traj ensemble_n/traj_n.xtc -top  ensemble_n/topol_n -ch gua_rc/guanine_rc_charges.dat -nm 9,11:25 -o TEST_N_GUA_RC

#PMM radical cation ensemble, neutral guanine
pymm run_pmm -g gua_n/geom_gua_n -gu bohr -dm gua_n/dipmat_gua_n -e gua_n/energies_gua_n -traj ensemble_rc/traj_rc.xtc -top  ensemble_rc/topol_rc -ch gua_n/guanine_n_charges.dat -nm 9,11:25 -o TEST_RC_GUA_N

#PMM radical cation ensemble, radical cation guanine
pymm run_pmm -g gua_rc/geom_gua_rc -gu bohr -dm gua_rc/dipmat_gua_rc -e gua_rc/energies_gua_rc -traj ensemble_rc/traj_rc.xtc -top  ensemble_rc/topol_rc -ch gua_rc/guanine_rc_charges.dat -nm 9,11:25 -o TEST_RC_GUA_RC

pymm free_en -T 300 -eii TEST_N_GUA_N_eigvals.dat -efi TEST_N_GUA_RC_eigvals.dat -eif TEST_RC_GUA_N_eigvals.dat -eff TEST_RC_GUA_RC_eigvals.dat

## Risultato dato dal programma in Fortran in results/free_en.log