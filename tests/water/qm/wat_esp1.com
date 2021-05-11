%NProcShared=24
%Chk=wat_esp.chk
#P B3LYP/6-31G(d) Pop=(MK,dipole) Density=Current TD=(Root=1,NStates=3)

 esp test water py

0 1
O          0.00000        0.00000        0.11969
H          0.00000        0.76143       -0.47878
H          0.00000       -0.76143       -0.47878

