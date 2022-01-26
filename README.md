# PyMM

PyMM is a program that allows you to easily run MD-PMM [[1]](#1) (an _a posteriori_ QM/MM method) simulations. Furthermore, it includes tools for the analysis of the dynamical behaviour of the electronic states and properties during the simulation, and for the prediction of experimental properties such as the UV-Vis spectrum and the redox potential.

## Installation

We recommend the installation of the required libraries via [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
Choose the suitable installer file for your system under the section _Latest Miniconda Installer Links_.
Once the installer is downloaded you can proceed with the installation by following [these instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

If you are using Linux just run the following command:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Now create a new environment (e.g. "pmm"):
```
conda create --name pmm
```

Activate the new environment:
```
conda activate pmm
```

Now we are ready to install all the packages required to run PyMM:
```
conda install numpy matplotlib scipy numba
```
```
conda config --add channels conda-forge
```
```
conda install mdanalysis
```

Move to the folder that contains setup.py and use the command:
```
cd pymm-path/pymm/
```

```
pip install .
```
<br></br>
## Usage

Before trying to launch the program, be sure that you have activated the correct environment (the one where the required libraries are installed). You can then proceed using the PyMM program.

If you need help use:

```
pymm -h
```

It will list the modules currently implemented in PyMM. To launch these programs type:

```
pymm <command>
```

or:

```
pymm <program> -h
```

to obtain info on how to use them.

### run_pmm

The main feature of PyMM is to perform MD-PMM calculations. This can be done using the **run_pmm** program (arguments in "[ ]" are optional):

```
PyMM run_pmm -g geometry -gu [units] -dm dipole_matrix -e energies -ch [QC_atomic_charges] -traj traj.xtc -top traj.tpr -q [QC_charge] -nm 1:3 -o [eigvals.txt] -oc [eigvecs]
```

#### INPUTS

* **-g** : QC reference geometry input file. Each line should be formatted as:
 
```
<atom symbol> <x> <y> <z> 
```

> **NOTE**: The program also accept other types of formatting and will check and correct some common errors. This is done mainly for compatibility with legacy software.

* -gu : units used in the QC geometry input file (default: Angstrom).

* **-dm** : electric dipole matrix file (in a.u.). Each line should be formatted as:

```
<state_n> <state_m> <x> <y> <z>
```

* **-e** : electronic energies file (in a.u.). Each line should be formatted as:

```
<state_n energy>
```

* -ch : file containing QC atomic charges for each electronic state. It should be formatted as follows:

```
<atom1 state1 charge> <atom2 state1 charge> ...
<atom1 state2 charge> <atom2 state2 charge> ...
```
>**NOTE**: -ch is optional. When it is provided, the MD-PMM calculation will be performed by expanding the perturbation operator on each of atom of the QC [[2]](#2). If the QC atomic charges are not provided, **run_pmm** will run the calculation using the dipole approximation [[1]](#1).

* **-traj** : XTC file containing the MD simulation trajectory.

* **-top** : TPR file corresponding to the XTC file. It's necessary to obtain the topology and the MM charges.
> **NOTE**: As an alternative, a text file listing the indexes and the charge of the atoms in the MD simulation can be used as the topology file.

* -q : QC total charge. By default is 0.

* **-nm** : select the atoms of the QC in the MD simulation according to the atom indexes (indexes taken from the MD simulation). Indexes start from 1. 

>**Selection Algebra**:
>1. Use "**:**" to indicate a range.
\
Example 1: select atoms from 1 to 5: 
\
-nm 1:5
>2. Use "**,**" as the logical *or* operator.
\
> Example 2: select atoms 1 or 5:
\
-nm 1,5
>3. Example 3: select atoms from 1 to 5 and 7:
\
-nm 1:5,7

> **NOTE**: Don't leave any whitespaces between the indexes. 

> **NOTE**: the selection system follows the <em>bysum index-range</em> selection mode of [MDAnalysis](https://docs.mdanalysis.org/stable/documentation_pages/selections.html), a library that PyMM relies on.

* --match: reorder the QC reference geometry to match the atoms order in the MD simulation. It is useful when the geometry used for the quantum mechanical calculation doesn't match the order of the MD simulation.
> **NOTE**: not matching the QC atoms is a common source of error.

#### OUTPUT
* -o: job name that will be used as a prefix to all the output files. (Default: pymm). The programs will output: QC geometry as an xyz file, eigenvalues and eigenvectors. 

<br></br>

### eig
**eig** is a tool to analyse the perturbed electronic states corresponding eigenvectors after the MD-PMM calculation. Given the dynamic nature of these states during the simulation, a straightforward description can be difficult. Therefore we provide three graphical representations which combined can present a clearer picture to the user.

#### **ESSENTIAL ANALYSIS**
The contributions arising from two unperturbed states, to a selected unperturbed state, are plotted against each other.

Inputs

* -i: eigenvectors trajectory (eigenvecs.npy).
* -first: select the index of the an unperturbed state.
* -last: select the index of another unperturbed state.
* -state: select the perturbed state to be considered.
* -oc: select to plot the essential analysis.


#### **COMPLETE ESSENTIAL ANALYSIS**
The contributions to a selected perturbed state, arising from all the unperturbed states, are plotted against each other.

Inputs

* -i: eigenvectors trajectory (eigenvecs.npy).
* -state: perturbed state to be considered. 
* -ot: select to plot the complete essential analysis.

#### **CUMULATIVE HISTOGRAMS**
The cumulative histograms of the average of contribution of each unperturbed state to each perturbed state is plotted (i.e. calculate the squared mean coefficients with each unperturbed state contributes to each perturbed state).

Inputs

* -i: eigenvectors trajectory (eigenvecs.npy).
* -state: select number of states to consider
* -oh: select to calculate the cumulative histograms

<br></br>
> **NOTE**: Despite its semplicity in depicting the overall trends during the simulation, considering the average contribution can be quite deceiving. This representation is preferably to be paired with an essential analysis.

### calc_abs

**calc_abs** is a program used to calculate the absorption spectrum of the perturbed QC. After running the MD-PMM calculation you can obtain the absorption spectrum from the output files using:

```
pymm calc_abs -dm dipmat -el eigvals.dat -ev eigvecs.npy -sigma [0.0003] -ot [abs_spectrum]
```

#### INPUTS
* **-dm**: electric dipole matrix file (in a.u.) used in the MD-PMM calculation.
* -el: trajectory of the perturbed eigenvalues.
* -ev: trajectory of the perturbed eigenvectors.
* -sigma: value of the gaussian broadening applied to the signal (expressed in eV). Default: 0.05 eV.

#### OUTPUTS
* -ot: calculated absorption spectra names prefix (default: abs_spectrum). Both the total spectrum and each individual transitions will be printed.
<br></br>
### free_en

Calculate the free energy difference between two states (initial and final) each considered in the two ensembles (of the initial and final state) using the Zwanzig formula.

```
pymm free_en -T 298 -eii file1 -efi file2 -eif file3 -eff file4
```

#### INPUTS
* -T: temperature (in Kelvin) at which the simulation is carried. Default: 298 K.
* **-eii**: eigenvalues obtained from the MD-PMM calculation performed in the initial state ensemble and considering the QC in its initial state.
* **-efi**: eigenvalues obtained from the MD-PMM calculation performed in the initial state ensemble and considering the QC in its final state.
* **-eif**: eigenvalues obtained from the MD-PMM calculation performed in the final state ensemble and considering the QC in its initial state.
* **-eff**: eigenvalues obtained from the MD-PMM calculation performed in the final state ensemble and considering the QC in its final state.

| ensemble\state | initial |  final  | 
|----------------|---------|---------|
|    initial     |   eii   |   efi   |
|    final       |   eif   |   eff   |


Three different approaches can be adopted to the calculation of free energy in the MD-PMM framework.

1. Considering only the initial ensemble: -eii and -efi need to be provided.
2. Considering only the final ensemble: -eif and -eff need to be provided.
3. Considering the average between the two ensembles: -eii, -efi, -eif and -eff are all necessary for the calculation.

According to the inputs provided, the corresponding model will be selected. For a better comprehension of the implemented approches, please consult previously published literature [[3]] [[4]].


<br></br>

## Glossary

* **Quantum center** (QC): portion of the system to be described on a quantum-mechanical level.
<br></br>

## References
<a id="1">[1]</a>
Aschi, M., Spezia, R., Di Nola, A., & Amadei, A. (2001). A first-principles method to model perturbed electronic wavefunctions: the effect of an external homogeneous electric field. _Chemical physics letters_, 344(3-4), 374-380, https://doi.org/10.1016/S0009-2614(01)00638-8.

<a id="2">[2]</a>
Zanetti-Polzi, L., Del Galdo, S., Daidone, I., D'Abramo, M., Barone, V., Aschi, M., & Amadei, A. (2018). Extending the perturbed matrix method beyond the dipolar approximation: comparison of different levels of theory. _Physical Chemistry Chemical Physics_, 20(37), 24369-24378, https://doi.org/10.1039/C8CP04190C.

