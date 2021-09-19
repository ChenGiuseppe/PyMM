# PyMM

Computational chemistry package to perform MD-PMM calculations

## Installation of required libraries

We recommend the installation of the required libraries via [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
Under the section 'Latest Miniconda Installer Links' choose the suitable installer file for your OS.
Once the installer is downloaded you can proceed with the [installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### For Linux:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Create a new environment (e.g. "pmm"):

```
conda create --name pmm
```

Activate the new environment:

```
conda activate pmm
```

Now we are ready to install all the packages required for the PyMM program:

```
conda install numpy matplotlib scipy numba MDAnalysis
```

Move to the folder that contains setup.py and use the command:

```
pip install .
```

Since the program is in its early stages and changes are frequent, if you are a Git user we reccomend installing using this command instead:

```
pip install -e .
```

## Usage

Before trying to launch the program be sure that you have activated the correct environment (the one where the required libraries are installed).
You can then proceed usign the PyMM program.
If you need help use:

```
pymm -h
```

It will list the commands currently supported. To launch these programs type:

```
pymm <command>
```

or:

```
pymm <program> -h
```

to obtain info on how to use them.

### run_pmm

The main feature of PyMM is to perform MD-PMM calculations (arguments in "[]" are optional):

```
PyMM run_pmm -g geometry -gu [units] -dm dipole_matrix -e energies -ch [QC_atomic_charges] -traj traj.xtc -top traj.tpr -q [QC_charge] -nm 1:3 -o [eigvals.txt] -oc [eigvecs]
```

* **-g** : QC geometry input file. Each line should be formatted as:
 
```
<atom symbol> <x> <y> <z> 
```

> **NOTE**: The program also accept other types of formatting and will check and correct some common errors.

* -gu : units used in the QC geometry input file.

* **-dm** : electric dipole matrix file. Each line should be formatted as:

```
<state_n> <state_m> <x> <y> <z>
```

* **-e** : electronic energies file. Each line should be formatted as:

```
<state_n energy>
```

* -ch : file containing QC atomic charges for each electronic state. It should be formatted as follows:

```
<atom1 state1 charge> <atom2 state1 charge> ...
<atom1 state2 charge> <atom2 state2 charge> ...
```

* **-traj** : XTC file containing the MD simulation trajectory.

* **-top** : TPR file corresponding to the XTC file. It's necessary to obtain the topology and the MM charges.

* -q : QC total charge. By default is 0.

* **-nm** : select the atoms of the QC in the MD simulation according to the atom indexes. Indexes start from 1. 

> Example 1: select atoms from 1 to 5: 

```
-nm 1:5
```

> Example 2: select atoms 1 and 5:

```
-nm "1 5"
```

> Example 3: select atoms from 1 to 5 and 7:

```
-nm "1:5 7"
```

> **NOTE**: after -nm, when more than one atom range selection is specified, it's necessary to write the selection between the quotation marks "" (see Example 2 and 3). 

> **NOTE**: the selection system follows the <em>bysum index-range</em> selection mode of [MDAnalysis](https://docs.mdanalysis.org/stable/documentation_pages/selections.html), a library that PyMM relies on.

* -o : filename of the output containing the eigenvalues. By default it is "eigvals.txt".

* -oc : filename of the output containing the eigenvectors (saved in the .npy format). By default it is "eigvecs" (i.e. <em>eigvecs.npy</em>).


### calc_abs

Program used to calculate the absorption spectrum of the QC. After running the MD-PMM calculation you can obtain the absorption spectrum using:

```
pymm calc_abs -dm dipmat -el eigvals.txt -ev eigvecs.npy -sigma [0.0003] -ot [abs_spectrum]
```

### free_en

Calculate the free energy difference between two states (initial and final) each considered in the two ensembles (of the initial and final state). The Zwanzig formula was used.

```
pymm free_en -T 298 -eii file1 -efi file2 -eif file3 -eff file4
```

> **NOTE**: The energy files are to be provided according to these scheme:
>
> | ensemble\state | initial |  final  | 
> |----------------|---------|---------|
> |    initial     |   eii   |   efi   |
> |    final       |   eif   |   eff   |
