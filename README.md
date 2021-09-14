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
conda install numpy matplotlib mdtraj MDAnalysis scipy
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
PyMM -h
```

It will list the commands currently supported. To launch these programs type:

```
PyMM <command>
```

or:

```
PyMM <program> -h
```

to obtain info on how to use them.

### run_pmm

The main feature of PyMM is to perform MD-PMM calculations (options in "[]" are optional):

```
PyMM run_pmm -g geometry -gu [units] -dm dipole_matrix -e energies -ch [QC_atomic_charges] -traj traj.xtc -top traj.tpr -q [QC_charge] -nm 1:3 -o [eigenval.txt] -oc [eigvecs]
```

* -g : QC geometry input file. Each line should be formatted as:
 
``` 
[atom symbol] [x] [y] [z] 
```

> The program also accept other types of formatting and will check and correct some errors.