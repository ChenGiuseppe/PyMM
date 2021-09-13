# pmm
Computational chemistry package to perform MD-PMM calculations

## Installation of required libraries

We recommend the installation of the required libraries via [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
Under the section 'Latest Miniconda Installer Links' chose the suitable installer file for your os.
Once the installer is downloaded you can proceed with the [installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
### For Linux:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

For creating a new environment called pmm:

```
conda create --name pmm
```

To use, or "activate" the new environment, type the following:

```
conda activate pmm
```

Now we are ready to install all the packages required for the pmm program

```
conda install numpy mdtraj MDAnalysis ...
```

