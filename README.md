# CMEE MRes Project

This directory contains the codes and data for producing the results in the author's MRes research project thesis. All scripts were written and tested using **R 4.4.1**, **Python 3.12.3**, **Visual Studio Code 1.94.2**, and the **bash terminal** on **Ubuntu 24.04.1 LTS**.

## Project structure

- `code/`: Contains codes for running on HPC and Jupyter Notebook for analysis. 
- `data/`: Statistics outputs transfered from HPC that will be use for analysis in the Notebooks.
- `outputs/`: An empty directory for storage of images produced from the codes.

---

## Usage

### Running the scripts
- **aestivation3.c**

This is not standalone C but rather a "plugin" to be mounted and recognised by R via an API. 

In terminal, run `R CMD SHLIB aestivation3.c -fopenmp`. This will create a shared object with the same name but in .so extension (.dll in Windows). The shared object can then be mounted to R. The code to do it is already included in `LD.R`.

- **hpc_*.sh**

All the shell scripts with the hpc_ prefix are for batch jobs on PBS HPC system. Transfer them along with all the Python and R scripts to the `$HOME` directory on HPC, then submit the jobs by `qsub` command. E.g. `qsub hpc_LD.sh`.

### Dependencies
To ensure Python and R are set up on your HPC system, you need to first enter the command lines on the remote system:

```
module load anaconda3/personal
anaconda-setup
```

After anaconda is installed, install R by

```
conda install R
```

Then set up a virtual environment by
```
conda create -n msprime_env msprime tskit numpy pandas r-tidyverse r-reticulate 
```
This will create an environment `msprime_env` with all necessary modules and packages installed.

---

## Author
Laiyin Zhou
l.zhou24@imperial.ac.uk


## Acknowledgment
Tin-Yu Hui provided codes for unlinked loci forward simulator (`aestivation3.c` and the `sim_pop` function in `LD.R`), as well as analysis codes adapted by the author into `Season_SampleSize.ipynb`.
