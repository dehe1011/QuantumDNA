# Quantum_DNA

**Author: Dennis Herb**

Release date: hopefully in 2024 

## To Do

* Upload on Github and add a .gitignore file
* Write tests, add data types, add docstrings with examples
* create tutorials 

## What's new

* This is where I want to keep track of major changes to the project. 

## Getting started 

**NOTE**: These set-up instructions have only been tested on Windows and may not work on macOS. 

### Installation 

Copy the following ``commands`` into Anaconda Power Shell 

1. Clone the Github repository (downloads all files and folders from the Github project):\
`` git clone https://github.com/dehe1011/quantum_DNA.git ``\
`` cd quantum_DNA ``
3. Create a virtual environment from a ``.yml`` file that contains name, channels and dependencies:\
`` conda env create -f environment.yml ``
4. Activate the virtual environment that you just created:\
`` conda activate qDNA ``\
Optional: check if the virtual environment was successfully installed:\
`` conda info --envs ``
5. Create a new kernel that can be selected inside Jupyter notebooks:\
`` python -m ipykernel install --name qDNA --display-name "Python (qDNA)" ``
6. Run all the tests to make sure that everything works:\
`` .\run_tests.ps1 ``
7. Open a new Jupyter notebook and select the kernel to "Python (qDNA)":
`` jupyter notebook ``

To remove the virtual environment use `` conda remove --name qDNA --all ``

### Shortcuts
 
* ham: hamiltonian
* dm: density matrix
* tb: tight-binding
* eigv: eigenvalue/ eigenenergy
* eigs: eigenstates/ eigenvectors
* dim: dimension
* fig: figure
* op: operator
* loc: local
* glob: global
* deph: dephasing
* therm: thermalizing

### Tutorials

In the ``quantum_DNA_1.0/doc/tutorials/`` directory, there exists the following tutorials:

* **0_TB_Model.ipynb**: tutorials on predefined and custom tiught-binding models


## References

* [R. Siebert, O. Ammerpohl, M. Rossini et al. A quantum physics layer of epigenetics: a hypothesis deduced from charge transfer and chirality-induced spin selectivity of DNA. *Clin Epigenet 15*, 145 (2023).](https://doi.org/10.1186/s13148-023-01560-3)
* [D. Herb, M. Rossini and J. Ankerhold, Ultrafast excitonic dynamics in DNA: Bridging correlated quantum dynamics and sequence dependence.](https://arxiv.org/abs/2402.16892)

## Notes

Create a link from the Python environment's site-packages directory (where packages like numpy or matplotlib are saved) to the /path/to/my_package. You can edit the package's source code (indicated by -e):
`` pip install -e path/to/package ``



