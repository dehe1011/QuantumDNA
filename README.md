# Quantum_DNA

**Author: Dennis Herb**

Release date: hopefully in 2024 

## To Do

* Add more tests
* Improve tutorials and docstrings
* Improve README file
* Publish package on PyPI?
* Implement disorder and correlations (Mirko)
* Include conductance
* Reproduce plots from other papers

## What's new

* 09.07.: Added a user interface to the package such that it is easier accessible for users that are less familiar with computer science.

## Getting started 

**NOTE**: These set-up instructions have only been tested on Windows and may not work on macOS. 

### Installation 

Copy the following ``commands`` into Anaconda Power Shell 

1. Clone the Github repository (downloads all files and folders from the Github project):\
`` git clone https://github.com/dehe1011/quantum_DNA.git ``\
`` cd quantum_DNA ``

2. Create and active a virtual environment and create a new kernel that can be selected in Jupyter Notebooks:\
`` powershell -ExecutionPolicy Bypass -File ./activate.ps1 ``

If all tests worked you have successfully installed the package and the user interface opens automatically. You can access all the implemented functionalities. Enjoy :)

After you have already installed the package, there are two possibilities to access the code:

1. Reopen the user interface: navigate to the directory of the package and use the following command:
`` powershell -ExecutionPolicy Bypass -File ./activate.ps1 ``
2. Open a new Jupyter notebook and select the kernel to "Python (qDNA)":
`` jupyter notebook ``

To remove the virtual environment use `` conda remove --name qDNA --all ``

### Shortcuts

The following shortcuts are frequently used in the code: 

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
* seq: sequence
* calc: calculate

### Tutorials

In the ``quantum_DNA_1.0/doc/tutorials/`` directory, there exist the following tutorials:

* **0_TB_Model.ipynb**: tutorials on predefined and custom tight-binding models.
* **1_Plotting_Routines**: tutorial on the predefined plotting routines.
* **2_Open_System**: tutorial on different ways to treat DNA relaxation the DNA environment.
* **Exploration_Notebook**: if you are interested how single functions behave, this is the place where you can test the functionality isolated from the rest of the code. 

In the ``quantum_DNA_1.0/doc/`` directory you can find:
* **Paper**: Notebook to reproduce all the plots contained in [D. Herb, M. Rossini and J. Ankerhold, Ultrafast excitonic dynamics in DNA: Bridging correlated quantum dynamics and sequence dependence.](https://arxiv.org/abs/2402.16892)
* **Produce_Data**: In this notebook longer simulations are performed to obtain data (e.g., the lifetimes and average charge separation of up to seven base pairs)

## References

* [R. Siebert, O. Ammerpohl, M. Rossini et al. A quantum physics layer of epigenetics: a hypothesis deduced from charge transfer and chirality-induced spin selectivity of DNA. *Clin Epigenet 15*, 145 (2023).](https://doi.org/10.1186/s13148-023-01560-3)
* [D. Herb, M. Rossini and J. Ankerhold, Ultrafast excitonic dynamics in DNA: Bridging correlated quantum dynamics and sequence dependence.](https://arxiv.org/abs/2402.16892)

## Notes

* If there occur any unexpected errors or problems please contact the author via dennis.herb@uni-ulm.de



