<p align="center">
  <picture>
    <img src="stored_data/doc_images/qDNA_logo.png">
  </picture>
</p>

<div align="center">
  <span>
    <img src="https://img.shields.io/badge/license-MIT-blue" alt="PyPI - License">
  </span>
  <span>
    <a href="https://doi.org/10.5281/zenodo.12734027">
      <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.12734027.svg" alt="DOI">
    </a>
  </span>
</div>

---

# QuantumDNA

**Author: Dennis Herb**

This Python package can be cited as:

> *QuantumDNA (github.com/dehe1011/QuantumDNA)*, D. Herb, 2024, DOI: [10.5281/zenodo.12734027](https://doi.org/10.5281/zenodo.12734027)


## To Do

* Add more tests
* Improve tutorials and docstrings
* Improve README file
* Publish package on PyPI?
* Implement disorder and correlations (Mirko)
* Include conductance
* Reproduce plots from other papers


## What's new

**July 2024**

* Added a graphical user interface (GUI) to the package such that it is easier accessible for users that are less familiar with computer science. The user interface is based on the 
* Added a Jupyter Notebook ```quantum_DNA_1.0/doc/Paper.ipynb``` that reproduces all figures contained in our paper (and the supplementary) [D. Herb, M. Rossini and J. Ankerhold, *Physical Review E 109*, 064413 (2024).](https://doi.org/10.1103/PhysRevE.109.064413)


## Introduction

### Abstract

After photoexcitation of DNA, the excited electron (in the LUMO) and the remaining hole (in the HOMO) localized on the same DNA base form a bound pair, called the Frenkel exciton, due to their mutual Coulomb interaction. In this study, we demonstrate that a tight-binding (TB) approach, using TB parameters for electrons and holes available in the literature, allows us to correlate relaxation properties, average charge separation, and dipole moments to a large ensemble of double-stranded DNA sequences (all 16384 possible sequences with 14 nucleobases). This way, we are able to identify a relatively small subensemble of sequences responsible for long-lived excited states, high average charge separation, and high dipole moment. Further analysis shows that these sequences are particularly T rich. By systematically screening the impact of electron-hole interaction (Coulomb forces), we verify that these correlations are relatively robust against finite-size variations of the interaction parameter, not directly accessible experimentally. This methodology combines simulation methods from quantum physics and physical chemistry with statistical analysis known from genetics and epigenetics, thus representing a powerful bridge to combine information from both fields.


## Getting started 

**NOTE**: These set-up instructions have only been tested on Windows and may not work on macOS. 

**Pre-requisites**

* Conda ([Download Anaconda](https://www.anaconda.com/download) or [Download Miniconda](https://docs.anaconda.com/miniconda/))
* Git ([Download](https://gitforwindows.org/) )

### Installation 

Open the Anconda Powershell Prompt. Copy and execute the following ```commands```.

1. Clone the Github repository (downloads all files and folders from the Github project):\
``` git clone https://github.com/dehe1011/quantum_DNA.git```\
```cd quantum_DNA ```

2. Create and active a virtual environment and create a new kernel that can be selected in Jupyter Notebooks:\
``` powershell -ExecutionPolicy Bypass -File ./activate.ps1 ```

If all tests worked you have successfully installed the package and the user interface opens automatically. You can access all the implemented functionalities. Enjoy :)

After you have already installed the package, you can either access the code via the user interface or in a Jupyter Notebook. Open the Anconda Powershell Prompt again and use the following ``` commands ```:

1. Open the user interface: 

    (i) Navigate to the directory of the package (insert the location of the package): 
    ``` Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA ``` \
    (ii) ``` powershell -ExecutionPolicy Bypass -File ./activate.ps1 ```

2. Open a new Jupyter notebook: 

    (i) ``` jupyter notebook ``` \
    (ii) Select the kernel to "Python (qDNA)"

To remove the virtual environment use ``` conda remove --name qDNA --all ```

### Shortcuts

To increase the readability of the code I collected some of the shortcuts that are used requently: 

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



## Documentation

### Example Program

To test QuantumDNA you can run the following simple example where the exciton lifetime and the average charge separation of a double-stranded GCG DNA sequence are calculated. You can try different sequences, tight-binding models and keyword arguments to investigate if and how strong the exciton lifetime and average charge separation are affected. For example you might find that in general more uniform sequences show higher values. This is because their energy landscape is less distorted leading to less localization.  

```python

from DNA import calc_lifetime, calc_dipole

# input
upper_strand = 'GCG'
tb_model_name = 'ELM'
kwargs = dict(unit='rad/ps', relax_rate=3, source='Hawke2010')

# calculation 
lifetime = calc_lifetime(upper_strand, tb_model_name, **kwargs)
dipole = calc_dipole(upper_strand, tb_model_name, **kwargs)

# output 
print(f"Exciton lifetime {lifetime} fs")
print(f"Average charge separation {dipole} A")

```

### Tutorials

The code contains some tutorials and demostrations to better understand and explore the functionalities. 

In the ```quantum_DNA_1.0/doc/``` folder you can find the notebook **Paper.ipynb** that contains and reproduces all the figures contained in [D. Herb, M. Rossini and J. Ankerhold, Ultrafast excitonic dynamics in DNA: Bridging correlated quantum dynamics and sequence dependence.](https://arxiv.org/abs/2402.16892)

In the ```quantum_DNA_1.0/doc/tutorials/ ``` folder, there exist the following tutorials:

* **0_TB_Model.ipynb**: tutorials on predefined and custom tight-binding models.
* **1_Plotting_Routines.ipynb**: tutorial on the predefined plotting routines.
* **2_Open_System.ipynb**: tutorial on different ways to treat DNA relaxation the DNA environment.
* **Exploration_Notebook.ipynb**: if you are interested how single functions behave, this is the place where you can test the functionality isolated from the rest of the code. 


### Graphical user interface

The usage of the graphical user interface is demonstrated in the following images. Many functionalities of the code can be accessed in a very user friendly manner from the menu window:

![](stored_data/doc_images/user_interface_doc/menu_1.png)
> _Screenshot of the menu of the user interface._

![](stored_data/doc_images/user_interface_doc/plot_1.png)
> _Plot obtained after pressing the submit button on the menu (see image above)._

![](stored_data/doc_images/user_interface_doc/menu3.png)
> _Screenshot of the menu of the user interface with calculations of the exciton lifetime and average charge separation displayed in the frame on the bottom right._


## References

Papers from our group:

* [R. Siebert, O. Ammerpohl, M. Rossini et al. A quantum physics layer of epigenetics: a hypothesis deduced from charge transfer and chirality-induced spin selectivity of DNA. *Clin Epigenet 15*, 145 (2023).](https://doi.org/10.1186/s13148-023-01560-3)
* [D. Herb, M. Rossini and J. Ankerhold, Ultrafast excitonic dynamics in DNA: Bridging correlated quantum dynamics and sequence dependence. *Physical Review E 109*, 064413 (2024).](https://doi.org/10.1103/PhysRevE.109.064413)

Tight-binding parameters:

* [L.G.D. Hawke, G. Kalosakas and C. Simserides, Electronic parameters for charge transfer along DNA. *The European Physical Journal E 32*, 291 (2010)](https://doi.org/10.1140/EPJE/I2010-10650-Y)
* [C. Simserides, A systematic study of electron or hole transfer along DNA dimers, trimers and polymers. *Chemical Physics 440*, 31 (2014)](https://doi.org/10.1016/j.chemphys.2014.05.024)
* [M. Mantela, C. Simserides and R. Di Felice, LCAO electronic structure of nucleic acid bases and other heterocycles and transfer integrals in B-DNA, including structural variability. *Materials 14*, 4930 (2021)]

Tight-binding models:

* [K. Lambropoulos and C. Simserides, Tight-binding modeling of nucleic acid sequences: Interplay between various types of order or disorder and charge transport. *Symmetry 11*, 968 (2019)](https://doi.org/10.3390/sym11080968)

Excitons and electron-hole Coulomb interaction:

* [C. Crespo-Hernandez, B. Cohen and B. Kohler, Base stacking controls excited-state dynamics in A·T DNA. *Nature 436*, 1141 (2005)]()
* [E.R. Bittner, Lattice theory of ultrafast excitonic and charge-transfer dynamics in DNA. *Journal of Chemical Physics 125*, 094909 (2006)](https://doi.org/10.1063/1.2335452)
* [E.M. Conwell, P.M. McLaughlin and S.M. Bloch, Charge-Transfer Excitons in DNA. *The Journal of Physical Chemistry B 112*, 2268 (2008)](https://doi.org/10.1021/jp077344x)
* [S. Tornow, R. Bulla, F.B. Anders and G. Zwicknagl, Multiple-charge transfer and trapping in DNA dimers. *Physical Review B 82*, 195106 (2010)](https://doi.org/10.1103/PhysRevB.82.195106)

Biological relevance of DNA charge transfer:

* [J.C. Genereux and J.K. Barton, Mechanisms for  DNA charge transport. *Chemical Reviews 110*, 1642 (2010)](https://doi.org/10.1021/cr900228f)
* [A.R. Arnold, M.A. Grodick and J.K. Barton, DNA Charge Transport: from Chemical Principles to the Cell. *Cell Chemical Biology 23*, 183 (2016)](https://doi.org/10.1016/j.chembiol.2015.11.010)

Simulation of open quantum systems:

* [quantum_HEOM (github.com/jwa7/quantum_HEOM), J.W. Abbott, 2022](https://doi.org/10.5281/zenodo.7230160)



## Notes

* If there occur any unexpected errors or problems please contact the author via dennis.herb@uni-ulm.de
