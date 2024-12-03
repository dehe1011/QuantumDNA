.. _guide:

**************************
Jupyter Notebook Tutorials
**************************

Welcome to the **QuantumDNA Jupyter Notebook Tutorials**! These tutorials and demonstrations are designed to help users explore and
understand the functionalities of the `qDNA` package through practical examples.

For a detailed description of the classes and functions used in these tutorials, visit the :ref:`API documentation <apidoc>`.

Overview
========

The Jupyter notebooks in the `tutorials` folder provide step-by-step demonstrations of various features and functionalities of the `qDNA` package.
Each tutorial focuses on a specific aspect, helping you get started with quantum physical description of DNA and applications of quantum biology.

Tutorials
=========

Below is a list of available tutorials. Click on each to access the detailed notebook:

.. toctree::
   :maxdepth: 1

   PRE_2024 <tutorials/PRE_2024>
   Tight Binding Parameters <tutorials/1_Tight_Binding_Parameters>
   Tight Binding Method <tutorials/2_Tight_Binding_Method>
   Environment Simulation <tutorials/3_Environment_Simulation>
   Visualization <tutorials/4_Visualization>
   Evaluation <tutorials/5_Evaluation>
   Reproduce Papers <tutorials/6_Reproduce_Papers>

Descriptions
------------

**PRE2024**
  Reproduces all the figures presented in the reference paper :cite:`Herb2024`.
  This serves as a comprehensive example of `qDNA`'s visualization and analysis features.

**Tight Binding Parameters**
  Learn the Linear Combination of Atomic Orbitals (LCAO) approach using Slater–Koster two-center
  integrals and Harrison-type expressions. Ideal for tight-binding model parameterization.

**Tight_Binding_Method**
  Explore predefined and custom tight-binding models. Includes calculating time-averaged exciton
  populations in the Fishbone Ladder Model (FLM) and simulating charge transfer in the
  Fenna-Matthews-Olson (FMO) complex.

**Environment_Simulation**
  Model DNA excited-state relaxation and environmental interactions. This tutorial covers dephasing
  and thermalization models inspired by Quantum Biology.

**Visualization**
  Use `qDNA`'s built-in plotting routines for effective result visualization. Learn to create custom
  visualizations tailored to your data.

**Evaluation**
  Perform calculations for observables like exciton lifetimes, average charge separation, and dipole
  moments. Includes parallelization features for efficient computation.

**Reproduce Papers**
  Reproduce the plots from the papers :cite:`Giese1999`, :cite:`Giese2001`, :cite:`Bittner2006`,
  :cite:`Bittner2007`, :cite:`Simserides2014` and :cite:`Mantela2023` using the `qDNA` package. This tutorial
  demonstrates the package's capabilities in generating results that have already been published elsewhere.


Getting Started
===============

These tutorials provide hands-on examples designed to guide you through using the `qDNA` package effectively. To get started:

1. Navigate to the `tutorials` folder and open the desired `.ipynb` file in Jupyter Notebook or JupyterLab.
2. Follow the instructions provided in the notebook to run the cells and explore the package's features interactively.
3. Refer to the :ref:`API documentation <apidoc>` for deeper insights into the functions and classes used in the tutorials.


Tips for Using the Tutorials
----------------------------

- **Run in a Jupyter Environment**: Ensure you have Jupyter Notebook or JupyterLab installed to execute the tutorials interactively.
- **Dependencies**: Before starting, confirm that all dependencies for `qDNA` are installed. Check the `requirements.txt` file in the repository for details.
- **Explore Further**: Modify and experiment with the code to deepen your understanding of the concepts.


We hope these tutorials help you leverage the full potential of the `qDNA` package for your quantum biology research!
