Installation
============

**NOTE**: These set-up instructions have only been tested on Windows.

Installation via PyPI
---------------------

We recommend creating a new virtual environment and Jupyter notebook kernel to avoid conflicts with existing packages.

1. Open the Windows Powershell and navigate to your project folder. Create a new virtual environment:

   .. code-block:: bash

       python -m venv .venv

2. Activate the virtual environment:

   .. code-block:: bash

       .venv/Scripts/activate.ps1

3. Install the qDNA package:

   .. code-block:: bash

       pip install qDNA

*Optional*: Use `qDNA` inside a Jupyter notebook:

   .. code-block:: bash

       pip install ipykernel
       jupyter notebook

Installation via Cloning the GitHub Repository
----------------------------------------------

If you want to make changes to the source code, you can clone the project's GitHub repository and install it in editable mode. On Windows, you first need to install Git (`Download <https://gitforwindows.org/>`_).

1. Clone the GitHub repository:

   .. code-block:: bash

       git clone https://github.com/dehe1011/QuantumDNA.git

2. Navigate to the directory of the cloned repository:

   .. code-block:: bash

       cd QuantumDNA

3. Use the provided activation script to finish the installation. Optionally, you can open the Graphical User Interface or a Jupyter Notebook:

   .. code-block:: bash

       powershell -ExecutionPolicy Bypass -File scripts/activate.ps1

If all tests passed, the package has been successfully installed, and you can access all the implemented functionalities. Enjoy!

Usage
-----

After installing the package, you can access the code via the Graphical User Interface or in a Jupyter Notebook simply by running the activation script again:

.. code-block:: bash

    powershell -ExecutionPolicy Bypass -File scripts/activate.ps1

Uninstallation
--------------

Uninstall the package by running:

.. code-block:: bash

    pip uninstall qDNA

If you cloned the GitHub repository, remember to manually delete the `QuantumDNA` folder from your computer.
