# Installation Guide

## Installation via PyPI:

We recommend creating a new virtual environment and Jupyter notebook kernel to avoid conflicts with existing packages like `qutip`. 
1. Create a new virtual environment (using conda):\
```conda create -n qDNA```

2. Create a new Jupyter kernel (if you are using Jupyter notebook):\
 ```ipykernel install --name qDNA --display-name "Python (qDNA)"```

3. Install the qDNA package:\
```pip install qDNA```


## Installation via Cloning the Github Repository

If you want to make changes to the source code or try the example notebooks, you can clone the project's GitHub repository.

**NOTE**: These set-up instructions have only been tested on Windows and may not work on macOS. 

### Pre-requisites

* Conda ([Download Anaconda](https://www.anaconda.com/download) or [Download Miniconda](https://docs.anaconda.com/miniconda/))
* Git ([Download](https://gitforwindows.org/))
* Python ([Download](https://www.python.org/downloads/))

### Installation procedure

Open the Anconda Powershell Prompt. Copy and execute the following ```commands```.

1. Clone the Github repository: \
```git clone https://github.com/dehe1011/QuantumDNA.git```

2. Navigate to the project directory: \
```cd QuantumDNA ```

2. Create and active a virtual environment (using a provided script): \
```powershell -ExecutionPolicy Bypass -File tools/scripts/activate.ps1```

If all tests passed, the package has been successfully installed, and the user interface opens automatically. You can access all the implemented functionalities. Enjoy!

### Usage

After installing the package, you can access the code via the user interface or in a Jupyter Notebook.

Using the User Interface: 

1. Open the Anaconda PowerShell Prompt and navigate to the package directory:\
 ```Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA``` 

2. Run the activation script: \
```powershell -ExecutionPolicy Bypass -File tools/scripts/activate.ps1```

Using Jupyter Notebook: 

1. Open a new Jupyter Notebook: \
```jupyter notebook```

2. Select the kernel "Python (qDNA)":
    * In the Jupyter Notebook interface, go to Kernel > Change kernel > Python (qDNA)


## Uninstallation

To remove the virtual environment and Jupyter kernel:

1. Remove the virtual environment:\
 ```conda remove --name qDNA --all``` 

2. Remove the Jupyter kernel:\
 ```jupyter kernelspec remove qDNA``` 
 
3. Delete the project folder:
    * Manually delete the `QuantumDNA` folder that contains the cloned GitHub repository.