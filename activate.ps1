# Use the following two commands on the Anaconda Power Shell
# Set-Location -Path "C:\Users\Dennis Herb\OneDrive\Dokumente\2. Uni\Doktor\Python Scripts\Quantum_DNA_1.0"
# .\activate.ps1

# Check if the conda environment 'qDNA' exists
$envExists = conda info --envs | Select-String -Pattern "qDNA"

if ($envExists) {
    # Activate the existing conda environment
    conda activate qDNA
} else {
    # Create a virtual environment from a .yml file that contains name, channels and dependencies:
    conda env create -f environment.yml
    # Activate the virtual environment that you just created:
    conda activate qDNA
    # Create a new kernel that can be selected inside Jupyter notebooks:
    python -m ipykernel install --name qDNA --display-name "Python (qDNA)"
    # Run all the tests to make sure that everything works:
    powershell -ExecutionPolicy Bypass -File ./run_tests.ps1
}

# Open the user interface
python user_interface\qdna_app.py

# Open a new Jupyter notebook and select the kernel to "Python (qDNA)"
# jupyter notebook

# To remove the virtual environment use conda remove --name qDNA --all