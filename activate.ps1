# Use the following two commands on the Anaconda Power Shell
# Set-Location -Path "C:\Users\Dennis Herb\OneDrive\Dokumente\2. Uni\Doktor\Python Scripts\Quantum_DNA_1.0"
# .\activate.ps1

# git clone https://github.com/dehe1011/quantum_DNA.git
# cd quantum_DNA

# Check if the conda environment 'qDNA' exists
$envExists = conda info --envs | Select-String -Pattern "qDNA"

if ($envExists) {
    # Activate the existing conda environment
    conda activate qDNA
} else {
    # Clone the repository and create the conda environment
    conda env create -f environment.yml
    conda activate qDNA
    python -m ipykernel install --name qDNA --display-name "Python (qDNA)"
    .\run_tests.ps1
}

python qdna_app.py
