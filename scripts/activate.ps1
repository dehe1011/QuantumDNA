# Use the following two commands on the Anaconda Powershell Prompt to run this file:

# 1. Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA"
# For example: Set-Location -Path "C:\Users\Dennis Herb\OneDrive\2_Uni\Doktor\python_projects\QuantumDNA"
# 2. powershell -ExecutionPolicy Bypass -File scripts\activate.ps1

# -----------------------------------------------------------

# Check if the conda environment 'qDNA' exists
$envExists = conda info --envs | Select-String -Pattern "qDNA"

if ($envExists) {

    # Activate the existing conda environment
    conda activate qDNA

} else {

    # Create a virtual environment from a .yml file that contains name, channels and dependencies:
    conda env create -f configs/environment.yml

    # Activate the virtual environment that you just created:
    conda activate qDNA

    # Create a new kernel that can be selected inside Jupyter notebooks:
    python -m ipykernel install --name qDNA --display-name "Python (qDNA)"

    # Run all the tests to make sure that everything works:
    powershell -ExecutionPolicy Bypass -File scripts\run_tests.ps1
}

# Ask the user if they want to open the graphical user interface
$openGUI = Read-Host "Open Graphical User Interface ([Y/N])?"

if ($openGUI -eq "Y" -or $openGUI -eq "y") {
    # Open the graphical user interface
    python qDNA\gui\qdna_app.py
}

# remove virtual env with
# conda info --envs
# conda remove --name qDNA --all

# remove jupyter kernely
# jupyter kernelspec list
# jupyter kernelspec remove qdna
