# Use the following two commands on the Anaconda Powershell Prompt to run this file:

# 1. Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA"
# 2. powershell -ExecutionPolicy Bypass -File scripts\activate.ps1

# -----------------------------------------------------------

# Check if the conda environment 'qDNA' exists
$envExists = conda info --envs | Select-String -Pattern "qDNA"

if ($envExists) {

    # Activate the existing conda environment
    conda activate qDNA

} else {

    conda create --name qDNA --file your_file.txt

    conda activate qDNA
    conda install pip
    pip install jupyter
    pip install e .

    # Create a new kernel that can be selected inside Jupyter notebooks:
    python -m ipykernel install --name qDNA --display-name "Python (qDNA)"

    # Run all the tests to make sure that everything works:
    Write-Output "Running tests"
    pip install pytest
    python -m pytest -vv tests/ --disable-pytest-warnings
}

# Ask the user if they want to open the graphical user interface
$openGUI = Read-Host "Open Graphical User Interface ([Y/N])?"

if ($openGUI -eq "Y" -or $openGUI -eq "y") {
    # Open the graphical user interface
    python qDNA\gui\qdna_app.py
}
