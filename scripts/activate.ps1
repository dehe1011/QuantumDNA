# Use the following two commands on the Anaconda Powershell Prompt to run this file:

# 1. Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA"
# 2. powershell -ExecutionPolicy Bypass -File scripts\activate.ps1

# -----------------------------------------------------------

Write-Output "Checking if the conda environment 'qDNA' exists..."
$envExists = conda info --envs | Select-String -Pattern "qDNA"

if ($envExists) {

    Write-Output "Conda environment 'qDNA' found. Activating the environment..."
    conda activate qDNA

} else {

    Write-Output "Conda environment 'qDNA' not found. Creating a new environment from requirements.txt..."
    conda create --name qDNA --file requirements/requirements.txt
    Write-Output "Activating the newly created environment..."
    conda activate qDNA

    Write-Output "Installing package in editable mode..."
    pip install e .

    Write-Output "Checking if Jupyter Notebook is installed globally..."
    $jupyterExists = (jupyter --version) -ne $null

    if (-not $jupyterExists) {
        Write-Output "Jupyter Notebook not found globally. Installing Jupyter Notebook inside the virtual environment..."
        pip install notebook
    } else {
        Write-Output "Jupyter Notebook is available globally."
    }

    Write-Output "Installing Jupyter kernel for the 'qDNA' environment..."
    pip install ipykernel
    python -m ipykernel install --name qDNA --display-name "Python (qDNA)"

    Write-Output "Running tests to verify everything is working..."
    pip install pytest
    python -m pytest -vv tests/ --disable-pytest-warnings
}

$openGUI = Read-Host "Open Graphical User Interface ([Y/N])?"

if ($openGUI -eq "Y" -or $openGUI -eq "y") {
    Write-Output "Opening the qDNA graphical user interface..."
    python qDNA\gui\qdna_app.py
} else {
    Write-Output "Graphical User Interface launch skipped."
}
