# Use the following two commands on the Anaconda Powershell Prompt to run this file:

# 1. Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA"
# For example: Set-Location -Path "C:\Users\Dennis Herb\OneDrive\2_Uni\Doktor\python_projects\QuantumDNA"
# 2. powershell -ExecutionPolicy Bypass -File tools\scripts\run_tests.ps1

# -----------------------------------------------------------

Write-Output "Running tests"
# -m: module, -s: start of search, -v: verbose 
# python -m unittest discover -s tests -v
python -m pytest -vv tests/ --disable-pytest-warnings

# -------------------------------------------------------

Write-Output "Removing unnecessary files and directories created in setup"

# Specify the directories that should be deleted
$paths = @(
    "tests/__pycache__",
    "tests/.ipynb_checkpoints"
)

foreach ($path in $paths) {

    if (Test-Path -Path $path) {
    
        Remove-Item -Path $path -Recurse -Force
        Write-Output "Removed $path"
        
    } else {
    
        Write-Output "$path does not exist"
        
    }
}

# ------------------------------------------------------