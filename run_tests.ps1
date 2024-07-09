# Use the following two commands on the Anaconda Power Shell
# Set-Location -Path "C:\Users\Dennis Herb\OneDrive\Dokumente\2. Uni\Doktor\Python Scripts\Quantum_DNA_1.0"
# .\run_tests.ps1

# -------------------------------------------------------

echo "Running tests"
# -m: module, -s: start of search, -v: verbose 
# python -m unittest discover -s tests -v
python -m pytest -vv tests/ --disable-pytest-warnings

# -------------------------------------------------------

echo "Removing unnecessary files and directories created in setup"

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