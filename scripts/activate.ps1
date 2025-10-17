# Function to prompt with a default value
function Prompt-Default {
    param(
        [string]$Message,
        [string]$Default
    )
    $response = Read-Host "$Message [$Default]"
    if ([string]::IsNullOrWhiteSpace($response)) {
        return $Default
    }
    return $response
}

# Check if Python 3 is installed
if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    Write-Error "Python is not installed. Please install Python and try again."
    exit 1
}

# Check if pip is installed
if (-not (Get-Command pip -ErrorAction SilentlyContinue)) {
    Write-Error "pip is not installed. Please install pip and try again."
    exit 1
}

# Check if the virtual environment exists
if (Test-Path ".venv/Scripts/Activate.ps1") {
    Write-Output "Virtual environment found. Activating the environment..."
    . .\.venv\Scripts\Activate.ps1
} else {
    Write-Output "Virtual environment not found."

    # Ask the user if they want to create a new virtual environment (default: Yes)
    $createVenv = Prompt-Default "Do you want to create a new virtual environment (Y/N)" "Y"
    if ($createVenv -match "^(Y|y)$") {
        Write-Output "Creating a new virtual environment..."
        python -m venv .venv
        Write-Output "Activating the newly created environment..."
        . .\.venv\Scripts\Activate.ps1

        Write-Output "Installing package in editable mode..."
        pip install --upgrade pip
        pip install -e .

        # Ask the user if they want to install Jupyter Notebook
        $installJupyter = Prompt-Default "Install Jupyter Notebook (Y/N)" "N"
        if ($installJupyter -match "^(Y|y)$") {
            Write-Output "Installing Jupyter Notebook..."
            pip install ipykernel
            pip install notebook
        } else {
            Write-Output "Jupyter Notebook installation skipped."
        }

        Write-Output "Running tests to verify everything is working..."
        pip install pytest
        python -m pytest --disable-pytest-warnings
    } else {
        Write-Output "Virtual environment creation skipped. Exiting..."
        exit 1
    }
}

# Option to open the GUI
$openGUI = Prompt-Default "Open Graphical User Interface (Y/N)" "N"
if ($openGUI -match "^(Y|y)$") {
    Write-Output "Opening the qDNA graphical user interface..."
    python Scripts/open_gui.py
} else {
    Write-Output "Graphical User Interface launch skipped."
}

# Option to open Jupyter Notebook
$openJupyter = Prompt-Default "Open Jupyter Notebook (Y/N)" "N"
if ($openJupyter -match "^(Y|y)$") {
    Write-Output "Opening Jupyter Notebook..."
    jupyter notebook
} else {
    Write-Output "Jupyter Notebook launch skipped."
}
