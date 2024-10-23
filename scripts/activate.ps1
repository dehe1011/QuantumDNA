Write-Output "Checking if the virtual environment already exists..."
$envExists = Test-Path -Path ".venv_qdna/Scripts/activate"

if ($envExists) {

    Write-Output "Virtual environment found. Activating the environment..."
    & .venv_qdna/Scripts/activate

} else {

    Write-Output "Virtual environment not found. Creating a new virtual environment..."
    python -m venv .venv_qdna
    Write-Output "Activating the newly created environment..."
    & .venv_qdna/Scripts/activate

    Write-Output "Installing package in editable mode..."
    pip install -e .


    $installJupyter = Read-Host "Install Jupyter Notebook ([Y/N])?"

    if ($installJupyter -eq "Y" -or $openGUI -eq "y") {
        Write-Output "Installing Jupyter Notebook..."
        pip install ipykernel
        pip install notebook
    } else {
        Write-Output "Installation skipped."
    }

    Write-Output "Running tests to verify everything is working..."
    pip install pytest
    python -m pytest --disable-pytest-warnings
}


$openGUI = Read-Host "Open Graphical User Interface ([Y/N])?"

if ($openGUI -eq "Y" -or $openGUI -eq "y") {
    Write-Output "Opening the qDNA graphical user interface..."
    python qDNA/gui/qdna_app.py
} else {
    Write-Output "Graphical User Interface launch skipped."
}

$openJupyter = Read-Host "Open Jupyter Notebook ([Y/N])?"

if ($openJupyter -eq "Y" -or $openGUI -eq "y") {
    Write-Output "Opening Jupyter Notebook..."
    jupyter notebook
} else {
    Write-Output "Jupyter Notebook launch skipped."
}
