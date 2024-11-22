#!/bin/bash

# Function to prompt with a default value
prompt_default() {
    local message=$1
    local default=$2
    read -p "$message [$default]: " response
    echo "${response:-$default}"
}

# Check if Python 3 is installed
if ! command -v python3 &>/dev/null; then
    echo "Python 3 is not installed. Please install Python and try again."
    exit 1
fi

# Check if pip is installed
if ! command -v pip3 &>/dev/null; then
    echo "pip is not installed. Please install pip and try again."
    exit 1
fi

# Check if the virtual environment exists
if [ -f ".venv/bin/activate" ]; then
    echo "Virtual environment found. Activating the environment..."
    source .venv/bin/activate
else
    echo "Virtual environment not found."
    create_venv=$(prompt_default "Do you want to create a new virtual environment (Y/N)" "Y")
    if [[ "$create_venv" =~ ^[Yy]$ ]]; then
        echo "Creating a new virtual environment..."
        python3 -m venv .venv
        source .venv/bin/activate

        echo "Installing packages..."
        pip install --upgrade pip
        pip install -e .

        install_jupyter=$(prompt_default "Install Jupyter Notebook (Y/N)" "N")
        if [[ "$install_jupyter" =~ ^[Yy]$ ]]; then
            pip install ipykernel notebook
        else
            echo "Skipping Jupyter Notebook installation."
        fi

        echo "Running tests..."
        pip install pytest
        python3 -m pytest --disable-warnings
    else
        echo "Skipping virtual environment setup. Exiting."
        exit 1
    fi
fi

open_gui=$(prompt_default "Open Graphical User Interface (Y/N)" "N")
if [[ "$open_gui" =~ ^[Yy]$ ]]; then
    if [ -f "Scripts/open_gui.py" ]; then
        python3 Scripts/open_gui.py
    else
        echo "GUI script not found."
    fi
else
    echo "Skipping GUI launch."
fi

open_jupyter=$(prompt_default "Open Jupyter Notebook (Y/N)" "N")
if [[ "$open_jupyter" =~ ^[Yy]$ ]]; then
    jupyter notebook
else
    echo "Skipping Jupyter Notebook launch."
fi
