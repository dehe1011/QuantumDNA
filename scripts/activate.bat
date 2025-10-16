@echo off

:: Check if Python 3 is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo Python 3 is not installed. Please install Python and try again.
    exit /b
)

:: Check if pip is installed
pip --version >nul 2>&1
if errorlevel 1 (
    echo pip is not installed. Please install pip and try again.
    exit /b
)

:: Check if the virtual environment exists
if exist ".venv\Scripts\activate.bat" (
    echo Virtual environment found. Activating the environment...
    call .venv\Scripts\activate.bat
) else (
    echo Virtual environment not found.

    :: Ask the user if they want to create a new virtual environment (default: Yes)
    set /p createVenv="Do you want to create a new virtual environment (Y/N)? [Y]: "
    if "%createVenv%"=="" set createVenv=Y
    if /i "%createVenv%"=="Y" (
        echo Creating a new virtual environment...
        python -m venv .venv
        echo Activating the newly created environment...
        call .venv\Scripts\activate.bat

        echo Installing packages...
        pip install --upgrade pip
        pip install -e .

        set /p installJupyter="Install Jupyter Notebook (Y/N)? [N]: "
        if "%installJupyter%"=="" set installJupyter=N
        if /i "%installJupyter%"=="Y" (
            pip install ipykernel notebook
        ) else (
            echo Skipping Jupyter Notebook installation.
        )

        echo Running tests...
        pip install pytest
        python -m pytest --disable-warnings
    ) else (
        echo Skipping virtual environment setup. Exiting.
        exit /b
    )
)

set /p openGUI="Open Graphical User Interface (Y/N)? [N]: "
if "%openGUI%"=="" set openGUI=N
if /i "%openGUI%"=="Y" (
    python Scripts/open_gui.py
) else (
    echo Skipping GUI launch.
)

set /p openJupyter="Open Jupyter Notebook (Y/N)? [N]: "
if "%openJupyter%"=="" set openJupyter=N
if /i "%openJupyter%"=="Y" (
    jupyter notebook
) else (
    echo Skipping Jupyter Notebook launch.
)
