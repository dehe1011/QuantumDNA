
# QuantumDNA Installation Guide

Welcome to the installation guide for `QuantumDNA`. Follow the steps below to install the package, set up a virtual environment, and start using `qDNA` either through a Graphical User Interface or a Jupyter Notebook.

---

## Prerequisites

Before proceeding, make sure the following prerequisites are met:

1. **Python Installation**:
   - Verify that Python is installed by running:

     ```bash
     python --version
     ```

   - If Python is not installed, download and install it from [python.org](https://www.python.org/downloads/).

2. **Linux/macOS Users**:
   - Ensure `pip` and `tkinter` are installed:

     ```bash
     python3 -m pip --version
     python3 -m tkinter
     ```

   - If they are not installed, install them using your package manager:
     - **Debian/Ubuntu**:

       ```bash
       sudo apt install python3-pip python3-tk
       ```

     - **Fedora**:

       ```bash
       sudo dnf install python3-pip python3-tkinter
       ```

3. **Windows Users**:
   - Ensure Git is installed by running:

     ```bash
     git --version
     ```

   - If Git is not installed, download and install it from [Git for Windows](https://gitforwindows.org/).

## Installation via PyPI

The easiest way to install `qDNA` is through PyPI. For best results, we recommend creating a new virtual environment to avoid package conflicts.

### Steps

1. **Create a New Virtual Environment**:

    Open your terminal and navigate to your project folder. Run:

    ```bash
    python -m venv .venv
    ```

2. **Activate the Virtual Environment**:
    - **Windows**:

      ```bash
      .venv\Scripts\activate
      ```

    - **macOS/Linux**:

      ```bash
      source .venv/bin/activate
      ```

3. **Install the `qDNA` Package**:

    ```bash
    pip install qDNA
    ```

4. **Optional: Use `qDNA` Inside a Jupyter Notebook**:
    Install Jupyter and launch the notebook:

    ```bash
    pip install ipykernel notebook
    jupyter notebook
    ```

---

## Installation via Cloning the GitHub Repository

If you plan to contribute to the development or make changes to the source code, install `qDNA` in editable mode by cloning its GitHub repository.

### Steps

1. **Clone the GitHub Repository**:

    ```bash
    git clone https://github.com/dehe1011/QuantumDNA.git
    ```

2. **Navigate to the Cloned Repository**:

    ```bash
    cd QuantumDNA
    ```

3. **Run the Activation Script**:
    Use the provided activation script to complete the installation. Instructions vary by platform (see below).

---

## Platform-Specific Instructions for Activation

### **Windows**

1. Navigate to the project directory:

    ```powershell
    Set-Location -Path "C:\Users\<YourUsername>\QuantumDNA"
    ```

2. Run the activation script:

    ```powershell
    powershell -ExecutionPolicy Bypass -File scripts\Activate.ps1
    ```

---

### **macOS**

1. Navigate to the project directory:

    ```bash
    cd /Users/<YourUsername>/QuantumDNA
    ```

2. Run the activation script:

    ```bash
    source scripts/activate
    ```

---

### **Linux**

1. Navigate to the project directory:

    ```bash
    cd /home/<YourUsername>/QuantumDNA
    ```

2. Run the activation script:

    ```bash
    source scripts/activate
    ```

---

## Post-Installation and Usage

If all tests pass, the package has been successfully installed! You can now:

- Launch the **Graphical User Interface** or
- Start using `qDNA` inside a **Jupyter Notebook**.

Run the activation script as mentioned in the platform-specific instructions to start the **Graphical User Interface** or a **Jupyter Notebook**. It is recommended to always run the activation script.

---

## Uninstallation

To remove the package:

```bash
pip uninstall qDNA
```

If you cloned the GitHub repository, manually delete the `QuantumDNA` folder from your computer.

---

🎉 **Congratulations!** You’ve successfully installed and set up `QuantumDNA`. Enjoy exploring the physics of DNA with this powerful tool.
