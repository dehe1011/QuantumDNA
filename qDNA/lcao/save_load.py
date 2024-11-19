import json
import os

import pandas as pd

from .. import DATA_DIR, DEFAULTS

__all__ = ["load_xyz", "convert_pdb_to_xyz", "convert_json_to_xyz"]

# ----------------------------- JSON -----------------------------


def load_json(filepath):
    """Loads a JSON file."""
    try:
        with open(filepath, "r") as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        message = f"File {filepath} does not exist."
        if DEFAULTS["verbose"]:
            print(message)


PARAMETRIZATION = load_json(
    os.path.join(
        DATA_DIR,
        "raw",
        "lcao_params",
        DEFAULTS["lcao_default"]["parametrization"] + ".json",
    )
)

# ----------------------------- XYZ -----------------------------


def find_xyz_files(directory):
    """List all .xyz files in the given directory."""
    try:
        files = os.listdir(directory)
        xyz_files = [
            os.path.splitext(file)[0] for file in files if file.endswith(".xyz")
        ]
        return xyz_files
    except Exception as e:
        print(f"An error occurred: {e}")


def load_xyz(filename, directory=os.path.join(DATA_DIR, "geometries")):
    """Loads a XYZ file."""

    filepath = os.path.join(directory, filename + ".xyz")

    with open(filepath, "r") as file:
        # Skip the first line (number of atoms) and second line (comment)
        lines = file.readlines()

    xyz_identifier = lines[1].split()[0]
    xyz_data = pd.DataFrame(
        [line.split() for line in lines[2:]], columns=["Atom", "X", "Y", "Z"]
    )
    xyz_data[["X", "Y", "Z"]] = xyz_data[["X", "Y", "Z"]].astype(
        float
    )  # Convert coordinates to float
    return xyz_identifier, xyz_data


def convert_pdb_to_xyz(filepath_pdb):
    """Split a PDB file into multiple XYZ files based on unique bases."""

    # Extract the filename without extension to create a folder
    filename = os.path.splitext(os.path.basename(filepath_pdb))[0]
    output_dir = os.path.join(os.path.dirname(filepath_pdb), filename)
    os.makedirs(output_dir, exist_ok=True)

    # Variables to hold atom data and track the current base
    elements = []
    coordinates = []
    current_base = None
    xyz_files = {}

    # Read the PDB file
    with open(filepath_pdb, "r") as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Parse element symbol and coordinates
                element = line[76:78].strip()  # Extract element symbol
                x = float(line[30:38].strip())  # Extract X coordinate
                y = float(line[38:46].strip())  # Extract Y coordinate
                z = float(line[46:54].strip())  # Extract Z coordinate

                # Parse base identifier
                residue_name = line[17:20].strip()
                residue_sequence_number = line[22:26].strip()
                base_identifier = (
                    str(residue_sequence_number).zfill(2) + residue_name[1]
                )

                # If base changes, write the previous base to a new .xyz file
                if current_base is not None and base_identifier != current_base:
                    write_xyz_file(output_dir, current_base, elements, coordinates)
                    elements = []
                    coordinates = []

                # Add current atom data
                elements.append(element)
                coordinates.append((x, y, z))
                current_base = base_identifier

        # Write the last base to a .xyz file
        if current_base is not None:
            write_xyz_file(output_dir, current_base, elements, coordinates)

        if DEFAULTS["verbose"]:
            print(f"XYZ files created in directory: {output_dir}")


def write_xyz_file(output_dir, base_identifier, elements, coordinates):
    """Write a single .xyz file for a specific base."""

    filepath_xyz = os.path.join(output_dir, f"{base_identifier}.xyz")
    num_atoms = len(elements)

    # Prepare the content for the XYZ file
    xyz_content = f"{num_atoms}\n{base_identifier}\n"
    for element, (x, y, z) in zip(elements, coordinates):
        xyz_content += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"

    # Write the content to the file
    with open(filepath_xyz, "w") as file:
        file.write(xyz_content)


def convert_json_to_xyz(filename, directory):
    """Converts a JSON file from Pubchem https://pubchem.ncbi.nlm.nih.gov/ to an XYZ file."""

    # Input file path
    filepath_json = os.path.join(directory, filename + ".json")

    # Read the JSON file
    json_file = load_json(filepath_json)

    num_atoms = len(json_file["PC_Compounds"][0]["atoms"]["aid"])
    elements = json_file["PC_Compounds"][0]["atoms"]["element"]
    coordinates = json_file["PC_Compounds"][0]["coords"][0]["conformers"][0]
    x_coords = coordinates["x"]
    y_coords = coordinates["y"]
    z_coords = coordinates["z"]
    atomic_symbols = {1: "H", 6: "C", 7: "N", 8: "O"}

    # Remove the JSON file
    os.remove(filepath_json)

    # Prepare the content for the XYZ file
    xyz_content = f"{num_atoms}\n{filename}\n"  # Header
    for i in range(num_atoms):
        element_symbol = atomic_symbols.get(
            elements[i], "X"
        )  # Default to 'X' if not found
        x = x_coords[i]
        y = y_coords[i]
        z = z_coords[i]
        xyz_content += f"{element_symbol} {x:.4f} {y:.4f} {z:.4f}\n"

    # Output file path
    filepath_xyz = os.path.join(directory, filename + ".xyz")

    # Write the XYZ file
    with open(filepath_xyz, "w") as file:
        file.write(xyz_content)

    if DEFAULTS["verbose"]:
        print(f"File successfully converted and saved at {filepath_xyz}")
