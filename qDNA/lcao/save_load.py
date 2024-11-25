import json
import os

import pandas as pd

from .. import DATA_DIR, DEFAULTS

__all__ = [
    "load_xyz",
    "PARAMETRIZATION",
    "find_xyz_files",
    "lcao_load_json",
    "convert_json_to_xyz",
    "convert_pdb_to_xyz",
]

# ----------------------------- JSON -----------------------------


def lcao_load_json(filepath):
    """Loads a JSON file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        message = f"File {filepath} does not exist."
        if DEFAULTS["verbose"]:
            print(message)
    return data


PARAMETRIZATION = lcao_load_json(
    os.path.join(
        DATA_DIR,
        "raw",
        "lcao_params",
        DEFAULTS["lcao_default"]["parametrization"] + ".json",
    )
)

# ----------------------------- XYZ -----------------------------


def find_xyz_files(directory):
    """
    Find all .xyz files in the given directory.
    Parameters
    ----------
    directory : str
        The path to the directory where .xyz files are to be searched.
    Returns
    -------
    list of str
        A list of filenames (without extensions) of all .xyz files found in the directory.
    Raises
    ------
    Exception
        If an error occurs while accessing the directory or reading its contents.
    """

    try:
        files = os.listdir(directory)
        xyz_files = [
            os.path.splitext(file)[0] for file in files if file.endswith(".xyz")
        ]
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    return xyz_files


def load_xyz(filename, directory=os.path.join(DATA_DIR, "geometries")):
    """
    Load atomic coordinates from an XYZ file.

    Parameters
    ----------
    filename : str
        The name of the XYZ file (without the .xyz extension).
    directory : str, optional
        The directory where the XYZ file is located. Default is a subdirectory "geometries" within DATA_DIR.

    Returns
    -------
    tuple
        A tuple containing:
        - xyz_identifier : str
            The identifier from the second line of the XYZ file.
        - xyz_data : pandas.DataFrame
            A DataFrame containing the atomic coordinates with columns ["Atom", "X", "Y", "Z"].

    Notes
    -----
    The XYZ file is expected to have the following format:
    - The first line contains the number of atoms (ignored).
    - The second line contains a comment or identifier.
    - Subsequent lines contain atomic coordinates in the format: Atom X Y Z.
    """

    filepath = os.path.join(directory, filename + ".xyz")

    with open(filepath, "r", encoding="utf-8") as file:
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
    """
    Converts a PDB file to multiple XYZ files, one for each base.

    Parameters
    ----------
    filepath_pdb : str
        The path to the input PDB file.

    Notes
    -----
    - The function creates a directory named after the input PDB file (without extension)
      to store the generated XYZ files.
    - Each base in the PDB file is written to a separate XYZ file.
    - If the chain identifier changes, the base numbering is adjusted to continue from
      the previous chain's last base number.
    - The function assumes that the base counter starts from one if the chain changes.
    - The function prints the directory where the XYZ files are created if the verbose
      mode is enabled in the DEFAULTS dictionary.
    """

    # Extract the filename without extension to create a folder
    filename = os.path.splitext(os.path.basename(filepath_pdb))[0]
    output_dir = os.path.join(os.path.dirname(filepath_pdb), filename)
    os.makedirs(output_dir, exist_ok=True)

    # Variables to hold atom data and track the current base
    elements = []
    coordinates = []
    current_base = None

    current_chain = ""
    current_base_number = 0

    # Read the PDB file
    lower_strand = False
    start_from_one = False
    num_bases_per_strand = 0

    with open(filepath_pdb, "r", encoding="utf-8") as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Parse element symbol and coordinates
                element = line[76:78].strip()  # Extract element symbol
                x = float(line[30:38].strip())  # Extract X coordinate
                y = float(line[38:46].strip())  # Extract Y coordinate
                z = float(line[46:54].strip())  # Extract Z coordinate

                # Parse base identifier
                base_identifier = line[17:20].strip()
                chain_identifier = line[21].strip()
                base_number = int(line[22:26].strip())

                chain_changes = current_chain != chain_identifier
                base_changes = current_base_number != base_number

                # checks if the base counter starts from one if the chain changes
                if chain_changes:
                    lower_strand = True
                    num_bases_per_strand = current_base_number
                    start_from_one = bool(base_number == 1)

                if base_changes and lower_strand and start_from_one:
                    base_number += num_bases_per_strand

                current_chain = chain_identifier
                current_base_number = base_number

                base_identifier = str(base_number).zfill(2) + base_identifier[1]

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
    """
    Write atomic coordinates to an XYZ file.

    Parameters
    ----------
    output_dir : str
        The directory where the XYZ file will be saved.
    base_identifier : str
        The base name for the XYZ file.
    elements : list of str
        A list of atomic element symbols.
    coordinates : list of tuple of float
        A list of tuples, each containing the x, y, and z coordinates of an atom.

    Returns
    -------
    None
    """

    filepath_xyz = os.path.join(output_dir, f"{base_identifier}.xyz")
    num_atoms = len(elements)

    # Prepare the content for the XYZ file
    xyz_content = f"{num_atoms}\n{base_identifier}\n"
    for element, (x, y, z) in zip(elements, coordinates):
        xyz_content += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"

    # Write the content to the file
    with open(filepath_xyz, "w", encoding="utf-8") as file:
        file.write(xyz_content)


def convert_json_to_xyz(filename, directory):
    """
    Converts a JSON file to an XYZ file, e.g., from Pubchem https://pubchem.ncbi.nlm.nih.gov/.

    Parameters
    ----------
    filename : str
        The name of the JSON file (without extension) to be converted.
    directory : str
        The directory where the JSON file is located and where the XYZ file will be saved.

    Raises
    ------
    FileNotFoundError
        If the JSON file does not exist in the specified directory.
    KeyError
        If the expected keys are not found in the JSON file.
    ValueError
        If the JSON file contains invalid data.

    Notes
    -----
    The JSON file is expected to follow a specific structure with atomic data under
    "PC_Compounds" -> "atoms" and coordinates under "coords" -> "conformers".
    The atomic symbols are mapped from atomic numbers using a predefined dictionary.
    The JSON file is removed after conversion.

    Examples
    --------
    >>> convert_json_to_xyz("molecule", "/path/to/directory")
    File successfully converted and saved at /path/to/directory/molecule.xyz
    """

    # Input file path
    filepath_json = os.path.join(directory, filename + ".json")

    # Read the JSON file
    json_file = lcao_load_json(filepath_json)

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
    with open(filepath_xyz, "w", encoding="utf-8") as file:
        file.write(xyz_content)

    if DEFAULTS["verbose"]:
        print(f"File successfully converted and saved at {filepath_xyz}")
