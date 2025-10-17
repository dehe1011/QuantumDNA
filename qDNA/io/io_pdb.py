import os

from .io_xyz import write_xyz

# ----------------------------------------------------------------------


def load_pdb(filepath):
    """
    Load atomic data from a PDB file.
    Parses the ATOM and HETATM records in a PDB file and extracts relevant
    information such as atom type, residue, chain, residue ID, coordinates,
    and element type.

    Parameters
    ----------
    filepath : str
        Path to the PDB file to be loaded.

    Returns
    -------
    list of dict
        A list of dictionaries, each containing atomic data with keys:
        'atom', 'residue', 'chain', 'res_id', 'x', 'y', 'z', and 'element'.
    """

    pdb_content = []

    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            # if line.startswith("TER"):
            #     break
            if line.startswith("ATOM") or line.startswith("HETATM"):  # careful with HETATM
                pdb_data = {
                    "atom": line[12:16].strip(),
                    "residue": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "res_id": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": line[76:78].strip(),
                }
                pdb_content.append(pdb_data)

    return pdb_content


def modify_base_idx(base_idx, start_idx, n, **kwargs):

    lower_idx_list = kwargs.get("lower_idx_list", None)
    lower_offset = kwargs.get("lower_offset", 0)
    lower_direction = kwargs.get("lower_direction", "5-3")
    lower_continued = kwargs.get("lower_continued", True)

    if lower_idx_list is not None:
        return lower_idx_list[base_idx]

    base_idx_mod = None
    i = abs(base_idx)
    if lower_direction == "5-3" and lower_continued:
        base_idx_mod = i
    if lower_direction == "5-3" and not lower_continued:
        base_idx_mod = n + i
    if lower_direction == "3-5" and lower_continued:
        base_idx_mod = 2 * n - 1 - (i - start_idx) + n + start_idx
    if lower_direction == "3-5" and not lower_continued:
        base_idx_mod = 2 * n - 1 - (i - start_idx) + start_idx
    base_idx_mod += lower_offset
    return base_idx_mod


def pdb_to_xyz(filepath, **kwargs):
    """
    Converts a PDB file to XYZ format and writes the output to separate files
    for each base and backbone in the structure.

    Parameters
    ----------
    filepath : str
        Path to the input PDB file.

    Notes
    -----
    .. note::
        - The function creates a directory named after the input file (without extension) to store the output XYZ files.
        - Each base and backbone is written to separate XYZ files.
        - Base indices are adjusted for the lower strand if applicable.

    """

    filename = os.path.splitext(os.path.basename(filepath))[0]
    directory = os.path.join(os.path.dirname(filepath), filename)
    os.makedirs(directory, exist_ok=True)

    elements, elements_backbone = [], []
    coordinates, coordinates_backbone = [], []

    start_idx = 0
    old_base_id = None
    old_chain_id = None
    old_base_idx = 0

    lower_strand = False
    n = 0

    pdb_content = load_pdb(filepath)

    for i, entry in enumerate(pdb_content):

        element_id = entry["atom"]
        chain_id = entry["chain"]  # e.g. A
        base_idx = entry["res_id"]  # e.g. 1
        x, y, z = entry["x"], entry["y"], entry["z"]
        element = entry["element"]  # e.g. N, C, O

        if element == "":
            element = element_id[0]

        # increase base_idx for the lower strand
        first_entry = i == 0
        if first_entry:
            start_idx = base_idx

        if kwargs.get('no_chain_id', False):
            chain_changes = base_idx == 1 and old_base_idx == 27
        else:
            chain_changes = chain_id != old_chain_id  #  # changed!!
        if chain_changes and not first_entry:
            lower_strand = True
            n = old_base_idx + 1 - start_idx

        if lower_strand:
            base_idx = modify_base_idx(base_idx, start_idx, n, **kwargs)

        base_id = entry["residue"]  # e.g. DC
        base_id = str(base_idx).zfill(2) + base_id[1]  # e.g. 01C
        backbone_id = str(base_idx).zfill(2) + "B"  # e.g. 01B

        info = None  # for debugging

        if base_id != old_base_id and old_base_id is not None:
            write_xyz(directory, old_base_id, elements, coordinates, info=info)
            write_xyz(
                directory, old_backbone_id, elements_backbone, coordinates_backbone, info=info
            )
            elements, elements_backbone = [], []
            coordinates, coordinates_backbone = [], []

        if "'" in element_id or "P" in element_id:
            elements_backbone.append(element)
            coordinates_backbone.append((x, y, z))
        else:
            elements.append(element)
            coordinates.append((x, y, z))

        old_base_idx = base_idx
        old_base_id = base_id
        old_chain_id = chain_id
        old_backbone_id = backbone_id

    if old_base_id is not None:
        write_xyz(directory, old_base_id, elements, coordinates, info=info)
        write_xyz(directory, old_backbone_id, elements_backbone, coordinates_backbone, info=info)


def find_pdb(directory):
    files = os.listdir(directory)
    return [os.path.splitext(file)[0] for file in files if file.endswith(".pdb")]


# ----------------------------------------------------------------------
