import os
from .helpers import get_non_overwriting_path

# ----------------------------------------------------------------------


def create_xyz(elements, coordinates, base_id):
    num_atoms = len(elements)
    xyz_content = f"{num_atoms}\n{base_id}\n"
    for element, (x, y, z) in zip(elements, coordinates):
        xyz_content += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"
    return xyz_content


def save_xyz(xyz_content, filepath):
    filepath = get_non_overwriting_path(filepath)
    with open(filepath, "w", encoding="utf-8") as file:
        file.write(xyz_content)


def write_xyz(directory, base_id, elements, coordinates):
    xyz_content = create_xyz(elements, coordinates, base_id)
    filepath = os.path.join(directory, f"{base_id}.xyz")
    save_xyz(xyz_content, filepath)


def load_xyz(filepath):
    with open(filepath, "r", encoding="utf-8") as file:
        lines = file.readlines()

    identifier = lines[1].strip()

    atoms = []
    coordinates = []

    for line in lines[2:]:
        parts = line.split()
        atom = parts[0]
        x, y, z = map(float, parts[1:4])

        atoms.append(atom)
        coordinates.append((x, y, z))

    return identifier, atoms, coordinates


def find_xyz(directory):
    files = os.listdir(directory)
    return [os.path.splitext(file)[0] for file in files if file.endswith(".xyz")]


# ----------------------------------------------------------------------
