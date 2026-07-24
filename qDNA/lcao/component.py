from itertools import combinations
import numpy as np

from ..io import load_xyz

# ----------------------------------------------------------------------


class Component:
    """
    Represents a component of DNA, e.g. a nucleobase, the sugar-phosphate backbone, etc.

    Attributes
    ----------
    filepath : str
        Path to the XYZ file.
    kwargs : dict
        Keyword arguments for LCAO configuration.
    atoms : list of str
        List of atom types in the molecule.
    atoms_coordinates : ndarray
        Array of atomic coordinates.
    atoms_id : list of str
        Unique identifiers for each atom.
    num_atoms : int
        Total number of atoms in the molecule.
    orbitals : list of str
        List of orbital identifiers.
    orbitals_coordinates : ndarray
        Array of orbital coordinates.
    num_orbitals : int
        Total number of orbitals in the molecule.
    num_electrons : int
        Total number of valence electrons in the molecule.

    Parameters
    ----------
    filepath : str
        Path to the XYZ file containing atomic coordinates.

    Methods
    -------
    get_orbital_distance_matrix()
        Calculate the distance matrix between orbitals.
    get_orbital_bond_matrix(cutoff=1.59)
        Generate a bond matrix based on orbital distances and a cutoff value.

    Notes
    -----
    .. note::
        Coordinates are in Angstroms (Å), not nanometers (nm).

    """

    def __init__(self, filepath):

        # load component (XYZ file)
        self.filepath = filepath
        _, atoms, coordinates = load_xyz(self.filepath)

        # component atoms
        self.atoms = atoms
        self.atoms_coordinates = np.array(coordinates)
        self.atoms_id = [f"{atom}_{atom_idx}" for atom_idx, atom in enumerate(self.atoms)]
        self.num_atoms = len(self.atoms)

        # component orbitals
        self.orbitals, self.orbitals_coordinates = self._get_orbitals()
        self.num_orbitals = len(self.orbitals)

        # valence electrons
        num_atom_electrons = {"H": 1, "C": 4, "N": 5, "O": 6, "P": 5, "X": 0}
        self.num_electrons = sum(num_atom_electrons[atom] for atom in self.atoms)

    def __repr__(self):
        return f"Component({self.filepath})"

    def _get_orbitals(self):
        """Generate orbital identifiers and coordinates based on the atoms in the component."""

        orbitals = []
        orbitals_coordinates = []

        for i in range(self.num_atoms):
            atom = self.atoms[i]
            atom_coordinates = self.atoms_coordinates[i]
            if atom == "H":
                orbital_types = ["s"]
            elif atom in ["C", "N", "O", "P"]:
                orbital_types = ["s", "px", "py", "pz"]
            else:
                orbital_types = []
                raise ValueError(f"Unknown atom type: {atom}")

            for orbital_type in orbital_types:
                orbitals.append(atom + "_" + orbital_type)
                orbitals_coordinates.append(atom_coordinates)
        return orbitals, orbitals_coordinates

    def get_orbital_distance_matrix(self):
        """Calculate the distance matrix between atomic orbitals."""

        orbital_distance_matrix = np.zeros((self.num_orbitals, self.num_orbitals))
        for i, j in combinations(range(self.num_orbitals), r=2):
            vector = self.orbitals_coordinates[i] - self.orbitals_coordinates[j]
            distance = np.linalg.norm(vector)
            orbital_distance_matrix[i, j] = distance
            orbital_distance_matrix[j, i] = distance
        return orbital_distance_matrix

    def get_orbital_bond_matrix(self, cutoff=1.59):
        """Generate a bond matrix based between atomic orbitals. Two atomic orbitals are
        bonded whenever their distance is bigger then zero but smaller than a cutoff value."""

        orbital_distance_matrix = self.get_orbital_distance_matrix()
        orbital_bond_matrix = (orbital_distance_matrix > 0) & (orbital_distance_matrix < cutoff)
        return orbital_bond_matrix.astype(int)


# ----------------------------------------------------------------------
