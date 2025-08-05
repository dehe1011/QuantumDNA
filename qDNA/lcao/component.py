import copy
from itertools import combinations
import numpy as np

from ..utils import check_lcao_kwargs
from ..io import load_xyz, load_lcao_param, DEFAULTS

# ----------------------------------------------------------------------


class Component:
    """
    Represents a component of DNA, e.g. a nucleobase, the sugar-phosphate backbone, etc.

    Attributes
    ----------
    param_id : str
        Identifier for the LCAO parameters.
    lcao_param : dict
        Loaded LCAO parameters.
    num_atom_orbitals : dict
        Number of orbitals for each atom type.
    num_atom_electrons : dict
        Number of electrons for each atom type.
    filepath : str
        Path to the XYZ file.
    xyz_id : str
        Identifier for the XYZ file.
    atoms : list of str
        List of atom types in the molecule.
    atoms_coordinates : ndarray
        Array of atomic coordinates.
    atoms_id : list of str
        Unique identifiers for each atom.
    num_atoms : int
        Total number of atoms in the molecule.
    num_electrons : int
        Total number of electrons in the molecule.
    orbitals : list of str
        List of orbital identifiers.
    orbitals_coordinates : ndarray
        Array of orbital coordinates.
    num_orbitals : int
        Total number of orbitals in the molecule.

    Methods
    -------
    get_orbital_distance_matrix()
        Calculates the distance matrix between orbitals.
    get_orbital_bond_matrix(cutoff=1.59)
        Generates a bond matrix based on orbital distances and a cutoff value.

    Notes
    -----
    .. note::
        The coordinates are given in Angstroms (Å), not in nanometers (nm).

    """

    def __init__(self, filepath, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(DEFAULTS["lcao_kwargs_default"])
        self.kwargs.update(kwargs)
        check_lcao_kwargs(**self.kwargs)

        self.param_id = self.kwargs.get("param_id")
        self.lcao_param = load_lcao_param(self.param_id)

        self.num_atom_orbitals = {"H": 1, "C": 4, "N": 4, "O": 4, "P": 4, "X": 0}
        self.num_atom_electrons = {"H": 1, "C": 4, "N": 5, "O": 6, "P": 5, "X": 0}

        self.filepath = filepath
        self.xyz_id, atoms, coordinates = load_xyz(self.filepath)

        self.atoms = atoms
        self.atoms_coordinates = np.array(coordinates)
        self.atoms_id = [f"{atom}_{atom_idx}" for atom_idx, atom in enumerate(self.atoms)]
        self.num_atoms = len(self.atoms)
        self.num_electrons = sum(self.num_atom_electrons[atom] for atom in self.atoms)

        self.orbitals, self.orbitals_coordinates = self._get_orbitals()
        self.num_orbitals = len(self.orbitals)

    def _get_orbitals(self):
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
        orbital_distance_matrix = np.zeros((self.num_orbitals, self.num_orbitals))
        for i, j in combinations(range(self.num_orbitals), r=2):
            vector = self.orbitals_coordinates[i] - self.orbitals_coordinates[j]
            distance = np.linalg.norm(vector)
            orbital_distance_matrix[i, j] = distance
            orbital_distance_matrix[j, i] = distance
        return orbital_distance_matrix

    def get_orbital_bond_matrix(self, cutoff=1.59):
        orbital_distance_matrix = self.get_orbital_distance_matrix()
        orbital_bond_matrix = (orbital_distance_matrix > 0) & (orbital_distance_matrix < cutoff)
        return orbital_bond_matrix.astype(int)


# ----------------------------------------------------------------------
