import os
import json
from itertools import combinations

import numpy as np

from .slater_koster import calc_orbital_overlap

# ----------------------------------------------


class BasePair:
    """
    A class to represent a base pair in a DNA molecule.

    Attributes
    ----------
    base1 : object
        The first base in the pair.
    base2 : object
        The second base in the pair.
    identifier : str, optional
        An identifier for the base pair (default is None).
    atoms : list
        List of atoms in the base pair.
    num_atoms : int
        Total number of atoms in the base pair.
    atom_identifiers : list
        List of atom identifiers in the base pair.
    atom_coordinates : ndarray
        Coordinates of atoms in the base pair.
    atom_distance_matrix : ndarray
        Matrix of distances between atoms in the base pair.
    orbitals : list
        List of orbitals in the base pair.
    num_orbitals : int
        Total number of orbitals in the base pair.
    orbital_identifiers : list
        List of orbital identifiers in the base pair.
    orbital_coordinates : ndarray
        Coordinates of orbitals in the base pair.
    orbital_distance_matrix : ndarray
        Matrix of distances between orbitals in the base pair.
    mass : float
        Mass of the base pair.
    center_of_mass : ndarray
        Center of mass of the base pair.
    relative_atom_coordinates : ndarray
        Atom coordinates relative to the center of mass.
    relative_orbital_coordinates : ndarray
        Orbital coordinates relative to the center of mass.
    num_electrons : int
        Total number of valence electrons in the base pair.
    H : ndarray
        LCAO Hamiltonian matrix for the base pair.
    eigv : ndarray
        Eigenvalues of the Hamiltonian matrix.
    eigs : ndarray
        Eigenvectors of the Hamiltonian matrix.
    E_HOMO : float
        Energy of the highest occupied molecular orbital (HOMO).
    HOMO : ndarray
        Highest occupied molecular orbital (HOMO).
    HOMO_occupation : ndarray
        Occupation of the HOMO.
    HOMO_type : str
        Type of the HOMO (sigma, pi, or n).
    HOMO_type_occupation : list
        Occupation of the HOMO by type.
    E_HOMO_1 : float
        Energy of the HOMO-1.
    HOMO_1 : ndarray
        HOMO-1.
    HOMO_occupation_1 : ndarray
        Occupation of the HOMO-1.
    HOMO_type_1 : str
        Type of the HOMO-1 (sigma, pi, or n).
    HOMO_type_occupation_1 : list
        Occupation of the HOMO-1 by type.
    E_HOMO_2 : float
        Energy of the HOMO-2.
    HOMO_2 : ndarray
        HOMO-2.
    HOMO_occupation_2 : ndarray
        Occupation of the HOMO-2.
    HOMO_type_2 : str
        Type of the HOMO-2 (sigma, pi, or n).
    HOMO_type_occupation_2 : list
        Occupation of the HOMO-2 by type.
    E_LUMO : float
        Energy of the lowest unoccupied molecular orbital (LUMO).
    LUMO : ndarray
        Lowest unoccupied molecular orbital (LUMO).
    LUMO_occupation : ndarray
        Occupation of the LUMO.
    LUMO_type : str
        Type of the LUMO (sigma, pi, or n).
    LUMO_type_occupation : list
        Occupation of the LUMO by type.
    E_LUMO_1 : float
        Energy of the LUMO+1.
    LUMO_1 : ndarray
        LUMO+1.
    LUMO_occupation_1 : ndarray
        Occupation of the LUMO+1.
    LUMO_type_1 : str
        Type of the LUMO+1 (sigma, pi, or n).
    LUMO_type_occupation_1 : list
        Occupation of the LUMO+1 by type.
    E_LUMO_2 : float
        Energy of the LUMO+2.
    LUMO_2 : ndarray
        LUMO+2.
    LUMO_occupation_2 : ndarray
        Occupation of the LUMO+2.
    LUMO_type_2 : str
        Type of the LUMO+2 (sigma, pi, or n).
    LUMO_type_occupation_2 : list
        Occupation of the LUMO+2 by type.
    MO_types : dict
        Dictionary of molecular orbital types.
    E_exc : float
        Excitation energy (HOMO-LUMO gap).

    Methods
    -------
    save_results(directory="results"):
        Saves the results to a JSON file.
    calc_MO(identifier, deviation=0):
        Returns properties of the molecular orbital (MO).
    """

    def __init__(self, base1, base2, identifier=None):
        self.base1 = base1
        self.base2 = base2

        # Identifier
        self.identifier = identifier
        if self.identifier is None:
            self.identifier = base1.identifier  # + "_" + base2.identifier

        # Atoms
        self.atoms = self.base1.atoms + self.base2.atoms
        self.num_atoms = self.base1.num_atoms + self.base2.num_atoms
        self.atom_identifiers = (
            self.base1.atom_identifiers + self.base2.atom_identifiers
        )
        self.atom_coordinates = np.block(
            [[self.base1.atom_coordinates], [self.base2.atom_coordinates]]
        )
        self.atom_distance_matrix = self._get_atom_distance_matrix()

        # Orbitals
        self.orbitals = self.base1.orbitals + self.base2.orbitals
        self.num_orbitals = self.base1.num_orbitals + self.base2.num_orbitals
        self.orbital_identifiers = (
            self.base1.orbital_identifiers + self.base2.orbital_identifiers
        )
        self.orbital_coordinates = np.block(
            [[self.base1.orbital_coordinates], [self.base2.orbital_coordinates]]
        )
        self.orbital_distance_matrix = self._get_orbital_distance_matrix()

        # Coordinates relative to the center of mass
        self.mass, self.center_of_mass = self._calc_center_of_mass()
        self.relative_atom_coordinates = np.array(
            [coordinates - self.center_of_mass for coordinates in self.atom_coordinates]
        )
        self.relative_orbital_coordinates = np.array(
            [
                coordinates - self.center_of_mass
                for coordinates in self.orbital_coordinates
            ]
        )

        # Valence electrons
        self.num_electrons = self.base1.num_electrons + self.base2.num_electrons

        # Molecular Orbitals
        # ------------------

        # LCAO Hamiltonian
        self.H = self._calc_H()

        # Solve the eigenvalue problem
        self.eigv, self.eigs = np.linalg.eigh(self.H)

        # HOMO and LUMO
        (
            self.E_HOMO,
            self.HOMO,
            self.HOMO_occupation,
            self.HOMO_type,
            self.HOMO_type_occupation,
        ) = self.calc_MO("HOMO")
        (
            self.E_HOMO_1,
            self.HOMO_1,
            self.HOMO_occupation_1,
            self.HOMO_type_1,
            self.HOMO_type_occupation_1,
        ) = self.calc_MO("HOMO", deviation=1)
        (
            self.E_HOMO_2,
            self.HOMO_2,
            self.HOMO_occupation_2,
            self.HOMO_type_2,
            self.HOMO_type_occupation_2,
        ) = self.calc_MO("HOMO", deviation=2)
        (
            self.E_LUMO,
            self.LUMO,
            self.LUMO_occupation,
            self.LUMO_type,
            self.LUMO_type_occupation,
        ) = self.calc_MO("LUMO")
        (
            self.E_LUMO_1,
            self.LUMO_1,
            self.LUMO_occupation_1,
            self.LUMO_type_1,
            self.LUMO_type_occupation_1,
        ) = self.calc_MO("LUMO", deviation=1)
        (
            self.E_LUMO_2,
            self.LUMO_2,
            self.LUMO_occupation_2,
            self.LUMO_type_2,
            self.LUMO_type_occupation_2,
        ) = self.calc_MO("LUMO", deviation=2)

        self.MO_types = {
            "LUMO_2": self.LUMO_type_2,
            "LUMO_1": self.LUMO_type_1,
            "LUMO": self.LUMO_type,
            "HOMO": self.HOMO_type,
            "HOMO_1": self.HOMO_type_1,
            "HOMO_2": self.HOMO_type_2,
        }

        # Excitation energy (HOMO-LUMO gap)
        self.E_exc = self.E_LUMO - self.E_HOMO

    # ---------------------------------------------------------

    def __repr__(self):
        return f'BasePair({self.base1}, {self.base2}, identifier = "{self.identifier}")'

    def save_results(self, directory="results"):
        """
        Save the calculation results to a JSON file.

        Parameters
        ----------
        directory : str, optional
            The directory where the results file will be saved. Default is "results".

        Notes
        -----
        The results are saved in a JSON file named "results_<identifier>.json" where <identifier>
        is a unique identifier for the calculation. The file contains the following keys:
        - "E_HOMO": The energy of the highest occupied molecular orbital (HOMO), rounded to 4 decimal places.
        - "E_LUMO": The energy of the lowest unoccupied molecular orbital (LUMO), rounded to 4 decimal places.
        - "HOMO": The list of HOMO values.
        - "LUMO": The list of LUMO values.
        """

        dictionary = {
            "E_HOMO": round(float(self.E_HOMO), 4),
            "E_LUMO": round(float(self.E_LUMO), 4),
            "HOMO": list(self.HOMO),
            "LUMO": list(self.LUMO),
        }

        filename = "results_" + self.identifier + ".json"
        filepath = os.path.join(directory, filename)

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(dictionary, f, indent=2)
        print("Results saved at" + filepath)

    def calc_MO(self, identifier, deviation=0):
        """
        Calculate the Molecular Orbital (MO) properties.

        Parameters
        ----------
        identifier : str
            The identifier for the MO to calculate. Should be either "HOMO" (Highest Occupied Molecular Orbital)
            or "LUMO" (Lowest Unoccupied Molecular Orbital).
        deviation : int, optional
            The deviation from the HOMO or LUMO level to calculate. Default is 0.

        Returns
        -------
        E_MO : float
            The energy of the specified MO.
        MO : numpy.ndarray
            The molecular orbital eigenvector.
        MO_occupation : numpy.ndarray
            The occupation of the molecular orbital.
        MO_type : str
            The type of the molecular orbital.
        MO_type_occupation : float
            The occupation of the molecular orbital type.
        """

        start = 0
        if identifier == "HOMO":
            start = -1 - deviation
        if identifier == "LUMO":
            start = deviation

        E_MO = self.eigv[self.num_electrons // 2 + start]
        MO = self.eigs[:, self.num_electrons // 2 + start]
        MO_occupation = MO.conj() * MO
        MO_type, MO_type_occupation = self._get_MO_type(MO_occupation)

        return E_MO, MO, MO_occupation, MO_type, MO_type_occupation

    # --------------------------------------------------------

    def _get_atom_distance_matrix(self):
        """Returns a matrix of distances between atoms in the molecule."""
        atom_distance_matrix = np.zeros((self.num_atoms, self.num_atoms))
        for i, j in combinations(range(self.num_atoms), r=2):
            value = np.linalg.norm(self.atom_coordinates[i] - self.atom_coordinates[j])
            atom_distance_matrix[i, j] = value
            atom_distance_matrix[j, i] = value
        return atom_distance_matrix

    def _get_orbital_distance_matrix(self):
        """Returns a matrix of distances between orbitals in the molecule."""
        orbital_distance_matrix = np.zeros((self.num_orbitals, self.num_orbitals))
        for i, j in combinations(range(self.num_orbitals), r=2):
            value = np.linalg.norm(
                self.orbital_coordinates[i] - self.orbital_coordinates[j]
            )
            orbital_distance_matrix[i, j] = value
            orbital_distance_matrix[j, i] = value
        return orbital_distance_matrix

    def _calc_center_of_mass(self):
        """Calculates the center of mass of the molecule."""
        molecule_mass = self.base1.mass + self.base2.mass
        center_of_mass = (
            self.base1.mass * self.base1.center_of_mass
            + self.base2.mass * self.base2.center_of_mass
        )
        center_of_mass /= molecule_mass
        return molecule_mass, center_of_mass

    def _calc_H(self):
        """Calculates the LCAO Hamiltonian matrix for the molecule."""
        H_int = np.zeros((self.base1.num_orbitals, self.base2.num_orbitals))
        for i, orbital1 in enumerate(self.base1.orbitals):
            for j, orbital2 in enumerate(self.base2.orbitals):
                value = calc_orbital_overlap(orbital1, orbital2, "interbase")
                H_int[i, j] = value
        H = np.block([[self.base1.H, H_int], [H_int.conj().T, self.base2.H]])
        return H

    # Molecular Orbitals
    def _get_MO_type(self, MO_occupation):
        """Returns the nature of the molecular orbitals: sigma, pi, n (non-bonding)."""

        identifiers = [
            identifier.split("_")[0] + "_" + identifier.split("_")[2]
            for identifier in self.orbital_identifiers
        ]
        s = ["C_s", "C_px", "C_py", "N_s", "O_s", "O_px", "O_py", "H_s"]
        s_mask = np.array([int(identifier in s) for identifier in identifiers])
        pi = ["C_pz", "N_pz", "O_pz"]
        pi_mask = np.array([int(identifier in pi) for identifier in identifiers])
        n = ["N_s", "N_px", "N_py"]
        n_mask = np.array([int(identifier in n) for identifier in identifiers])

        s_pop = sum(s_mask * MO_occupation)
        pi_pop = sum(pi_mask * MO_occupation)
        n_pop = sum(n_mask * MO_occupation)

        MO_type = ""
        if s_pop >= max(pi_pop, n_pop):
            MO_type = "sigma"

        if pi_pop >= max(s_pop, n_pop):
            MO_type = "pi"

        if n_pop >= max(s_pop, pi_pop):
            MO_type = "n"

        MO_type_occupation = [s_pop, pi_pop, n_pop]
        return MO_type, MO_type_occupation
