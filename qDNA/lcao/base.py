import os
import json
from itertools import combinations

import numpy as np
import scipy.constants as c

from .slater_koster import calc_orbital_energy, calc_orbital_overlap

# --------------------------------------------------


class Orbital:  # pylint: disable=too-few-public-methods
    """
    A class used to represent an atomic orbital.

    Parameters
    ----------
    coordinates : tuple
        A tuple representing the coordinates of the orbital.
    atom_identifier : str
        A string identifier for the atom, in the format "atom_atom_idx".
    orbital_type : str
        A string representing the type of the orbital.

    Attributes
    ----------
    coordinates : tuple
        The coordinates of the orbital.
    atom_identifier : str
        The identifier for the atom.
    atom : str
        The atom part of the identifier.
    atom_idx : str
        The index part of the identifier.
    orbital_type : str
        The type of the orbital.
    identifier : str
        A unique identifier for the orbital, combining atom_identifier and orbital_type.
    """

    def __init__(self, coordinates, atom_identifier, orbital_type):
        self.coordinates = coordinates
        self.atom_identifier = atom_identifier
        self.atom, self.atom_idx = atom_identifier.split("_")
        self.orbital_type = orbital_type
        self.identifier = f"{atom_identifier}_{orbital_type}"

    def __repr__(self):
        return (
            f"Orbital({self.coordinates}, {self.atom_identifier}, {self.orbital_type})"
        )


class Base:
    """
    A class to represent a molecular system and perform various calculations on it.

    Parameters
    ----------
    xyz_identifier : str
        Identifier for the xyz data.
    xyz_data : pandas.DataFrame
        DataFrame containing the xyz coordinates and atom types.

    Attributes
    ----------
    xyz_data : pandas.DataFrame
        DataFrame containing the xyz coordinates and atom types.
    identifier : str
        Identifier for the xyz data.
    atoms : list of str
        List of atoms in the molecule.
    num_atoms : int
        Number of atoms in the molecule.
    atom_identifiers : list of str
        List of atom identifiers (e.g., 'C_0', 'H_1').
    atom_coordinates : numpy.ndarray
        Array of atom coordinates.
    atom_distance_matrix : numpy.ndarray
        Matrix of distances between atoms.
    atom_bond_matrix : numpy.ndarray
        Matrix indicating bonds between atoms.
    num_atom_orbitals : dict
        Dictionary mapping atom types to the number of orbitals.
    orbitals : list of Orbital
        List of orbitals in the molecule.
    num_orbitals : int
        Number of orbitals in the molecule.
    orbital_identifiers : list of str
        List of orbital identifiers (e.g., 'C_0_pz', 'H_1_s').
    orbital_coordinates : numpy.ndarray
        Array of orbital coordinates.
    orbital_distance_matrix : numpy.ndarray
        Matrix of distances between orbitals.
    orbital_bond_matrix : numpy.ndarray
        Matrix indicating bonds between orbitals.
    atom_masses : dict
        Dictionary mapping atom types to their masses.
    mass : float
        Total mass of the molecule.
    center_of_mass : numpy.ndarray
        Coordinates of the center of mass of the molecule.
    relative_atom_coordinates : numpy.ndarray
        Atom coordinates relative to the center of mass.
    relative_orbital_coordinates : numpy.ndarray
        Orbital coordinates relative to the center of mass.
    num_atom_electrons : dict
        Dictionary mapping atom types to the number of valence electrons.
    num_electrons : int
        Number of valence electrons in the molecule.
    H : numpy.ndarray
        LCAO Hamiltonian matrix.
    eigv : numpy.ndarray
        Eigenvalues of the Hamiltonian matrix.
    eigs : numpy.ndarray
        Eigenvectors of the Hamiltonian matrix.
    E_HOMO : float
        Energy of the highest occupied molecular orbital (HOMO).
    HOMO : numpy.ndarray
        Highest occupied molecular orbital (HOMO).
    HOMO_occupation : numpy.ndarray
        Occupation of the HOMO.
    HOMO_type : str
        Type of the HOMO (sigma, pi, n).
    HOMO_type_occupation : list of float
        Occupation of the HOMO by type.
    E_HOMO_1 : float
        Energy of the HOMO-1.
    HOMO_1 : numpy.ndarray
        HOMO-1.
    HOMO_occupation_1 : numpy.ndarray
        Occupation of the HOMO-1.
    HOMO_type_1 : str
        Type of the HOMO-1 (sigma, pi, n).
    HOMO_type_occupation_1 : list of float
        Occupation of the HOMO-1 by type.
    E_HOMO_2 : float
        Energy of the HOMO-2.
    HOMO_2 : numpy.ndarray
        HOMO-2.
    HOMO_occupation_2 : numpy.ndarray
        Occupation of the HOMO-2.
    HOMO_type_2 : str
        Type of the HOMO-2 (sigma, pi, n).
    HOMO_type_occupation_2 : list of float
        Occupation of the HOMO-2 by type.
    E_LUMO : float
        Energy of the lowest unoccupied molecular orbital (LUMO).
    LUMO : numpy.ndarray
        Lowest unoccupied molecular orbital (LUMO).
    LUMO_occupation : numpy.ndarray
        Occupation of the LUMO.
    LUMO_type : str
        Type of the LUMO (sigma, pi, n).
    LUMO_type_occupation : list of float
        Occupation of the LUMO by type.
    E_LUMO_1 : float
        Energy of the LUMO+1.
    LUMO_1 : numpy.ndarray
        LUMO+1.
    LUMO_occupation_1 : numpy.ndarray
        Occupation of the LUMO+1.
    LUMO_type_1 : str
        Type of the LUMO+1 (sigma, pi, n).
    LUMO_type_occupation_1 : list of float
        Occupation of the LUMO+1 by type.
    E_LUMO_2 : float
        Energy of the LUMO+2.
    LUMO_2 : numpy.ndarray
        LUMO+2.
    LUMO_occupation_2 : numpy.ndarray
        Occupation of the LUMO+2.
    LUMO_type_2 : str
        Type of the LUMO+2 (sigma, pi, n).
    LUMO_type_occupation_2 : list of float
        Occupation of the LUMO+2 by type.
    MO_types : dict
        Dictionary mapping molecular orbital types to their identifiers.
    E_exc : float
        Excitation energy (HOMO-LUMO gap).
    dipole_moment : float
        Dipole moment of the molecule.
    oscillator_strength : float
        Oscillator strength of the molecule.

    Methods
    -------
    save_results(directory='results')
        Saves the results to a JSON file.
    calc_dipole_moment(MO_1=None, MO_2=None, unit='Coulomb*Angstrom')
        Calculates the dipole moment between two molecular orbitals.
    calc_oscillator_strength(dipole_moment)
        Calculates the dimensionless oscillator strength.
    calc_MO(identifier, deviation=0)
        Returns properties of the molecular orbital.
    """

    def __init__(self, xyz_identifier, xyz_data):
        # load the xyz file
        self.xyz_data = xyz_data

        # Identifier
        self.identifier = xyz_identifier

        # Atoms, Orbitals, Center of Mass, Valence Electrons
        # --------------------------------------------------

        # Atoms
        self.atoms = self._get_atoms()
        self.num_atoms = len(self.atoms)
        self.atom_identifiers = [
            f"{atom}_{atom_idx}" for atom_idx, atom in enumerate(self.atoms)
        ]  # e.g. C_0, H_1, O_2
        for atom_id in self.atom_identifiers:
            assert atom_id[0] in [
                "H",
                "C",
                "N",
                "O",
            ], "Your File contains atoms other than ['H', 'C', 'N', 'O']. Maybe you forgot to remove the sugar-phosphate backbone?"
        self.atom_coordinates = self._get_atom_coordinates()
        self.atom_distance_matrix = self._get_atom_distance_matrix()
        self.atom_bond_matrix = (
            (self.atom_distance_matrix > 0) & (self.atom_distance_matrix < 1.59)
        ).astype(int)

        # Orbitals
        self.num_atom_orbitals = {"H": 1, "C": 4, "N": 4, "O": 4}
        self.orbitals = self._get_orbitals()
        self.num_orbitals = self._get_num_orbitals()
        self.orbital_identifiers = [
            orbital.identifier for orbital in self.orbitals
        ]  # e.g. C_0_pz, H_1_s, O_2_px
        self.orbital_coordinates = np.array(
            [orbital.coordinates for orbital in self.orbitals]
        )
        self.orbital_distance_matrix = self._get_orbital_distance_matrix()
        self.orbital_bond_matrix = (
            (self.orbital_distance_matrix > 0) & (self.orbital_distance_matrix < 1.59)
        ).astype(int)

        # Coordinates relative to the center of mass
        self.atom_masses = {"H": 1, "C": 6, "N": 7, "O": 8}
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
        self.num_atom_electrons = {"H": 1, "C": 4, "N": 5, "O": 6}
        self.num_electrons = self._get_num_electrons()

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

        # Dipole moment and oscillator strength
        # -------------------------------------

        self.dipole_moment = self.calc_dipole_moment()  # in Coulomb*Angstrom
        self.oscillator_strength = self.calc_oscillator_strength(
            self.dipole_moment
        )  # no unit

    # ---------------------------------------------------------

    def __repr__(self):
        return f'Base("{self.identifier}")'

    def save_results(self, directory="results"):
        """
        Save the calculation results to a JSON file.

        Parameters
        ----------
        directory : str, optional
            The directory where the results file will be saved. Default is "results".

        Notes
        -----
        .. note::

            The results are saved in a JSON file named "results_<identifier>.json" where <identifier>
            is a unique identifier for the calculation. The file contains the following keys:

            - "E_HOMO": The energy of the highest occupied molecular orbital (HOMO), rounded to 4 decimal places.
            - "E_LUMO": The energy of the lowest unoccupied molecular orbital (LUMO), rounded to 4 decimal places.
            - "HOMO": The list of HOMO values.
            - "LUMO": The list of LUMO values.
            - "dipole_moment": The dipole moment, rounded to 4 decimal places.
            - "oscillator_strength": The oscillator strength, rounded to 4 decimal places.
        """

        dictionary = {
            "E_HOMO": round(float(self.E_HOMO), 4),
            "E_LUMO": round(float(self.E_LUMO), 4),
            "HOMO": list(self.HOMO),
            "LUMO": list(self.LUMO),
            "dipole_moment": round(float(self.dipole_moment), 4),
            "oscillator_strength": round(float(self.oscillator_strength), 4),
        }

        filename = "results_" + self.identifier + ".json"
        filepath = os.path.join(directory, filename)

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(dictionary, f, indent=2)
        print("Results saved at" + filepath)

    def calc_dipole_moment(self, MO_1=None, MO_2=None, unit="Coulomb*Angstrom"):
        r"""
        Calculate the dipole moment between two molecular orbitals: one from the HOMO, the other one from the LUMO.
        The function weigthens the overlap by the relative vector with respect to the center of mass.

        Parameters
        ----------
        MO_1 : numpy.ndarray, optional
            Molecular orbital coefficients for the first orbital. If None, defaults to HOMO.
        MO_2 : numpy.ndarray, optional
            Molecular orbital coefficients for the second orbital. If None, defaults to LUMO.
        unit : str, optional
            Unit of the dipole moment. Options are "Coulomb*Angstrom", "Debye", and "atomic_units".
            Default is "Coulomb*Angstrom".

        Returns
        -------
        float
            The dipole moment in the specified unit.

        Notes
        -----
        .. note::

            The dipole moment is calculated using :math:`d = \sqrt{d_x^2 + d_y^2 + d_z^2}`
            where dipole_x, dipole_y, and dipole_z are the components of the dipole moment
            along the x, y, and z axes, respectively.
            The conversion factors used are:

            - "Coulomb*Angstrom": No conversion needed.
            - "Debye": dipole * c.c / 1e-11
            - "atomic_units": dipole * 1e-10 / (c.physical_constants["Bohr radius"][0] * c.e)
        """

        if MO_1 is None and MO_2 is None:
            MO_1, MO_2 = self.HOMO, self.LUMO

        dipole_x = -c.e * self.relative_orbital_coordinates[:, 0] @ (MO_1 * MO_2)
        dipole_y = -c.e * self.relative_orbital_coordinates[:, 1] @ (MO_1 * MO_2)
        dipole_z = -c.e * self.relative_orbital_coordinates[:, 2] @ (MO_1 * MO_2)
        dipole = np.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)

        if unit == "Coulomb*Angstrom":
            return dipole
        if unit == "Debye":
            return dipole * c.c / 1e-11
        if unit == "atomic_units":
            return dipole * 1e-10 / (c.physical_constants["Bohr radius"][0] * c.e)
        return None

    def calc_oscillator_strength(self, dipole_moment):
        """
        Calculate the oscillator strength of a transition.

        Parameters
        ----------
        dipole_moment : float
            The dipole moment of the transition in Coulomb*Angstrom.

        Returns
        -------
        float
            The oscillator strength of the transition.
        """

        return (
            2
            / 3
            * c.m_e
            / (c.e**2 * c.hbar**2)
            * self.E_exc
            * c.e
            * (dipole_moment * 1e-10) ** 2
        )

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

    def _get_atoms(self):
        """Returns a list of atoms in the molecule as strings."""
        return list(self.xyz_data["Atom"])

    def _get_atom_coordinates(self):
        """ "Returns the coordinates of the atoms in the molecule."""
        atom_coordinates = np.zeros((self.num_atoms, 3))
        for atom_idx in range(self.num_atoms):
            xyz_entry = self.xyz_data.loc[atom_idx]
            atom_coordinates[atom_idx] = np.array(
                [xyz_entry["X"], xyz_entry["Y"], xyz_entry["Z"]]
            )
        return atom_coordinates

    def _get_orbitals(self):
        """Returns a list of orbitals in the molecule as instances of the Orbital
        class."""
        orbitals = []
        for atom_idx in range(self.num_atoms):
            coordinates = self.atom_coordinates[atom_idx]
            atom = self.atoms[atom_idx]
            orbital_types = []
            if atom == "H":
                orbital_types = ["s"]
            elif atom in ["C", "N", "O"]:
                orbital_types = ["s", "px", "py", "pz"]

            for orbital_type in orbital_types:
                atom_identifier = self.atom_identifiers[atom_idx]
                orbital = Orbital(coordinates, atom_identifier, orbital_type)
                orbitals.append(orbital)
        return orbitals

    def _get_num_orbitals(self):
        """Returns the number of orbitals in the molecule."""
        num_orbitals = 0
        for atom in self.atoms:
            num_orbitals += self.num_atom_orbitals[atom]
        return num_orbitals

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
        center_of_mass = np.zeros(3)
        molecule_mass = 0

        for atom_idx in range(self.num_atoms):
            coordinates = self.atom_coordinates[atom_idx]
            atom_mass = self.atom_masses[self.atoms[atom_idx]]
            center_of_mass += coordinates * atom_mass
            molecule_mass += atom_mass

        center_of_mass /= molecule_mass
        return molecule_mass, center_of_mass

    def _get_num_electrons(self):
        """ " Returns the number of valence electrons in the molecule."""
        num_electrons = 0
        for atom in self.atoms:
            num_electrons += self.num_atom_electrons[atom]
        return num_electrons

    def _calc_H(self):
        """Calculates the LCAO Hamiltonian matrix for the molecule."""
        H = np.zeros((self.num_orbitals, self.num_orbitals))

        # orbital overlap
        for i, j in combinations(range(self.num_orbitals), r=2):
            value = calc_orbital_overlap(
                self.orbitals[i], self.orbitals[j], connection="intrabase"
            )
            H[i, j] = value
            H[j, i] = value

        # orbital energy
        for i in range(self.num_orbitals):
            value = calc_orbital_energy(self.orbitals[i])
            H[i, i] = value

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
