import os
import json
from itertools import combinations

import numpy as np
import scipy.constants as c

from .slater_koster import calc_orbital_energy, calc_orbital_overlap

# --------------------------------------------------


class Orbital:
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
        self.calc_oscillator_strength = self.calc_oscillator_strength(
            self.dipole_moment
        )  # no unit

    # ---------------------------------------------------------

    def __repr__(self):
        return f'Base("{self.identifier}")'

    def save_results(self, directory="results"):
        dict = {
            "E_HOMO": round(float(self.E_HOMO), 4),
            "E_LUMO": round(float(self.E_LUMO), 4),
            "HOMO": list(self.HOMO),
            "LUMO": list(self.LUMO),
            "dipole_moment": round(float(self.dipole_moment), 4),
            "oscillator_strength": round(float(self.calc_oscillator_strength), 4),
        }

        filename = "results_" + self.identifier + ".json"
        filepath = os.path.join(directory, filename)

        with open(filepath, "w") as f:
            json.dump(dict, f, indent=2)
        print("Results saved at" + filepath)

    def calc_dipole_moment(self, MO_1=None, MO_2=None, unit="Coulomb*Angstrom"):
        """Returns the dipole moment between two molecular orbitals: one from the HOMO, the other one from the LUMO.
        The function weigthens the overlap by the relative vector with respect to the center of mass.
        """

        if MO_1 == None and MO_2 == None:
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

    def calc_oscillator_strength(self, dipole_moment):
        """Returns the dimensionless oscillator strength.
        Note: The dipole moment should be provided in Coulomb*Angstrom"""
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
        """Returns properties of the MO.
        Note: deviation indicates lower levels than the MO (lower levels for HOMO, higher levels for the LUMO).
        """

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
        """Returns a list of orbitals in the molecule as instances of the Orbital class."""
        orbitals = []
        for atom_idx in range(self.num_atoms):
            coordinates = self.atom_coordinates[atom_idx]
            atom = self.atoms[atom_idx]
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
