import numpy as np
import scipy.constants as c

from .slater_koster import calc_H_intra, calc_H_inter
from .component import Component

# ----------------------------------------------------------------------


class Monomer:
    """
    Monomer class representing a molecular system composed of multiple components.

    This class calculates properties of the monomer, including the LCAO Hamiltonian,
    molecular orbitals (HOMO and LUMO), center of mass, transition dipole moment, and
    orbital types.

    Parameters
    ----------
    filepaths : list of str
        List of filepaths for the components of the monomer.

    Attributes
    ----------
    lcao_param : dict
        s LCAO parametrization.
    filepaths : list of str
        List of file paths for the components of the monomer.
    components : list of Component
        List of Component objects representing the monomer's parts.
    num_components : int
        Number of components in the monomer.
    atoms : list
        List of atoms in the monomer.
    atoms_coordinates : list
        List of atomic coordinates in Angstroms.
    atoms_id : list
        List of atom identifiers.
    num_atoms : int
        Total number of atoms in the monomer.
    orbitals : list
        List of orbitals in the monomer.
    orbitals_coordinates : list
        List of orbital coordinates in Angstroms.
    num_orbitals : int
        Total number of orbitals in the monomer.
    num_electrons : int
        Total number of valence electrons in the monomer.
    H : ndarray
        LCAO Hamiltonian matrix.
    eigv : ndarray
        Eigenvalues of the Hamiltonian matrix.
    eigs : ndarray
        Eigenvectors of the Hamiltonian matrix.
    HOMO_idx : int
        Index of the HOMO orbital.
    E_HOMO : float
        Energy of the HOMO orbital.
    HOMO : ndarray
        Molecular orbital coefficients corresponding to HOMO.
    LUMO_idx : int
        Index of the LUMO orbital.
    E_LUMO : float
        Energy of the LUMO orbital.
    LUMO : ndarray
        Molecular orbital coefficients corresponding to LUMO.
    center_of_mass : ndarray
        Center of mass of the monomer in Angstroms.
    rel_orbitals_coordinates : ndarray
        Orbital coordinates relative to the center of mass in Angstroms.
    dipole_moment : ndarray
        Transition dipole moment between HOMO and LUMO.

    Notes
    -----
    All coordinates are given in Angstroms (Å), not nanometers (nm).

    """

    def __init__(self, filepaths, lcao_param):

        self.lcao_param = lcao_param

        # init components
        self.filepaths = filepaths
        self.components = [Component(filepath) for filepath in filepaths]
        self.num_components = len(filepaths)

        # monomer atoms
        self.atoms = [atom for comp in self.components for atom in comp.atoms]
        self.atoms_coordinates = [
            coord for comp in self.components for coord in comp.atoms_coordinates
        ]
        self.atoms_id = [aid for comp in self.components for aid in comp.atoms_id]
        self.num_atoms = sum(comp.num_atoms for comp in self.components)

        # monomer orbitals
        self.orbitals = [orbital for comp in self.components for orbital in comp.orbitals]
        self.orbitals_coordinates = [
            coord for comp in self.components for coord in comp.orbitals_coordinates
        ]
        self.num_orbitals = sum(comp.num_orbitals for comp in self.components)

        # valence electrons
        self.num_electrons = sum(comp.num_electrons for comp in self.components)

        # --------------------------------------------------------------

        # Calculate LCAO Hamiltonian
        self.H = self.calc_H()
        self.eigv, self.eigs = np.linalg.eigh(self.H)

        # Calculate HOMO and LUMO
        self.HOMO_idx = self.num_electrons // 2 - 1
        self.E_HOMO = self.eigv[self.HOMO_idx]
        self.HOMO = self.eigs[:, self.HOMO_idx]
        self.LUMO_idx = self.num_electrons // 2
        self.E_LUMO = self.eigv[self.LUMO_idx]
        self.LUMO = self.eigs[:, self.LUMO_idx]

        # Calculate transition dipole moment
        self.center_of_mass = self._calc_center_of_mass()
        self.rel_orbitals_coordinates = self._get_rel_orbitals_coordinates()
        self.dipole_moment = self._calc_dipole_moment()

    def __repr__(self):
        return f"Monomer({self.filepaths})"

    def _calc_center_of_mass(self):
        """Calculate the monomer's center of mass."""

        atom_masses = {"H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "X": 0}
        center_of_mass = np.zeros(3)
        molecule_mass = 0

        for atom_idx in range(self.num_atoms):
            atom = self.atoms[atom_idx]
            atom_coordinates = self.atoms_coordinates[atom_idx]
            atom_mass = atom_masses[atom]
            center_of_mass += atom_coordinates * atom_mass
            molecule_mass += atom_mass

        center_of_mass /= molecule_mass
        return center_of_mass

    def _get_rel_orbitals_coordinates(self):
        """Calculate orbital coordinates relative to the monomer's center of mass."""

        rel_orbitals_coordinates = []
        for orbital_coordinates in self.orbitals_coordinates:
            rel_orbitals_coordinates.append(orbital_coordinates - self.center_of_mass)
        return np.array(rel_orbitals_coordinates)

    # pylint: disable=inconsistent-return-statements
    def _calc_dipole_moment(self, unit="Coulomb*Angstrom"):
        """Calculate the monomer's transition dipole moment between HOMO and LUMO given
        the orbital coordinates relative to the monomer's center of mass."""

        MO_1, MO_2 = self.HOMO, self.LUMO

        x = self.rel_orbitals_coordinates[:, 0]
        y = self.rel_orbitals_coordinates[:, 1]
        z = self.rel_orbitals_coordinates[:, 2]

        dipole_x = -c.e * np.abs(MO_1) * x * np.abs(MO_2)
        dipole_y = -c.e * np.abs(MO_1) * y * np.abs(MO_2)
        dipole_z = -c.e * np.abs(MO_1) * z * np.abs(MO_2)

        dipole = np.array([np.sum(dipole_x), np.sum(dipole_y), np.sum(dipole_z)])

        if unit == "Coulomb*Angstrom":
            return dipole
        if unit == "Debye":
            return dipole * c.c / 1e-11
        if unit == "atomic_units":
            return dipole * 1e-10 / (c.physical_constants["Bohr radius"][0] * c.e)

    # ----------------------------------------

    def _build_block_matrix(self, D, U, L):
        """Build a block matrix from diagonal blocks (D), upper diagonal (U), lower diagonal (L)."""

        n = len(D)

        block_matrix = []
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(D[i])
                elif i == j - 1:
                    row.append(U[i])
                elif i == j + 1:
                    row.append(L[j])
                else:
                    shape = (D[i].shape[0], D[j].shape[1])
                    row.append(np.zeros(shape, dtype=D[0].dtype))
            block_matrix.append(row)

        return np.block(block_matrix)

    def calc_H(self):
        """Calculates the LCAO Hamiltonian matrix for the molecule."""

        # inter component Hamiltonian
        U, L = [], []
        for k in range(self.num_components - 1):
            comp1 = self.components[k]
            comp2 = self.components[k + 1]

            H_inter = calc_H_inter(self.lcao_param, comp1, comp2)
            U.append(H_inter)
            L.append(H_inter.conj().T)

        # intra component Hamiltonian
        D = []
        for k, comp in enumerate(self.components):
            H_intra = calc_H_intra(self.lcao_param, comp)
            D.append(H_intra)

        return self._build_block_matrix(D, U, L)

    def get_MO_type(self, MO):
        """Calculate the type of the molecular orbital: sigma, pi, n (non-bonding) together
        with the molecular orbital occupations."""

        # molecular orbital population
        MO_occupation = MO.conj() * MO

        # sigma orbital population
        s = ["C_s", "C_px", "C_py", "N_s", "O_s", "O_px", "O_py", "H_s"]
        s_mask = np.array([int(orbital in s) for orbital in self.orbitals])
        s_pop = sum(s_mask * MO_occupation)

        # pi orbital population
        pi = ["C_pz", "N_pz", "O_pz"]
        pi_mask = np.array([int(orbital in pi) for orbital in self.orbitals])
        pi_pop = sum(pi_mask * MO_occupation)

        # n (non-bonding) orbital population
        n = ["N_s", "N_px", "N_py"]
        n_mask = np.array([int(orbital in n) for orbital in self.orbitals])
        n_pop = sum(n_mask * MO_occupation)

        # molecular orbital type and occupation
        MO_type_occupation = [s_pop, pi_pop, n_pop]
        MO_type = ["sigma", "pi", "n"][np.argmax(MO_type_occupation)]

        return MO_type, MO_type_occupation

    # ------------------------------------------------------------------
