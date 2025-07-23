import copy
import numpy as np
import scipy.constants as c

from ..utils import check_lcao_kwargs
from ..io import load_lcao_param, DEFAULTS
from .slater_koster import calc_H_intra, calc_H_inter
from .component import Component

# ----------------------------------------------------------------------


class Monomer:
    """
    Monomer class representing a molecular system composed of multiple components.
    This class calculates various properties of the molecular system, including
    the LCAO Hamiltonian, molecular orbitals (HOMO and LUMO), center of mass,
    transition dipole moment, and orbital types.

    Attributes
    ----------
    filepaths : list of str
        List of file paths for the components of the monomer.
    kwargs : dict
        Keyword arguments for LCAO parameter configuration.
    param_id : str
        Identifier for the LCAO parameter set.
    lcao_param : dict
        Loaded LCAO parameters.
    atom_masses : dict
        Atomic masses for supported elements.
    num_components : int
        Number of components in the monomer.
    components : list of Component
        List of Component objects representing the monomer's parts.
    xyz_id : list
        List of XYZ identifiers for the components.
    atoms : list
        List of atoms in the monomer.
    atoms_coordinates : list
        List of atomic coordinates.
    atoms_id : list
        List of atom identifiers.
    num_atoms : int
        Total number of atoms in the monomer.
    orbitals : list
        List of orbitals in the monomer.
    orbitals_coordinates : list
        List of orbital coordinates.
    num_orbitals : int
        Total number of orbitals in the monomer.
    num_electrons : int
        Total number of electrons in the monomer.
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
        Molecular orbital corresponding to HOMO.
    LUMO_idx : int
        Index of the LUMO orbital.
    E_LUMO : float
        Energy of the LUMO orbital.
    LUMO : ndarray
        Molecular orbital corresponding to LUMO.
    center_of_mass : ndarray
        Center of mass of the monomer.
    rel_orbitals_coordinates : ndarray
        Orbital coordinates relative to the center of mass.
    dipole_moment : ndarray
        Transition dipole moment.

    Methods
    -------
    calc_H()
        Calculates the LCAO Hamiltonian matrix for the molecule.
    build_block_matrix(D, U, L)
        Constructs a block matrix from diagonal, upper, and lower matrices.
    get_MO_type(MO)
        Determines the type of molecular orbital (sigma, pi, or non-bonding).

    Notes
    -----
    .. note::
        The coordinates are given in Angstroms (Å), not in nanometers (nm).

    """

    def __init__(self, filepaths, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(DEFAULTS["lcao_kwargs_default"])
        self.kwargs.update(kwargs)
        check_lcao_kwargs(**self.kwargs)

        self.filepaths = filepaths

        self.param_id = self.kwargs.get("param_id")
        self.lcao_param = load_lcao_param(self.param_id)

        self.atom_masses = {"H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "X": 0}

        self.num_components = len(filepaths)
        self.components = [Component(filepath, **self.kwargs) for filepath in filepaths]

        # Initialize attributes for the monomer
        self.xyz_id = []
        self.atoms, self.atoms_coordinates = [], []
        self.atoms_id = []
        self.num_atoms = 0
        self.orbitals, self.orbitals_coordinates = [], []
        self.num_orbitals = 0
        self.num_electrons = 0
        for comp in self.components:
            self.xyz_id.append(comp.xyz_id)
            self.atoms.extend(comp.atoms)
            self.atoms_coordinates.extend(comp.atoms_coordinates)
            self.atoms_id.extend(comp.atoms_id)
            self.num_atoms += comp.num_atoms
            self.orbitals.extend(comp.orbitals)
            self.orbitals_coordinates.extend(comp.orbitals_coordinates)
            self.num_orbitals += comp.num_orbitals
            self.num_electrons += comp.num_electrons

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

    def _get_rel_orbitals_coordinates(self):
        rel_orbitals_coordinates = []
        for orbital_coordinates in self.orbitals_coordinates:
            rel_orbitals_coordinates.append(orbital_coordinates - self.center_of_mass)
        return np.array(rel_orbitals_coordinates)

    def _calc_center_of_mass(self):
        center_of_mass = np.zeros(3)
        molecule_mass = 0

        for atom_idx in range(self.num_atoms):
            atom = self.atoms[atom_idx]
            atom_coordinates = self.atoms_coordinates[atom_idx]
            atom_mass = self.atom_masses[atom]
            center_of_mass += atom_coordinates * atom_mass
            molecule_mass += atom_mass

        center_of_mass /= molecule_mass
        return center_of_mass

    # pylint: disable=inconsistent-return-statements
    def _calc_dipole_moment(self, unit="Coulomb*Angstrom"):

        MO_1, MO_2 = self.HOMO, self.LUMO

        dipole_x = (
            -c.e * np.abs(MO_1) * self.rel_orbitals_coordinates[:, 0] * np.abs(MO_2)
        )
        dipole_y = (
            -c.e * np.abs(MO_1) * self.rel_orbitals_coordinates[:, 1] * np.abs(MO_2)
        )
        dipole_z = (
            -c.e * np.abs(MO_1) * self.rel_orbitals_coordinates[:, 2] * np.abs(MO_2)
        )
        dipole = np.array([np.sum(dipole_x), np.sum(dipole_y), np.sum(dipole_z)])

        if unit == "Coulomb*Angstrom":
            return dipole
        if unit == "Debye":
            return dipole * c.c / 1e-11
        if unit == "atomic_units":
            return dipole * 1e-10 / (c.physical_constants["Bohr radius"][0] * c.e)

    # ----------------------------------------

    def _build_block_matrix(self, D, U, L):
        n = len(D)

        full_matrix = []

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
            full_matrix.append(row)

        return np.block(full_matrix)

    def calc_H(self):
        """Calculates the LCAO Hamiltonian matrix for the molecule."""

        D, U, L = [], [], []

        for k in range(self.num_components - 1):
            comp1 = self.components[k]
            comp2 = self.components[k + 1]

            H_inter = calc_H_inter(self.lcao_param, comp1, comp2)
            U.append(H_inter)
            L.append(H_inter.conj().T)

        for k, comp in enumerate(self.components):
            H_intra = calc_H_intra(self.lcao_param, comp)
            D.append(H_intra)

        return self._build_block_matrix(D, U, L)

    def get_MO_type(self, MO):
        """Returns the nature of the molecular orbitals: sigma, pi, n (non-bonding)."""

        MO_occupation = MO.conj() * MO

        s = ["C_s", "C_px", "C_py", "N_s", "O_s", "O_px", "O_py", "H_s"]
        s_mask = np.array([int(orbital in s) for orbital in self.orbitals])
        pi = ["C_pz", "N_pz", "O_pz"]
        pi_mask = np.array([int(orbital in pi) for orbital in self.orbitals])
        n = ["N_s", "N_px", "N_py"]
        n_mask = np.array([int(orbital in n) for orbital in self.orbitals])

        s_pop = sum(s_mask * MO_occupation)
        pi_pop = sum(pi_mask * MO_occupation)
        n_pop = sum(n_mask * MO_occupation)

        MO_type = None
        if s_pop >= max(pi_pop, n_pop):
            MO_type = "sigma"

        if pi_pop >= max(s_pop, n_pop):
            MO_type = "pi"

        if n_pop >= max(s_pop, pi_pop):
            MO_type = "n"

        MO_type_occupation = [s_pop, pi_pop, n_pop]
        return MO_type, MO_type_occupation

    # ------------------------------------------------------------------
