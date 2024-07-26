"""
This module provides classes and functions to construct Lindblad dissipators for master equations,
specifically for quantum DNA models.
"""

from itertools import product, permutations
import numpy as np
import qutip as q

from qDNA.tools import get_config, check_diss_kwargs
from qDNA.model import TB_Ham, add_groundstate, global_to_local
from qDNA.utils import get_conversion
from .lindblad_rates import rate_constant_redfield
from .observables import get_pop_particle, get_eh_observable

# Shortcuts
# diss: dissipator
# loc: local
# glob: global
# deph: dephasing
# therm: thermalizing
# reorg: reorganization
# eigv: eigenvalue
# eigs: eigensystem
# pop: population

__all__ = [
    "Lindblad_Diss",
    "get_relax_op",
    "get_glob_therm_op",
    "get_glob_therm_ops",
    "get_loc_therm_op",
    "get_loc_therm_ops",
    "get_loc_deph_ops",
    "get_glob_deph_ops",
    "get_loc_deph_p_ops",
    "get_glob_deph_p_ops",
]

# -------------------------------------- Lindblad operators -------------------------------------


class Lindblad_Diss:
    """
    Provides the operators of the form :math:`\sqrt{\\gamma} a` that are used in the Lindblad dissipator of the master equation.

    Parameters
    ----------
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    diss_kwargs : dict
        Additional keyword arguments for the dissipator.

    Attributes
    ----------
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    num_sites : int
        Number of sites in the tight-binding model.
    loc_deph_rate : float
        Local dephasing rate.
    glob_deph_rate : float
        Global dephasing rate.
    uniform_relaxation : bool
        Flag for uniform relaxation.
    relax_rates : dict
        Relaxation rates for each DNA base.
    loc_therm : bool
        Flag for local thermalization.
    glob_therm : bool
        Flag for global thermalization.
    deph_rate : float
        Dephasing rate.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature in Kelvin.
    spectral_density : str
        Type of spectral density ('debye' or 'ohmic').
    exponent : float
        Exponent for the Ohmic spectral density.
    relax_ops : list
        List of relaxation operators.
    deph_ops : list
        List of dephasing operators.
    therm_ops : list
        List of thermalizing operators.
    c_ops : list
        List of collapse operators.
    num_c_ops : int
        Number of collapse operators.
    e_ops : tuple
        Tuple containing population, coherence, and ground state population operators.
    pop_ops : dict
        Population operators.
    coh_ops : dict
        Coherence operators.
    groundstate_pop_ops : dict
        Ground state population operators.
    unit : str
        Unit of the operators.

    """

    def __init__(self, tb_ham, **diss_kwargs):
        assert isinstance(
            tb_ham, TB_Ham
        ), "tb_ham must be an instance of the class TB_Ham"
        self.diss_kwargs = get_config()["diss_kwargs_default"]
        self.diss_kwargs.update(diss_kwargs)
        check_diss_kwargs(**self.diss_kwargs)
        self.verbose = get_config()["verbose"]
        if self.verbose:
            print("Successfully checked all inputs for the Lindblad_Diss instance.")

        self.tb_ham = tb_ham
        self.num_sites = self.tb_ham.tb_model.num_sites

        self.loc_deph_rate = self.diss_kwargs.get("loc_deph_rate")
        self.glob_deph_rate = self.diss_kwargs.get("glob_deph_rate")

        self.uniform_relaxation = self.diss_kwargs.get("uniform_relaxation")
        if self.uniform_relaxation:
            DNA_BASES = get_config()["DNA_BASES"]
            DNA_SITES = DNA_BASES + ["B"]
            relax_rate = self.diss_kwargs["relax_rate"]
            self.relax_rates = dict(zip(DNA_SITES, [relax_rate] * len(DNA_SITES)))
        else:
            self.relax_rates = self.diss_kwargs.get("relax_rates")

        self.loc_therm = self.diss_kwargs.get("loc_therm")
        self.glob_therm = self.diss_kwargs.get("glob_therm")
        self.deph_rate = self.diss_kwargs.get("deph_rate")
        self.cutoff_freq = self.diss_kwargs.get("cutoff_freq")
        self.reorg_energy = self.diss_kwargs.get("reorg_energy")
        self.temperature = self.diss_kwargs.get("temperature")
        self.spectral_density = self.diss_kwargs.get("spectral_density")
        self.exponent = self.diss_kwargs.get("exponent")

        self.relax_ops = self._get_relax_ops()
        self.deph_ops = self._get_deph_ops()
        self.therm_ops = self._get_therm_ops()
        self._c_ops = self.relax_ops + self.deph_ops + self.therm_ops
        self.num_c_ops = len(self.c_ops)

        if tb_ham.description == "2P":
            self.e_ops = self._get_e_ops()
            self.pop_ops, self.coh_ops, self.groundstate_pop_ops = self.e_ops

        self._unit = self.tb_ham.unit

        if self.verbose:
            print("Successfully initialized the Lindblad_Diss instance.")

    def __vars__(self):
        """
        Returns the instance variables as a dictionary.
        """
        return vars(self)

    def __repr__(self):
        """
        Returns a string representation of the Lindblad_Diss instance.
        """
        return f"Lindblad_Diss({self.tb_ham}, {self.diss_kwargs})"

    def __eq__(self, other):
        """
        Compares two Lindblad_Diss instances for equality.
        """
        return self.__repr__() == other.__repr__()

    # --------------------------------------------------------------------

    @property
    def c_ops(self):
        return self._c_ops

    @c_ops.setter
    def c_ops(self, new_c_ops):
        assert isinstance(new_c_ops, list), "new_c_ops must be of type list"
        old_c_ops = self._c_ops
        self._c_ops = new_c_ops
        if new_c_ops != old_c_ops:
            self.num_c_ops = len(self.c_ops)

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, new_unit):
        assert isinstance(new_unit, str), "new_unit must be of type str"
        assert (
            new_unit in get_config()["UNITS"]
        ), f"new_unit must be in {get_config()['UNITS']}"
        old_unit = self._unit
        self._unit = new_unit
        if new_unit != old_unit:
            self.c_ops = [
                c_op * np.sqrt(get_conversion(old_unit, new_unit))
                for c_op in self.c_ops
            ]

    # -------------------------------------------------------------------

    def _get_relax_ops(self):
        relax_ops = []
        if self.tb_ham.relaxation and self.relax_rates:
            relax_ops = []
            for tb_site in self.tb_ham.tb_basis:
                relax_rate = self.relax_rates[self.tb_ham.tb_basis_sites_dict[tb_site]]
                if relax_rate != 0:
                    relax_ops.append(
                        np.sqrt(relax_rate)
                        * get_relax_op(self.tb_ham.tb_basis, tb_site)
                    )
        return relax_ops

    def _get_deph_ops(self):
        _, eigs = self.tb_ham.get_eigensystem()
        deph_p_ops = []
        if self.loc_deph_rate:
            if self.tb_ham.description == "2P":
                loc_deph_ops = get_loc_deph_ops(
                    self.tb_ham.tb_basis, self.loc_deph_rate, self.tb_ham.relaxation
                )
                assert len(loc_deph_ops) == 2 * self.num_sites
                deph_p_ops = loc_deph_ops
            if self.tb_ham.description == "1P":
                loc_deph_p_ops = get_loc_deph_p_ops(
                    self.tb_ham.tb_basis, self.loc_deph_rate
                )
                assert len(loc_deph_p_ops) == self.num_sites
                deph_p_ops = loc_deph_p_ops
        if self.glob_deph_rate:
            if self.tb_ham.description == "2P":
                glob_deph_ops = get_glob_deph_ops(
                    eigs, self.glob_deph_rate, self.tb_ham.relaxation
                )
                assert len(glob_deph_ops) == self.num_sites**2
                deph_p_ops = glob_deph_ops
            if self.tb_ham.description == "1P":
                glob_deph_p_ops = get_glob_deph_p_ops(eigs, self.glob_deph_rate)
                assert len(glob_deph_p_ops) == self.num_sites
                deph_p_ops = glob_deph_p_ops
        return deph_p_ops

    def _get_therm_ops(self):
        eigv, eigs = self.tb_ham.get_eigensystem()
        eigv *= get_conversion(self.tb_ham.unit, "rad/ps")
        therm_ops = []
        if self.loc_therm:
            loc_therm_ops = get_loc_therm_ops(
                eigv,
                eigs,
                self.tb_ham.relaxation,
                deph_rate=self.deph_rate,
                cutoff_freq=self.cutoff_freq,
                reorg_energy=self.reorg_energy,
                temperature=self.temperature,
                spectral_density=self.spectral_density,
                exponent=self.exponent,
            )
            therm_ops = loc_therm_ops
        if self.glob_therm:
            glob_therm_ops = get_glob_therm_ops(
                eigv,
                eigs,
                self.tb_ham.relaxation,
                deph_rate=self.deph_rate,
                cutoff_freq=self.cutoff_freq,
                reorg_energy=self.reorg_energy,
                temperature=self.temperature,
                spectral_density=self.spectral_density,
                exponent=self.exponent,
            )
            therm_ops = glob_therm_ops
        return therm_ops

    def _get_e_ops(self):
        assert (
            self.tb_ham.description == "2P"
        ), "Only defined for 2P description. To treat the 1P case it is more efficient to return the whole density matrices and work with them."
        pop_dict, coh_dict, groundstate_pop_dict = {}, {}, {}
        for particle in self.tb_ham.particles:
            for tb_site1, tb_site2 in product(self.tb_ham.tb_basis, repeat=2):
                observable = get_eh_observable(
                    self.tb_ham.tb_basis, particle, tb_site1, tb_site2
                )
                if self.tb_ham.relaxation:
                    observable = add_groundstate(observable)
                if tb_site1 == tb_site2:
                    pop_dict[particle + "_" + tb_site1] = q.Qobj(observable)
                else:
                    coh_dict[particle + "_" + tb_site1 + "_" + tb_site2] = q.Qobj(
                        observable
                    )
        if self.tb_ham.relaxation:
            groundstate_pop_dict["groundstate"] = q.fock_dm(self.tb_ham.matrix_dim, 0)
        return pop_dict, coh_dict, groundstate_pop_dict


# -------------------------------------- relaxation operator ------------------------------------


def get_relax_op(tb_basis, tb_site):
    """
    Annihilation operator of an exciton on a given tight-binding site. Relaxation of the DNA to its ground state.

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    tb_site : str
        Tight-binding site.

    Returns
    -------
    q.Qobj
        Relaxation operator.
    """
    tb_site_idx = tb_basis.index(tb_site)
    num_sites = len(tb_basis)
    relax_op = np.zeros((num_sites**2 + 1, num_sites**2 + 1))
    relax_op[0, 1 + tb_site_idx * (num_sites + 1)] = 1
    return q.Qobj(relax_op)


# ------------------------------------- thermalizing operators -----------------------------------


def get_glob_therm_op(eigs, eigenstate_i, eigenstate_j, relaxation, matrix_dim):
    """
    Global thermalizing operator.

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    eigenstate_i : int
        Index of the initial eigenstate.
    eigenstate_j : int
        Index of the final eigenstate.
    relaxation : bool
        Flag for relaxation.
    matrix_dim : int
        Dimension of the matrix.

    Returns
    -------
    q.Qobj
        Thermalizing operator.
    """
    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    op[eigenstate_i, eigenstate_j] = 1
    op = global_to_local(op, eigs)
    if relaxation:
        op = add_groundstate(op)
    return q.Qobj(op)


def get_glob_therm_ops(
    eigv,
    eigs,
    relaxation,
    deph_rate=7,
    cutoff_freq=20,
    reorg_energy=1,
    temperature=300,
    spectral_density="debye",
    exponent=1,
):
    """
    Global thermalizing operators.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalues.
    eigs : np.ndarray
        Eigensystem.
    relaxation : bool
        Flag for relaxation.
    deph_rate : float
        Dephasing rate.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature in Kelvin.
    spectral_density : str
        Type of spectral density.
    exponent : float
        Exponent for Ohmic spectral density.

    Returns
    -------
    list
        List of thermalizing operators.
    """
    matrix_dim = eigs.shape[0]
    c_ops = []
    for eigenstate_i, eigenstate_j in permutations(range(matrix_dim), 2):
        omega_i, omega_j = eigv[eigenstate_i], eigv[eigenstate_j]
        lind_rate = rate_constant_redfield(
            (omega_i - omega_j),
            deph_rate,
            cutoff_freq,
            reorg_energy,
            temperature,
            spectral_density,
            exponent,
        )
        lind_op = get_glob_therm_op(
            eigs, eigenstate_i, eigenstate_j, relaxation, matrix_dim
        )
        c_ops.append(np.sqrt(lind_rate) * lind_op)
    return c_ops


def get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim):
    """
    Local thermalizing operator.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalues.
    eigs : np.ndarray
        Eigensystem.
    unique : float
        Unique frequency gap.
    site_m : int
        Local site index.
    relaxation : bool
        Flag for relaxation.
    matrix_dim : int
        Dimension of the matrix.

    Returns
    -------
    q.Qobj
        Thermalizing operator.
    """
    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    for i, j in product(range(matrix_dim), repeat=2):
        omega_i, omega_j = eigv[i], eigv[j]
        state_i, state_j = eigs[:, i], eigs[:, j]
        if omega_i - omega_j == unique:
            op += (
                state_j[site_m].conjugate()
                * state_i[site_m]
                * np.outer(state_j, state_i)
            )
    if relaxation:
        op = add_groundstate(op)
    return q.Qobj(op)


def get_loc_therm_ops(
    eigv,
    eigs,
    relaxation,
    deph_rate=7,
    cutoff_freq=20,
    reorg_energy=1,
    temperature=300,
    spectral_density="debye",
    exponent=1,
):
    """
    Local thermalizing operators.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalues.
    eigs : np.ndarray
        Eigensystem.
    relaxation : bool
        Flag for relaxation.
    deph_rate : float
        Dephasing rate.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature in Kelvin.
    spectral_density : str
        Type of spectral density.
    exponent : float
        Exponent for Ohmic spectral density.

    Returns
    -------
    list
        List of thermalizing operators.
    """
    matrix_dim = len(eigv)
    c_ops = []
    gaps = eigv.reshape(matrix_dim, 1) - eigv
    unique = np.unique(gaps.flatten())
    for unique, site_m in product(unique, range(matrix_dim)):
        lind_rate = rate_constant_redfield(
            unique,
            deph_rate,
            cutoff_freq,
            reorg_energy,
            temperature,
            spectral_density,
            exponent,
        )
        lind_op = get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim)
        c_ops.append(np.sqrt(lind_rate) * lind_op)
    return c_ops


# ---------------------- local and global dephasing operators (for e-h and particle description) -------------------


def get_loc_deph_ops(tb_basis, dephasing_rate, relaxation):
    """
    Local dephasing operators.

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    dephasing_rate : float
        Dephasing rate.
    relaxation : bool
        Flag for relaxation.

    Returns
    -------
    list
        List of local dephasing operators.
    """
    c_ops = []
    for tb_site in tb_basis:
        op_electron = get_pop_particle(tb_basis, "electron", tb_site)
        op_hole = get_pop_particle(tb_basis, "hole", tb_site)
        if relaxation:
            op_electron = add_groundstate(op_electron)
            op_hole = add_groundstate(op_hole)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_electron))
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_hole))
    return c_ops


def get_glob_deph_ops(eigs, dephasing_rate, relaxation):
    """
    Global dephasing operators.

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    dephasing_rate : float
        Dephasing rate.
    relaxation : bool
        Flag for relaxation.

    Returns
    -------
    list
        List of global dephasing operators.
    """
    num_eigenstates = eigs.shape[0]
    c_ops = []
    for i in range(num_eigenstates):
        matrix = np.zeros((num_eigenstates, num_eigenstates))
        matrix[i, i] = 1
        op = global_to_local(matrix, eigs)
        if relaxation:
            op = add_groundstate(op)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op))
    return c_ops


def get_loc_deph_p_ops(tb_basis, dephasing_rate):
    """
    Local dephasing operators for particle description.

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    dephasing_rate : float
        Dephasing rate.

    Returns
    -------
    list
        List of local dephasing operators for particle description.
    """
    num_sites = len(tb_basis)
    return [np.sqrt(dephasing_rate) * q.fock_dm(num_sites, i) for i in range(num_sites)]


def get_glob_deph_p_ops(eigs, dephasing_rate):
    """
    Global dephasing operators for particle description.

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    dephasing_rate : float
        Dephasing rate.

    Returns
    -------
    list
        List of global dephasing operators for particle description.
    """
    num_eigenstates = eigs.shape[0]
    return [
        np.sqrt(dephasing_rate)
        * q.Qobj(global_to_local(q.fock_dm(num_eigenstates, i).full(), eigs))
        for i in range(num_eigenstates)
    ]
