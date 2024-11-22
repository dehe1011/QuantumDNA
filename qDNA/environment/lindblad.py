"""
This module provides classes and functions to construct Lindblad dissipators for master equations,
specifically for quantum DNA models.

Shortcuts
---------
diss: dissipator
loc: local
glob: global
deph: dephasing
therm: thermalizing
reorg: reorganization
eigv: eigenvalue
eigs: eigensystem
pop: population
"""

from itertools import product
import copy

import numpy as np
import qutip as q

from ..tools import DEFAULTS, UNITS, check_diss_kwargs
from ..utils import get_conversion
from ..hamiltonian import TB_Ham, add_groundstate

from .relax_ops import get_relax_ops
from .deph_ops import (
    get_loc_deph_ops,
    get_loc_deph_p_ops,
    get_glob_deph_ops,
    get_glob_deph_p_ops,
)
from .therm_ops import get_loc_therm_ops, get_glob_therm_ops
from .observables import get_eh_observable

__all__ = ["Lindblad_Diss"]

# -------------------------------------- Lindblad operators -------------------------------------


class Lindblad_Diss:
    """
    Provides the operators of the form :math:`\sqrt{\gamma} a` that are used in the Lindblad dissipator of the master equation.

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
        # Check inputs
        assert isinstance(
            tb_ham, TB_Ham
        ), "tb_ham must be an instance of the class TB_Ham"
        self.diss_kwargs = copy.copy(DEFAULTS["diss_kwargs_default"])
        self.diss_kwargs.update(diss_kwargs)
        check_diss_kwargs(**self.diss_kwargs)
        self.verbose = DEFAULTS["verbose"]

        if self.verbose:
            print("Successfully checked all inputs for the Lindblad_Diss instance.")

        # Tight-binding Hamiltonian
        self.tb_ham = tb_ham
        self.num_sites = self.tb_ham.tb_model.num_sites
        self._unit = self.tb_ham.unit

        # Dephasing
        self.loc_deph_rate = self.diss_kwargs.get("loc_deph_rate")
        self.glob_deph_rate = self.diss_kwargs.get("glob_deph_rate")

        # Dephasing operators
        self.deph_ops = self._get_deph_ops()

        # Relaxation
        self.uniform_relaxation = self.diss_kwargs.get("uniform_relaxation")
        if self.uniform_relaxation:
            relax_rate = self.diss_kwargs["relax_rate"]
            tb_sites = tb_ham.tb_sites_flattened
            self.relax_rates = dict(zip(tb_sites, [relax_rate] * len(tb_sites)))
        else:
            self.relax_rates = self.diss_kwargs.get("relax_rates")
            assert set(self.relax_rates.keys()) == set(
                tb_ham.tb_sites_flattened
            ), "relax_rates must have the same keys as the tight-binding sites"

        # Relaxation operators
        self.relax_ops = self._get_relax_ops()

        # Thermalization
        self.loc_therm = self.diss_kwargs.get("loc_therm")
        self.glob_therm = self.diss_kwargs.get("glob_therm")
        self.deph_rate = self.diss_kwargs.get("deph_rate")
        self.cutoff_freq = self.diss_kwargs.get("cutoff_freq")
        self.reorg_energy = self.diss_kwargs.get("reorg_energy")
        self.temperature = self.diss_kwargs.get("temperature")
        self.spectral_density = self.diss_kwargs.get("spectral_density")
        self.exponent = self.diss_kwargs.get("exponent")

        # Thermalization operators
        self.therm_ops = self._get_therm_ops()

        # Collapse operators
        self._c_ops = self.relax_ops + self.deph_ops + self.therm_ops
        self.num_c_ops = len(self.c_ops)

        # Observables: population, coherence, and ground state population operators
        if tb_ham.description == "2P":
            self.e_ops = self._get_e_ops()
            self.pop_ops, self.coh_ops, self.groundstate_pop_ops = self.e_ops

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

        # Update the number of collapse operators
        if new_c_ops != old_c_ops:
            self.num_c_ops = len(self.c_ops)

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, new_unit):
        assert isinstance(new_unit, str), "new_unit must be of type str"
        assert new_unit in UNITS, f"new_unit must be in {UNITS}"
        old_unit = self._unit
        self._unit = new_unit

        # Update the collapse operators
        if new_unit != old_unit:
            self.c_ops = [
                c_op * np.sqrt(get_conversion(old_unit, new_unit))
                for c_op in self.c_ops
            ]

    # -------------------------------------------------------------------

    def _get_relax_ops(self):
        """
        Generate relaxation operators based on relaxation rates and Hamiltonian.

        Returns
        -------
        list
            A list of relaxation operators, each scaled by the square root of the corresponding relaxation rate.
        """
        if not self.tb_ham.relaxation:
            return []

        tb_basis = self.tb_ham.tb_basis
        tb_basis_sites_dict = self.tb_ham.tb_basis_sites_dict

        return get_relax_ops(tb_basis, tb_basis_sites_dict, self.relax_rates)

    def _get_deph_ops(self):
        """
        Generate dephasing operators based on the system's Hamiltonian description and dephasing rates.
        This method computes local and global dephasing operators for both one-particle (1P) and
        two-particle (2P) descriptions of the system. The dephasing operators are determined by the
        eigensystem of the Hamiltonian and the specified local and global dephasing rates.

        Returns
        -------
        deph_p_ops : list
            A list of dephasing operators, scaled by the square root of the corresponding dephasing rate. The length and content of the list depend on the Hamiltonian
            description and the dephasing rates:
            - For local dephasing in a 2P description: 2 * num_sites operators.
            - For local dephasing in a 1P description: num_sites operators.
            - For global dephasing in a 2P description: num_sites^2 operators.
            - For global dephasing in a 1P description: num_sites operators.

        Raises
        ------
        AssertionError
            If the number of generated dephasing operators does not match the expected number based on
            the Hamiltonian description and the number of sites.
        """

        _, eigs = self.tb_ham.get_eigensystem()
        deph_ops = []

        # Local dephasing operators
        if self.loc_deph_rate:
            # Local dephasing operators for two particle description
            if self.tb_ham.description == "2P":
                loc_deph_ops = get_loc_deph_ops(
                    self.tb_ham.tb_basis, self.loc_deph_rate, self.tb_ham.relaxation
                )
                assert len(loc_deph_ops) == 2 * self.num_sites
                deph_ops = loc_deph_ops

            # Local dephasing operators for one particle description
            if self.tb_ham.description == "1P":
                loc_deph_p_ops = get_loc_deph_p_ops(
                    self.tb_ham.tb_basis, self.loc_deph_rate
                )
                assert len(loc_deph_p_ops) == self.num_sites
                deph_ops = loc_deph_p_ops

        # Global dephasing operators
        if self.glob_deph_rate:
            # Global dephasing operators for two particle description
            if self.tb_ham.description == "2P":
                glob_deph_ops = get_glob_deph_ops(
                    eigs, self.glob_deph_rate, self.tb_ham.relaxation
                )
                assert len(glob_deph_ops) == self.num_sites**2
                deph_ops = glob_deph_ops

            # Global dephasing operators for one particle description
            if self.tb_ham.description == "1P":
                glob_deph_p_ops = get_glob_deph_p_ops(eigs, self.glob_deph_rate)
                assert len(glob_deph_p_ops) == self.num_sites
                deph_ops = glob_deph_p_ops

        return deph_ops

    def _get_therm_ops(self):
        """
        Generate thermal operations based on the eigensystem of the Hamiltonian.
        This method calculates the thermal operations using either local or global
        thermalization processes, depending on the configuration of the instance.

        Returns
        -------
        therm_ops : list
            A list of thermal operations, scaled by the square root of the corresponding thermalization rate.

        Notes
        -----
        - The method first retrieves the eigensystem of the Hamiltonian and converts
          the eigenvalues to the desired units.
        - If `loc_therm` is set to True, local thermal operations are calculated.
        - If `glob_therm` is set to True, global thermal operations are calculated.
        - The thermal operations are determined by various parameters such as
          relaxation, dephasing rate, cutoff frequency, reorganization energy,
          temperature, spectral density, and exponent.
        """

        eigv, eigs = self.tb_ham.get_eigensystem()
        # Important: all parameters must be in the same units as the Hamiltonian. If you want to give them in "rad/ps" you should uncomment the following line.
        # eigv *= get_conversion(self.tb_ham.unit, "rad/ps")
        therm_ops = []

        # Local thermalizing operators
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

        # Global thermalizing operators
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
        """
        Generates the expectation value operators (e_ops) for the Lindblad master equation.
        This method constructs dictionaries of population and coherence operators for each particle
        in the tight-binding Hamiltonian (tb_ham). It also includes ground state population operators
        if relaxation is considered.

        Returns
        -------
        tuple of dict
            A tuple containing three dictionaries:
            - pop_dict: Population operators for each particle and site.
            - coh_dict: Coherence operators for each particle and pair of sites.
            - groundstate_pop_dict: Ground state population operators (only if relaxation is considered).

        Raises
        ------
        AssertionError
            If the Hamiltonian description is not "2P". This method is only defined for the "2P" description.
        """

        assert (
            self.tb_ham.description == "2P"
        ), "Only defined for 2P description. To treat the 1P case it is more efficient to return the whole density matrices and work with it."
        pop_dict, coh_dict, groundstate_pop_dict = {}, {}, {}

        # Population and coherence operators
        for particle in self.tb_ham.particles:
            for tb_site1, tb_site2 in product(self.tb_ham.tb_basis, repeat=2):
                # Create electron-hole observable
                observable = get_eh_observable(
                    self.tb_ham.tb_basis, particle, tb_site1, tb_site2
                )
                if self.tb_ham.relaxation:
                    observable = add_groundstate(observable)

                # Add observable as population operator
                if tb_site1 == tb_site2:
                    key = particle + "_" + tb_site1
                    pop_dict[key] = q.Qobj(observable)

                # Add observable as coherence operator
                else:
                    key = particle + "_" + tb_site1 + "_" + tb_site2
                    coh_dict[key] = q.Qobj(observable)

        # Ground state population operator
        if self.tb_ham.relaxation:
            groundstate_pop_dict["groundstate"] = q.fock_dm(self.tb_ham.matrix_dim, 0)

        return pop_dict, coh_dict, groundstate_pop_dict
