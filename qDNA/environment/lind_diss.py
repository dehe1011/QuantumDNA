from itertools import product
import copy

import qutip as q

from ..io import DEFAULTS
from ..utils import check_lind_diss_kwargs
from ..hamiltonian import TBHam, add_groundstate

from .relax_ops import get_relax_ops
from .deph_ops import (
    get_loc_deph_ops,
    get_glob_deph_ops,
)
from .therm_ops import get_loc_therm_ops, get_glob_therm_ops
from .observables import get_eh_observable, get_tb_observable

__all__ = ["LindDiss"]

# TODO: ensure all rates are in the same units as the Hamiltonian.

# ----------------------------------------------------------------------


class LindDiss(TBHam):
    r"""Provides the operators of the form :math:`\sqrt{\gamma} a` that are used in the
    Lindblad dissipator of the master equation.

    Parameters
    ----------
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    kwargs : dict
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

    def __init__(self, tb_sites, **kwargs):

        # Check kwargs
        self.kwargs = copy.copy(kwargs)
        self.kwargs.update(DEFAULTS["lind_diss_kwargs_default"])
        self.kwargs.update(kwargs)
        check_lind_diss_kwargs(**self.kwargs)

        # Initialize TBHam
        super().__init__(tb_sites, **self.kwargs)

        # Dephasing parameters
        self.loc_deph_rate = self.kwargs.get("loc_deph_rate")
        self.glob_deph_rate = self.kwargs.get("glob_deph_rate")

        # Relaxation raparameterstes
        self.uniform_relaxation = self.kwargs.get("uniform_relaxation")
        if self.uniform_relaxation:
            relax_rate = self.kwargs["relax_rate"]
            tb_sites = self.tb_sites_flattened
            self._relax_rates = dict(zip(tb_sites, [relax_rate] * len(tb_sites)))
        else:
            self._relax_rates = self.kwargs.get("relax_rates")
            assert set(self.relax_rates.keys()) == set(
                self.tb_sites_flattened
            ), "relax_rates must have the same keys as the tight-binding sites"

        # Thermalization parameters
        self.loc_therm = self.kwargs.get("loc_therm")
        self.glob_therm = self.kwargs.get("glob_therm")
        self.deph_rate = self.kwargs.get("deph_rate")
        self.cutoff_freq = self.kwargs.get("cutoff_freq")
        self.reorg_energy = self.kwargs.get("reorg_energy")
        self.temperature = self.kwargs.get("temperature")
        self.spectral_density = self.kwargs.get("spectral_density")
        self.exponent = self.kwargs.get("exponent")

        # Lindblad operators
        self.relax_ops = self._get_relax_ops()
        self.deph_ops = self._get_deph_ops()
        self.therm_ops = self._get_therm_ops()
        self._c_ops = self.relax_ops + self.deph_ops + self.therm_ops
        self.num_c_ops = len(self.c_ops)

        # Observables
        self.e_ops = self._get_e_ops()
        self.pop_ops, self.coh_ops, self.groundstate_pop_ops = self.e_ops

    # ------------------------------------------------------------------

    def __vars__(self):
        """Returns the instance variables as a dictionary."""
        return vars(self)

    def __repr__(self):
        """Returns a string representation of the LindDiss instance."""
        return f"LindDiss({self.tb_sites}, {self.kwargs})"

    def __eq__(self, other):
        """Compares two LindDiss instances for equality."""
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------

    @property
    def c_ops(self):  # pylint: disable=missing-function-docstring
        """List of Lindblad operators (collapse operators) used in the Lindblad master equation."""
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
    def relax_rates(self):
        """Dictionary of relaxation rates for each DNA base."""
        return self._relax_rates

    @relax_rates.setter
    def relax_rates(self, new_relax_rates):
        assert isinstance(new_relax_rates, dict), "new_relax_rates must be of type dict"
        old_relax_rates = self._relax_rates
        self._relax_rates = new_relax_rates

        # Update the relaxation operators if the rates have changed
        if new_relax_rates != old_relax_rates:
            self.relax_ops = self._get_relax_ops()
            self._c_ops = self.relax_ops + self.deph_ops + self.therm_ops
            self.num_c_ops = len(self.c_ops)

    # ------------------------------------------------------------------

    def _get_relax_ops(self):
        """Generate relaxation operators based on relaxation rates and Hamiltonian.

        Returns
        -------
        list
            A list of relaxation operators, each scaled by the square root of the corresponding relaxation rate.
        """

        if not self.relaxation:
            return []

        return get_relax_ops(self.tb_basis, self.tb_basis_sites_dict, self.relax_rates)

    # pylint: disable=inconsistent-return-statements
    def _get_deph_ops(self):
        """Generate dephasing operators based on the system's Hamiltonian description
        and dephasing rates. This method computes local and global dephasing operators
        for both one-particle (1P) and two-particle (2P) descriptions of the system. The
        dephasing operators are determined by the eigensystem of the Hamiltonian and the
        specified local and global dephasing rates.

        Returns
        -------
        deph_p_ops : list
            A list of dephasing operators, scaled by the square root of the corresponding dephasing rate.
            The length and content of the list depend on the Hamiltonian description and the dephasing rates:
            - For local dephasing in a 2P description: 2 * num_sites operators.
            - For local dephasing in a 1P description: num_sites operators.
            - For global dephasing in a 2P description: num_sites^2 operators.
            - For global dephasing in a 1P description: num_sites operators.

        Raises
        ------
        AssertionError
            If the number of generated dephasing operators does not match the expected number
            based on the Hamiltonian description and the number of sites.
        """

        if not (self.loc_deph_rate or self.glob_deph_rate):
            return []

        _, eigs = self.get_eigensystem()

        # Local dephasing
        if self.loc_deph_rate:
            dephasing_dict = {p: self.loc_deph_rate for p in self.particles}
            deph_ops = get_loc_deph_ops(
                self.tb_basis, self.description, dephasing_dict, self.relaxation
            )
            assert len(deph_ops) == self.num_sites * len(self.particles)
            return deph_ops

        # Global dephasing (only for 2-particle description)
        if self.glob_deph_rate:
            deph_ops = get_glob_deph_ops(eigs, self.glob_deph_rate, self.relaxation)
            if self.description == "1P":
                assert len(deph_ops) == self.num_sites
            elif self.description == "2P":
                assert len(deph_ops) == self.num_sites**2
            return deph_ops

    # pylint: disable=inconsistent-return-statements
    def _get_therm_ops(self):
        """Generate thermal operations based on the eigensystem of the Hamiltonian. This
        method calculates the thermal operations using either local or global
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

        if not (self.loc_therm or self.glob_therm):
            return []

        eigv, eigs = self.get_eigensystem()
        # Important: all parameters must be in the same units as the Hamiltonian.
        # If you want to give them in "rad/ps" you should uncomment the following line.
        # eigv *= get_conversion(self.unit, "rad/ps")

        therm_kwargs = {
            "deph_rate": self.deph_rate,
            "cutoff_freq": self.cutoff_freq,
            "reorg_energy": self.reorg_energy,
            "temperature": self.temperature,
            "spectral_density": self.spectral_density,
            "exponent": self.exponent,
        }

        if self.loc_therm:
            return get_loc_therm_ops(eigv, eigs, self.relaxation, **therm_kwargs)

        if self.glob_therm:
            return get_glob_therm_ops(eigv, eigs, self.relaxation, **therm_kwargs)

    def _get_e_ops(self):
        """Generates the expectation value operators (e_ops) for the Lindblad master
        equation. This method constructs dictionaries of population and coherence
        operators for each particle in the tight-binding Hamiltonian (tb_ham). It also
        includes ground state population operators if relaxation is considered.

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

        pop_dict, coh_dict, groundstate_pop_dict = {}, {}, {}

        # Population and coherence operators
        for particle in self.particles:
            for tb_site1, tb_site2 in product(self.tb_basis, repeat=2):
                if self.description == "2P":
                    obs = get_eh_observable(self.tb_basis, particle, tb_site1, tb_site2)
                    if self.relaxation:
                        obs = add_groundstate(obs)

                if self.description == "1P":
                    obs = get_tb_observable(self.tb_basis, tb_site1, tb_site2)

                # Add observable as population operator
                if tb_site1 == tb_site2:
                    key = particle + "_" + tb_site1
                    pop_dict[key] = q.Qobj(obs)

                # Add observable as coherence operator
                else:
                    key = particle + "_" + tb_site1 + "_" + tb_site2
                    coh_dict[key] = q.Qobj(obs)

        # Ground state population operator
        if self.relaxation:
            groundstate_pop_dict["groundstate"] = q.fock_dm(self.num_sites**2 + 1, 0)

        return pop_dict, coh_dict, groundstate_pop_dict


# ----------------------------------------------------------------------
