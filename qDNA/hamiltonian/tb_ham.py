from itertools import chain
import copy

import numpy as np

from ..io import DEFAULTS, load_tb_params, OPTIONS
from ..utils import get_conversion, get_conversion_dict, check_tb_ham_kwargs
from ..model import get_eh_basis, get_particle_eh_states, TBModel
from .ham_analysis import calc_amplitudes, calc_average_pop, calc_frequencies
from .tb_matrices import (
    tb_ham_1P,
    tb_ham_2P,
    add_groundstate,
    add_interaction,
    delete_groundstate,
)

# ----------------------------------------------------------------------


class TBHam(TBModel):
    """
    TBHam class for modeling tight-binding Hamiltonians.
    This class extends the TBModel to represent tight-binding Hamiltonians for DNA structures.
    It supports both single-particle (1P) and two-particle (2P) descriptions, includes an
    option for the DNA relaxed state and can perform a Fourier analysis.

    Attributes
    ----------
    tb_sites : np.ndarray
        Array representing the tight-binding sites.
    tb_sites_flattened : np.ndarray
        Flattened array of tight-binding sites.
    tb_basis_sites_dict : dict
        Mapping of tight-binding basis to sites.
    sequence_id : str
        DNA sequence identifier.
    description : str
        Description of the Hamiltonian ("1P" or "2P").
    tb_params : dict
        Tight-binding parameters.
    matrix : np.ndarray or None
        Hamiltonian matrix.
    matrix_dim : int or None
        Dimension of the Hamiltonian matrix.
    relaxation : bool
    coulomb_param : float or None
    exchange_param : float or None
    nn_cutoff : float or None
    unit : str
    source : str
    particles : list

    Methods
    -------
    get_tb_params()
        Loads and converts tight-binding parameters.
    get_matrix()
        Computes the Hamiltonian matrix.
    get_eigensystem()
        Computes the eigenvalues and eigenvectors of the Hamiltonian matrix.
    get_fourier(init_state, end_state, quantities)
        Computes Fourier analysis for transitions between states.
    get_amplitudes(init_state, end_state)
        Computes transition amplitudes between states.
    get_frequencies(init_state, end_state)
        Computes transition frequencies between states.
    get_average_pop(init_state, end_state)
        Computes average population between states.
    get_backbone_average_pop(init_state)
        Computes average population for backbone sites.

    Examples
    --------
    >>> tb_sites = [["A", "T", "C"], ["G", "C", "A"]]
    >>> tb_ham = TBHam(tb_sites, description="1P", particles=["electron"], unit="eV")
    >>> tb_ham.get_matrix()

    """

    def __init__(self, tb_sites, **kwargs):

        # Check kwargs
        self.kwargs = copy.copy(kwargs)
        self.kwargs.update(DEFAULTS["tb_ham_kwargs_default"])
        self.kwargs.update(kwargs)
        check_tb_ham_kwargs(**self.kwargs)

        # Initialize TBModel
        num_sites_per_strand = len(tb_sites[0])
        super().__init__(num_sites_per_strand, **self.kwargs)

        # assigns each element of the DNA sequence to the corresponding tight-binding site
        self.tb_sites = np.array(tb_sites)
        self.tb_sites_flattened = self.tb_sites.flatten()
        self.tb_basis_sites_dict = dict(zip(self.tb_basis, self.tb_sites_flattened))

        if self.backbone:
            self.sequence_id = "".join(self.tb_sites[1, :])
        else:
            self.sequence_id = "".join(self.tb_sites[0, :])

        # Hamiltonian parameters
        self.description = self.kwargs.get("description")
        self._particles = self.kwargs.get("particles")
        self._source = self.kwargs.get("source")
        self._unit = self.kwargs.get("unit")

        # tight-binding parameters
        self.tb_params = self.get_tb_params()

        self._relaxation = False
        if self.description == "2P":
            self._coulomb_param = self.kwargs.get("coulomb_param")
            self._exchange_param = self.kwargs.get("exchange_param")
            self._relaxation = self.kwargs.get("relaxation")
            self.eh_basis = get_eh_basis(self.tb_dims)
            self._nn_cutoff = self.kwargs.get("nn_cutoff")

        # sa
        self.matrix = None
        self.matrix_dim = None

    # ------------------------------------------------------------------

    def __vars__(self):
        """Returns the instance variables as a dictionary."""
        return vars(self)

    def __repr__(self):
        """Returns a string representation of the TBHam instance."""
        return f"TBHam({self.tb_sites}, {self.kwargs})"

    def __eq__(self, other):
        """Compares two TBHam instances for equality."""
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------

    @property
    def particles(self):
        """Returns the particles in the Hamiltonian."""
        return self._particles

    @particles.setter
    def particles(self, new_particles):  # pylint: disable=missing-function-docstring
        assert isinstance(new_particles, list), "new_particles must be of type list"
        assert all(
            isinstance(new_particle, str) for new_particle in new_particles
        ), "elements of new_particles must be of type str"
        self._particles = new_particles

    @property
    def coulomb_param(self):
        """Returns the Coulomb interaction parameter."""
        return self._coulomb_param

    @coulomb_param.setter
    def coulomb_param(self, new_coulomb_param):  # pylint: disable=missing-function-docstring
        """Sets the Coulomb interaction parameter and updates the matrix if changed."""

        assert isinstance(new_coulomb_param, float), "coulomb_param must be of type float"
        old_coulomb_param = self._coulomb_param
        self._coulomb_param = new_coulomb_param

        # update the matrix
        if old_coulomb_param != new_coulomb_param:
            if self.matrix is not None:
                self.matrix = self.get_matrix()

    @property
    def exchange_param(self):
        """Returns the Exchange interaction parameter."""
        return self._exchange_param

    @exchange_param.setter
    def exchange_param(self, new_exchange_param):  # pylint: disable=missing-function-docstring
        assert isinstance(new_exchange_param, float), "exchange_param must be of type float"
        old_exchange_param = self._exchange_param
        self._coulomb_param = new_exchange_param

        # update the matrix
        if old_exchange_param != new_exchange_param:
            self.matrix = self.get_matrix()

    @property
    def relaxation(self):
        """Returns whether the DNA relaxed state is included."""
        return self._relaxation

    @relaxation.setter
    def relaxation(self, new_relaxation):  # pylint: disable=missing-function-docstring
        assert isinstance(new_relaxation, bool), "new_relaxation must be of type bool"
        old_relaxation = self._relaxation
        self._relaxation = new_relaxation

        # update the matrix
        if new_relaxation != old_relaxation:
            if new_relaxation:
                # add the ground state
                self.matrix = add_groundstate(self.matrix)
            if not new_relaxation:
                # remove the ground state
                self.matrix = delete_groundstate(self.matrix)
            self.matrix_dim = self.matrix.shape[0]

    @property
    def nn_cutoff(self):
        """Returns the nearest neighbor cutoff for interactions."""
        return self._nn_cutoff

    @nn_cutoff.setter
    def nn_cutoff(self, new_nearest_neighbor_cutoff):  # pylint: disable=missing-function-docstring
        old_nn_cutoff = self._nn_cutoff
        self._nn_cutoff = new_nearest_neighbor_cutoff

        # update the matrix
        if old_nn_cutoff != new_nearest_neighbor_cutoff:
            self.matrix = self.get_matrix()

    @property
    def unit(self):
        """Returns the unit of the Hamiltonian."""
        return self._unit

    @unit.setter
    def unit(self, new_unit):  # pylint: disable=missing-function-docstring
        units = OPTIONS["units"]
        assert isinstance(new_unit, str), "new_unit must be of type str"
        assert new_unit in units, f"new_unit must be in {units}"
        old_unit = self._unit
        self._unit = new_unit

        # update the matrix and tight-binding parameters
        if new_unit != old_unit:
            self.matrix *= get_conversion(old_unit, new_unit)
            self.tb_params = self.get_tb_params()

    @property
    def source(self):
        """Returns the source of the tight-binding parameters."""
        return self._source

    @source.setter
    def source(self, new_source):  # pylint: disable=missing-function-docstring
        sources = OPTIONS["sources"]
        assert isinstance(new_source, str), "new_source must be of type str"
        assert new_source in sources, f"new_source must be in {sources}"
        old_source = self._source
        self._source = new_source

        # update the matrix and tight-binding parameters
        if new_source != old_source:
            self.tb_params = self.get_tb_params()
            self.matrix = self.get_matrix()

    # ------------------------------------------------------------------

    def get_tb_params(self):
        tb_params, metadata = load_tb_params(self.source, self.tb_model_name, load_metadata=True)

        # convert the parameters to the expected unit
        if self.unit != metadata["unit"]:
            for key, value in tb_params.items():
                tb_params[key] = get_conversion_dict(value, metadata["unit"], self.unit)

        return tb_params

    def get_matrix(self):

        # Don't include this because the matrix cannot be overwritten
        # if self.matrix is not None:
        #     return self.matrix

        if self.description == "2P":
            # generate the Hamiltonian matrix for independent electron and hole
            self.matrix = tb_ham_2P(
                self.tb_dims,
                self.tb_config,
                self.tb_basis,
                self.tb_params,
                self.tb_basis_sites_dict,
            )

            # add interaction terms
            if self.coulomb_param:
                self.matrix = add_interaction(
                    self.matrix,
                    self.eh_basis,
                    self.coulomb_param,
                    "Coulomb",
                    nn_cutoff=self.nn_cutoff,
                )

            if self.exchange_param:
                self.matrix = add_interaction(
                    self.matrix,
                    self.eh_basis,
                    self.exchange_param,
                    "Exchange",
                    nn_cutoff=self.nn_cutoff,
                )

            # add relaxation terms
            if self._relaxation:
                self.matrix = add_groundstate(self.matrix)

        if self.description == "1P":
            particle = self.particles[0]
            self.matrix = tb_ham_1P(
                self.tb_dims,
                self.tb_config,
                self.tb_basis,
                self.tb_params[particle],
                self.tb_basis_sites_dict,
            )

        self.matrix_dim = self.matrix.shape[0]
        return self.matrix

    def get_eigensystem(self):

        # Compute the matrix if it has not been computed yet
        if self.matrix is None:
            self.get_matrix()

        matrix = self.matrix.copy()

        # Remove the ground state if relaxation is enabled
        if self.description == "2P" and self.relaxation:
            matrix = delete_groundstate(matrix)

        return np.linalg.eigh(matrix)

    def _get_fourier_1P(self, init_state, end_state, quantities):

        assert init_state in self.tb_basis, f"Initial state {init_state} must be in tb_basis."

        eigv, eigs = self.get_eigensystem()
        init_state_idx = self.tb_basis.index(init_state)
        end_state_idx = self.tb_basis.index(end_state)

        particle = self.particles[0]

        amplitudes_dict, frequencies_dict, average_pop_dict = {}, {}, {}
        if "amplitude" in quantities:
            val = calc_amplitudes(eigs, init_state_idx, end_state_idx)
            amplitudes_dict[particle] = val

        if "frequency" in quantities:
            val = calc_frequencies(eigv)
            frequencies_dict[particle] = val

        if "average_pop" in quantities:
            val = calc_average_pop(eigs, init_state_idx, end_state_idx)
            average_pop_dict[particle] = val
        return amplitudes_dict, frequencies_dict, average_pop_dict

    def _get_fourier_2P(self, init_state, end_state, quantities):

        assert init_state in self.eh_basis, f"initial state {init_state} must be in tb_basis."

        eigv, eigs = self.get_eigensystem()
        init_state_idx = self.eh_basis.index(init_state)

        amplitudes_dict, frequencies_dict, average_pop_dict = {}, {}, {}

        for particle in self.particles:

            # Determine the end state indices for each particle
            eh_states = get_particle_eh_states(particle, end_state, self.tb_basis)
            end_states_idx = []
            for eh_state in eh_states:
                end_states_idx.append(self.eh_basis.index(eh_state))

            if "amplitude" in quantities:

                amplitudes = []
                for end_state_idx in end_states_idx:
                    val = calc_amplitudes(eigs, init_state_idx, end_state_idx)
                    amplitudes.append(val)
                amplitudes_dict[particle] = list(chain.from_iterable(amplitudes))

            if "frequency" in quantities:
                val = calc_frequencies(eigv)
                frequencies_dict[particle] = list(val) * len(end_states_idx)

            if "average_pop" in quantities:
                average_pop = []
                for end_state_idx in end_states_idx:
                    val = calc_average_pop(eigs, init_state_idx, end_state_idx)
                    average_pop.append(val)
                average_pop_dict[particle] = np.sum(average_pop)
        return amplitudes_dict, frequencies_dict, average_pop_dict

    # pylint: disable=inconsistent-return-statements
    def get_fourier(self, init_state, end_state, quantities):

        if quantities == "all":
            quantities = ["amplitude", "frequency", "average_pop"]

        # check if the end state is in the tight-binding basis
        assert end_state in self.tb_basis, f"end_state {end_state} must be in tb_basis."

        if self.description == "1P":
            return self._get_fourier_1P(init_state, end_state, quantities)

        if self.description == "2P":
            return self._get_fourier_2P(init_state, end_state, quantities)

    def get_amplitudes(self, init_state, end_state):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["amplitude"])[0]

    def get_frequencies(self, init_state, end_state):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["frequency"])[1]

    def get_average_pop(self, init_state, end_state):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["average_pop"])[2]

    def get_backbone_average_pop(self, init_state):

        assert self.backbone, "Backbone population can only be calculated for Fishbone models"

        # collect all backbone sites
        upper_backbone_sites, lower_backbone_sites = [], []
        for site in range(self.num_sites_per_strand):
            upper_backbone_sites.append(f"(0, {site})")
            lower_backbone_sites.append(f"({self.num_channels-1}, {site})")
        backbone_sites = upper_backbone_sites + lower_backbone_sites

        # calculate the backbone population
        backbone_pop = dict(zip(self.particles, [0] * len(self.particles)))

        for tb_site in backbone_sites:
            val = self.get_average_pop(init_state, tb_site)
            for particle in self.particles:
                backbone_pop[particle] += val[particle]
        return backbone_pop


# ----------------------------------------------------------------------
