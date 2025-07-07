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
    A class representing a tight-binding Hamiltonian. This class
    inherits from :class:`.TBModel` and provides methods for constructing and manipulating the
    Hamiltonian matrix, as well as calculating various properties such as eigenvalues,
    eigenvectors, and Fourier components.

    Attributes
    ----------
    tb_sites : numpy.ndarray
        Array representing the tight-binding sites for the DNA sequence.
    tb_sites_flattened : numpy.ndarray
        Flattened array of tight-binding sites.
    tb_basis_sites_dict : dict
        Dictionary mapping tight-binding basis states to corresponding sites.
    sequence_id : str
        String representation of the DNA sequence.
    description : str
        System description ("1P" or "2P").
    particles : list of str
        List of particle types.
    source : str
        Source of the tight-binding parameters.
    unit : str
        Unit of the tight-binding parameters.
    tb_params : dict
        Tight-binding parameters for the system.
    matrix : numpy.ndarray or None
        Hamiltonian matrix for the system.
    matrix_dim : int or None
        Dimension of the Hamiltonian matrix.
    relaxation : bool
        Whether relaxation terms are included in the Hamiltonian matrix.
    coulomb_param : float or None
        Coulomb interaction parameter.
    exchange_param : float or None
        Exchange interaction parameter.
    nn_cutoff : float or None
        Nearest-neighbor cutoff distance for interactions.

    Notes
    -----
    .. note::

        - For "2P" description, the Hamiltonian matrix may include interaction and relaxation terms if specified.
        - For "1P" description, the Hamiltonian matrix is generated for either electrons or holes based on the `particles` attribute.

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
    def particles(self):  # pylint: disable=missing-function-docstring
        return self._particles

    @particles.setter
    def particles(self, new_particles):
        assert isinstance(new_particles, list), "new_particles must be of type list"
        assert all(
            isinstance(new_particle, str) for new_particle in new_particles
        ), "elements of new_particles must be of type str"
        self._particles = new_particles

    @property
    def coulomb_param(self):  # pylint: disable=missing-function-docstring
        return self._coulomb_param

    @coulomb_param.setter
    def coulomb_param(self, new_coulomb_param):
        assert isinstance(
            new_coulomb_param, float
        ), "coulomb_param must be of type float"
        old_coulomb_param = self._coulomb_param
        self._coulomb_param = new_coulomb_param

        # update the matrix
        if old_coulomb_param != new_coulomb_param:
            if self.matrix is not None:
                self.matrix = self.get_matrix()

    @property
    def exchange_param(self):  # pylint: disable=missing-function-docstring
        return self._exchange_param

    @exchange_param.setter
    def exchange_param(self, new_exchange_param):
        assert isinstance(
            new_exchange_param, float
        ), "exchange_param must be of type float"
        old_exchange_param = self._exchange_param
        self._coulomb_param = new_exchange_param

        # update the matrix
        if old_exchange_param != new_exchange_param:
            self.matrix = self.get_matrix()

    @property
    def relaxation(self):  # pylint: disable=missing-function-docstring
        return self._relaxation

    @relaxation.setter
    def relaxation(self, new_relaxation):
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
    def nn_cutoff(self):  # pylint: disable=missing-function-docstring
        return self._nn_cutoff

    @nn_cutoff.setter
    def nn_cutoff(self, new_nearest_neighbor_cutoff):
        old_nn_cutoff = self._nn_cutoff
        self._nn_cutoff = new_nearest_neighbor_cutoff

        # update the matrix
        if old_nn_cutoff != new_nearest_neighbor_cutoff:
            self.matrix = self.get_matrix()

    @property
    def unit(self):  # pylint: disable=missing-function-docstring
        return self._unit

    @unit.setter
    def unit(self, new_unit):
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
    def source(self):  # pylint: disable=missing-function-docstring
        return self._source

    @source.setter
    def source(self, new_source):
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
        """Retrieves the tight-binding parameters for the selected particle. This method
        loads the tight-binding parameters from the specified source and model name. If
        the unit of the loaded parameters does not match the expected unit, it converts
        the parameters to the expected unit.

        Parameters
        ----------
        particle : str
            The particle for which to retrieve the tight-binding parameters.

        Returns
        -------
        tuple
            tb_params : dict
                The tight-binding parameters.
        """

        tb_params, metadata = load_tb_params(
            self.source, self.tb_model_name, load_metadata=True
        )

        # convert the parameters to the expected unit
        if self.unit != metadata["unit"]:
            for key, value in tb_params.items():
                tb_params[key] = get_conversion_dict(value, metadata["unit"], self.unit)

        return tb_params

    def get_matrix(self):
        """Generate the tight-binding Hamiltonian matrix.

        Returns
        -------
        matrix : numpy.ndarray
            The Hamiltonian matrix for the system.

        Notes
        -----
        .. note::

            - For a "2P" description, the Hamiltonian matrix is generated using `tb_ham_2P` and may include interaction and relaxation terms if specified.
            - For a "1P" description, the Hamiltonian matrix is generated using `tb_ham_1P` for either electrons or holes based on the `particles` attribute.

        Raises
        ------
        ValueError
            If the `description` attribute is not "2P" or "1P".
        """

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
        """Compute the eigenvalues and eigenvectors of the matrix. This method computes
        the eigenvalues and eigenvectors of the matrix associated with the instance. If
        the description is "2P" and relaxation is enabled, the ground state is deleted
        from the matrix before computing the eigensystem.

        Returns
        -------
        tuple
            eigenvalues : ndarray
                The eigenvalues of the matrix.
            eigenvectors : ndarray
                The eigenvectors of the matrix.
        """

        # Compute the matrix if it has not been computed yet
        if self.matrix is None:
            self.get_matrix()

        matrix = self.matrix.copy()

        # Remove the ground state if relaxation is enabled
        if self.description == "2P" and self.relaxation:
            matrix = delete_groundstate(matrix)

        return np.linalg.eigh(matrix)

    def _get_fourier_1P(self, init_state, end_state, quantities):

        assert (
            init_state in self.tb_basis
        ), f"Initial state {init_state} must be in tb_basis."

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

        assert (
            init_state in self.eh_basis
        ), f"initial state {init_state} must be in tb_basis."

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
        """Calculate the Fourier components of the transition between initial and end
        states.

        Parameters
        ----------
        init_state : str
            The initial state from which the transition starts.
        end_state : str
            The end state to which the transition occurs.
        quantities : list of str
            List of quantities to calculate. Possible values are "amplitude", "frequency", and "average_pop".

        Returns
        -------
        amplitudes_dict : dict
            Dictionary containing the amplitudes for each particle.
        frequencies_dict : dict
            Dictionary containing the frequencies for each particle.
        average_pop_dict : dict
            Dictionary containing the average population for each particle.

        Raises
        ------
        AssertionError
            If `end_state` is not in `self.tb_basis`.
            If `init_state` is not in `self.eh_basis` or `self.tb_basis` depending on the description.
        """

        if quantities == "all":
            quantities = ["amplitude", "frequency", "average_pop"]

        # check if the end state is in the tight-binding basis
        assert end_state in self.tb_basis, f"end_state {end_state} must be in tb_basis."

        if self.description == "1P":
            return self._get_fourier_1P(init_state, end_state, quantities)

        if self.description == "2P":
            return self._get_fourier_2P(init_state, end_state, quantities)

    def get_amplitudes(
        self, init_state, end_state
    ):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["amplitude"])[0]

    def get_frequencies(
        self, init_state, end_state
    ):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["frequency"])[1]

    def get_average_pop(
        self, init_state, end_state
    ):  # pylint: disable=missing-function-docstring
        return self.get_fourier(init_state, end_state, ["average_pop"])[2]

    def get_backbone_average_pop(self, init_state):
        """Calculate the population of particles on the backbone sites of a Fishbone
        model.

        Parameters
        ----------
        init_state : str
            The initial state of the system.

        Returns
        -------
        dict
            A dictionary where keys are particle identifiers and values are the population
            of each particle on the backbone sites.

        Raises
        ------
        AssertionError
            If the model is not a Fishbone model.
        """

        assert (
            self.backbone
        ), "Backbone population can only be calculated for Fishbone models"

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
