"""
This module provides the implementation of a tight-binding Hamiltonian for DNA sequences.
It includes the `TB_Ham` class, which represents the Hamiltonian and provides methods for its construction, manipulation, and analysis.
The module also includes several utility functions for setting matrix elements, constructing Hamiltonian matrices for single and two-particle systems, and adding interaction terms.

Shortcuts
---------
- tb: tight-binding
- nn: nearest-neighbor
- ham: hamiltonian
- param: parameter
- dim: dimension

"""

from itertools import chain
import copy

import numpy as np

from .. import DNA_Seq
from ..tools import check_ham_kwargs, DEFAULTS, UNITS, SOURCES
from ..utils import (
    calc_amplitudes,
    calc_average_pop,
    calc_frequencies,
    get_conversion,
    get_conversion_dict,
)
from ..model.tb_basis import get_eh_basis, get_particle_eh_states
from ..model.tb_model import TB_Model
from .tb_matrices import (
    tb_ham_1P,
    tb_ham_2P,
    add_groundstate,
    add_interaction,
    delete_groundstate,
)

from .tb_params import wrap_load_tb_params

__all__ = ["TB_Ham"]

# ------------------------------------------------------


class TB_Ham:
    """
    A class used to represent the tight-binding Hamiltonian for DNA sequences.

    Parameters
    ----------
    dna_seq : DNA_Seq
        The DNA sequence object containing the sequence and model information.
    ham_kwargs : dict
        Additional keyword arguments for the Hamiltonian construction.

    Attributes
    ----------
    dna_seq : DNA_Seq
        The DNA sequence object containing the sequence and model information.
    tb_model : TB_Model
        The tight-binding model associated with the DNA sequence.
    tb_basis : list
        The tight-binding basis states.
    tb_sites_flattened : list
        The flattened list of tight-binding sites.
    tb_basis_sites_dict : dict
        Dictionary mapping the TB basis to the TB sites.
    tb_sites : np.ndarray
        The array of TB sites.
    description : str
        Description of the Hamiltonian (e.g., '1P' or '2P').
    particles : list
        List of particles considered (e.g., ['electron', 'hole']).
    source : str
        Source of the tight-binding parameters.
    unit : str
        Unit of the Hamiltonian parameters.
    tb_params_electron : dict
        Tight-binding parameters for electrons.
    tb_params_hole : dict
        Tight-binding parameters for holes.
    relaxation : bool
        Flag indicating if relaxation is considered.
    interaction_param : float
        Interaction parameter for the Hamiltonian.
    eh_basis : list
        Electron-hole basis states.
    nn_cutoff : bool
        Nearest-neighbor cutoff for interactions.
    matrix : np.ndarray
        The Hamiltonian matrix.
    matrix_dim : int
        Dimension of the Hamiltonian matrix.
    backbone : bool
        Flag indicating if the model has a backbone.

    Methods
    -------
    get_param_dicts()
        Retrieves the tight-binding parameters for electrons and holes.
    get_eigensystem()
        Computes and returns the eigenvalues and eigenvectors of the Hamiltonian.
    get_matrix()
        Computes and returns the Hamiltonian matrix.
    get_fourier(init_state, end_state, quantities)
        Computes Fourier components of the Hamiltonian for given states.
    get_amplitudes(init_state, end_state)
        Computes and returns the amplitudes for given states.
    get_frequencies(init_state, end_state)
        Computes and returns the frequencies for given states.
    get_average_pop(init_state, end_state)
        Computes and returns the average population for given states.
    get_backbone_pop(init_state)
        Computes and returns the population of particles on the backbone sites.
    """

    def __init__(self, dna_seq, **ham_kwargs):
        # check inputs
        assert isinstance(dna_seq, DNA_Seq), "dna_seq must be an instance of DNA_Seq"
        self.ham_kwargs = copy.copy(DEFAULTS["ham_kwargs_default"])
        print(DEFAULTS)
        self.ham_kwargs.update(ham_kwargs)
        print(DEFAULTS)
        check_ham_kwargs(**self.ham_kwargs)
        self.verbose = DEFAULTS["verbose"]
        if self.verbose:
            print("Successfully checked all inputs for the TB_Ham instance.")

        # instance of DNA_Seq and TB_Model
        self.dna_seq = dna_seq
        self.tb_model = TB_Model(dna_seq.tb_model_name, dna_seq.tb_dims)

        # assigns each element of the DNA sequence to the corresponding tight-binding site
        self.tb_basis = self.tb_model.tb_basis
        self.tb_sites_flattened = list(chain(*self.dna_seq.dna_seq))
        self.tb_basis_sites_dict = dict(zip(self.tb_basis, self.tb_sites_flattened))
        self.tb_sites = np.array(self.tb_sites_flattened).reshape(
            self.tb_model.num_strands, self.tb_model.num_sites_per_strand
        )

        # Hamiltonian parameters
        self.description = self.ham_kwargs.get("description")
        self._particles = self.ham_kwargs.get("particles")
        self._source = self.ham_kwargs.get("source")
        self._unit = self.ham_kwargs.get("unit")

        # tight-binding parameters
        self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()

        self._relaxation = False
        if self.description == "2P":
            self._interaction_param = self.ham_kwargs.get("interaction_param")
            self._relaxation = self.ham_kwargs.get("relaxation")
            self.eh_basis = get_eh_basis(self.tb_model.tb_dims)
            self._nn_cutoff = self.ham_kwargs.get("nn_cutoff")
        self.matrix = self.get_matrix()
        self.matrix_dim = self.matrix.shape[0]
        self.backbone = True if self.tb_model.num_strands in (3, 4) else False

        if self.verbose:
            print("Successfully initialized the TB_Ham instance.")

    def __vars__(self) -> dict:
        """
        Returns the instance variables as a dictionary.
        """
        return vars(self)

    def __repr__(self) -> str:
        """
        Returns a string representation of the TB_Ham instance.
        """
        return f"TB_Ham({self.dna_seq}, {self.ham_kwargs})"

    def __eq__(self, other) -> bool:
        """
        Compares two TB_Ham instances for equality.
        """
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------------

    @property
    def particles(self):
        return self._particles

    @particles.setter
    def particles(self, new_particles):
        assert isinstance(new_particles, list), "new_particles must be of type list"
        assert all(
            [isinstance(new_particle, str) for new_particle in new_particles]
        ), "elements of new_particles must be of type str"
        self._particles = new_particles

    @property
    def interaction_param(self):
        return self._interaction_param

    @interaction_param.setter
    def interaction_param(self, new_interaction_param):
        assert isinstance(
            new_interaction_param, float
        ), "interaction_param must be of type float"
        old_interaction_param = self._interaction_param
        self._interaction_param = new_interaction_param

        # update the matrix
        if old_interaction_param != new_interaction_param:
            self.matrix = self.get_matrix()

    @property
    def relaxation(self):
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
    def nn_cutoff(self):
        return self._nn_cutoff

    @nn_cutoff.setter
    def nn_cutoff(self, new_nearest_neighbor_cutoff):
        old_nn_cutoff = self._nn_cutoff
        self._nn_cutoff = new_nearest_neighbor_cutoff

        # update the matrix
        if old_nn_cutoff != new_nearest_neighbor_cutoff:
            self.matrix = self.get_matrix()

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, new_unit):
        assert isinstance(new_unit, str), "new_unit must be of type str"
        assert new_unit in UNITS, f"new_unit must be in {UNITS}"
        old_unit = self._unit
        self._unit = new_unit

        # update the matrix and tight-binding parameters
        if new_unit != old_unit:
            self.matrix *= get_conversion(old_unit, new_unit)
            self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, new_source):
        assert isinstance(new_source, str), "new_source must be of type str"
        assert new_source in SOURCES, f"new_source must be in {SOURCES}"
        old_source = self._source
        self._source = new_source

        # update the matrix and tight-binding parameters
        if new_source != old_source:
            self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()
            self.matrix = self.get_matrix()

    # ---------------------------------------------------------------

    def get_param_dicts(self):
        """
        Retrieves the tight-binding parameters for electrons and holes.
        This method loads the tight-binding parameters for both electrons and holes
        from the specified source and model name. If the unit of the loaded parameters
        does not match the expected unit, it converts the parameters to the expected unit.

        Returns
        -------
        tuple
            A tuple containing two dictionaries:
              tb_params_electron : dict
                The tight-binding parameters for electrons.
              tb_params_hole : dict
                The tight-binding parameters for holes.
        """

        model_name = self.tb_model.tb_model_name

        # load the tight-binding parameters for electrons
        tb_params_electron, metadata = wrap_load_tb_params(
            self.source, "electron", model_name, load_metadata=True
        )
        # load the tight-binding parameters for holes
        tb_params_hole, metadata = wrap_load_tb_params(
            self.source, "hole", model_name, load_metadata=True
        )
        # convert the parameters to the expected unit
        if self.unit != metadata["unit"]:
            tb_params_electron = get_conversion_dict(
                tb_params_electron, metadata["unit"], self.unit
            )
            tb_params_hole = get_conversion_dict(
                tb_params_hole, metadata["unit"], self.unit
            )
        return tb_params_electron, tb_params_hole

    def get_eigensystem(self):
        """
        Compute the eigenvalues and eigenvectors of the matrix.
        This method computes the eigenvalues and eigenvectors of the matrix
        associated with the instance. If the description is "2P" and relaxation
        is enabled, the ground state is deleted from the matrix before computing
        the eigensystem.

        Returns
        -------
        tuple
            A tuple containing:
              eigenvalues : ndarray
                The eigenvalues of the matrix.
              eigenvectors : ndarray
                The eigenvectors of the matrix.
        """

        matrix = self.matrix

        # remove the ground state if relaxation is enabled
        if self.description == "2P":
            if self.relaxation:
                matrix = delete_groundstate(matrix)

        return np.linalg.eigh(matrix)

    def get_matrix(self):
        """
        Generate the tight-binding Hamiltonian matrix based on the system description.

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

        if self.description == "2P":
            # generate the Hamiltonian matrix for independent electron and hole
            matrix = tb_ham_2P(
                self.tb_model,
                self.tb_params_electron,
                self.tb_params_hole,
                self.tb_basis_sites_dict,
            )

            # add interaction terms
            if self.interaction_param:
                matrix = add_interaction(
                    matrix,
                    self.eh_basis,
                    self.interaction_param,
                    nn_cutoff=self.nn_cutoff,
                )

            # add relaxation terms
            if self._relaxation:
                matrix = add_groundstate(matrix)

        if self.description == "1P":
            # generate the Hamiltonian matrix for the electron
            if self.particles == ["electron"]:
                matrix = tb_ham_1P(
                    self.tb_model, self.tb_params_electron, self.tb_basis_sites_dict
                )

            # generate the Hamiltonian matrix for the hole
            if self.particles == ["hole"]:
                matrix = tb_ham_1P(
                    self.tb_model, self.tb_params_hole, self.tb_basis_sites_dict
                )

        return matrix

    def get_fourier(self, init_state, end_state, quantities):
        """
        Calculate the Fourier components of the transition between initial and end states.

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

        eigv, eigs = self.get_eigensystem()

        # check if the end state is in the tight-binding basis
        assert (
            end_state in self.tb_basis
        ), f"end_state must be in tb_basis {self.tb_basis}"

        # check if the initial state is in the electron-hole basis or tight-binding basis
        if self.description == "2P":
            assert (
                init_state in self.eh_basis
            ), f"init_state must be in tb_basis {self.eh_basis}"
            init_state_idx = self.eh_basis.index(init_state)

        if self.description == "1P":
            assert (
                init_state in self.tb_basis
            ), f"init_state must be in tb_basis {self.tb_basis}"
            init_state_idx = self.tb_basis.index(init_state)
            end_state_idx = self.tb_basis.index(end_state)

        amplitudes_dict, frequencies_dict, average_pop_dict = {}, {}, {}
        for particle in self.particles:
            if self.description == "2P":
                # get the end states for the particle
                end_states_idx = [
                    self.eh_basis.index(end_state)
                    for end_state in get_particle_eh_states(
                        particle, end_state, self.tb_basis
                    )
                ]

                # calculate amplitude
                if "amplitude" in quantities:
                    amplitudes_dict[particle] = list(
                        chain.from_iterable(
                            [
                                calc_amplitudes(eigs, init_state_idx, es)
                                for es in end_states_idx
                            ]
                        )
                    )

                # calculate frequency
                if "frequency" in quantities:
                    frequencies_dict[particle] = list(calc_frequencies(eigv)) * len(
                        end_states_idx
                    )

                # calculate average population
                if "average_pop" in quantities:
                    average_pop_dict[particle] = np.sum(
                        [
                            calc_average_pop(eigs, init_state_idx, es)
                            for es in end_states_idx
                        ]
                    )
            if self.description == "1P":
                # calculate amplitude
                if "amplitude" in quantities:
                    amplitudes_dict[particle] = calc_amplitudes(
                        eigs, init_state_idx, end_state_idx
                    )

                # calculate frequency
                if "frequency" in quantities:
                    frequencies_dict[particle] = calc_frequencies(eigv)

                # calculate average population
                if "average_pop" in quantities:
                    average_pop_dict[particle] = calc_average_pop(
                        eigs, init_state_idx, end_state_idx
                    )
        return amplitudes_dict, frequencies_dict, average_pop_dict

    def get_amplitudes(self, init_state, end_state):
        return self.get_fourier(init_state, end_state, ["amplitude"])[0]

    def get_frequencies(self, init_state, end_state):
        return self.get_fourier(init_state, end_state, ["frequency"])[1]

    def get_average_pop(self, init_state, end_state):
        return self.get_fourier(init_state, end_state, ["average_pop"])[2]

    def get_backbone_pop(self, init_state):
        """
        Calculate the population of particles on the backbone sites of a Fishbone model.

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

        # check if the model is a Fishbone model
        assert (
            self.tb_model.tb_model_name[0] == "F"
        ), "Backbone population can only be calculated for Fishbone models"

        # collect all backbone sites
        upper_backbone_sites = [
            f"(0, {site})" for site in range(self.tb_model.num_sites_per_strand)
        ]
        lower_backbone_sites = [
            f"({self.tb_model.num_strands-1}, {site})"
            for site in range(self.tb_model.num_sites_per_strand)
        ]
        backbone_sites = upper_backbone_sites + lower_backbone_sites

        # calculate the backbone population
        backbone_pop = dict(zip(self.particles, [0] * len(self.particles)))
        for particle in self.particles:
            for tb_site in backbone_sites:
                backbone_pop[particle] += self.get_average_pop(init_state, tb_site)[
                    particle
                ]
        return backbone_pop
