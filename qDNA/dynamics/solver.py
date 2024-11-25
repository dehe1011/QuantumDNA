"""Module for solving master equations using the ME_Solver class.

Shortcuts
---------
- me: master equation
- diss: dissipator
- t: time
- init: initial
- pop: population
- coh: coherence
"""

from itertools import permutations
import copy

import numpy as np
import qutip as q

from .. import DNA_Seq
from ..environment import Lindblad_Diss, get_eh_observable
from ..hamiltonian import TB_Ham, add_groundstate
from ..tools import check_me_kwargs, DEFAULTS

from .reduced_dm import get_reduced_dm

__all__ = ["ME_Solver", "get_me_solver"]

# --------------------------------------------------------------------------------------


class ME_Solver:
    """A class used to solve master equations using the tight-binding Hamiltonian and
    Lindblad dissipator.

    This class provides methods to initialize the solver, set up the Hamiltonian and Lindblad dissipator,
    and compute various properties such as populations, coherences, and ground state populations over time.

    Parameters
    ----------
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    lindblad_diss : Lindblad_Diss
        The Lindblad dissipator.
    me_kwargs : dict, optional
        Additional keyword arguments for the master equation solver.

    Attributes
    ----------
    me_kwargs : dict
        Keyword arguments for the master equation solver.
    verbose : bool
        Flag for verbose output.
    t_steps : int
        Number of time steps.
    t_end : int
        End time.
    times : np.ndarray
        Array of time points.
    t_unit : str
        Time unit.
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    tb_model : TB_Model
        The tight-binding model.
    lindblad_diss : Lindblad_Diss
        The Lindblad dissipator.
    init_state : tuple or str
        Initial state depending on the Hamiltonian description.
    init_matrix : qutip.Qobj
        Initial density matrix.
    result : list
        List to store the results.
    groundstate_pop : dict
        Dictionary to store ground state population.
    pop : dict
        Dictionary to store population of states.
    coh : dict
        Dictionary to store coherence of states.

    Methods
    -------
    get_init_matrix()
        Generate the initial state matrix for the quantum system based on the Hamiltonian description.
    get_result()
        Calculate and return the result of the master equation solver.
    get_result_particle(particle)
        Retrieve the reduced density matrix for a specified particle.
    get_pop()
        Calculate and return the population of particles in the system.
    get_coh()
        Calculate and return the coherence of the system.
    get_groundstate_pop()
        Calculate and return the ground state population.
    reset()
        Resets the solver's state by clearing results and initializing dictionaries for populations and coherences.
    """

    def __init__(self, tb_ham, lindblad_diss, **me_kwargs):
        # Check inputs
        assert isinstance(
            tb_ham, TB_Ham
        ), "tb_ham must be an instance of the class TB_Ham"
        assert isinstance(
            lindblad_diss, Lindblad_Diss
        ), "lindblad_diss must be an instance of the class Lindblad_Diss"

        self.me_kwargs = copy.copy(DEFAULTS["me_kwargs_default"])
        self.me_kwargs.update(me_kwargs)
        check_me_kwargs(**self.me_kwargs)
        self.verbose = DEFAULTS["verbose"]
        if self.verbose:
            print("Successfully checked all inputs for the ME_Solver instance.")

        # set the simulation time
        self._t_steps = int(self.me_kwargs.get("t_steps"))
        self._t_end = int(self.me_kwargs.get("t_end"))
        self.times = np.linspace(0, self.t_end, self.t_steps)
        self.t_unit = self.me_kwargs.get("t_unit")
        assert self.t_steps / self.t_end > 1 / 2, (
            f"t_end {self.t_end} cannot be sufficiently resolved by t_steps {self.t_steps}. "
            "Please increase the number of steps or reduce the timespan. "
            "Alternative: change the unit of time from fs to ps (the mesolver does not know about the unit, but you do ;) )"
        )

        # initialize the Hamiltonian, Model and Lindblad dissipator
        self.tb_ham = tb_ham
        self.tb_ham.unit = "rad/" + self.t_unit
        self.tb_model = self.tb_ham.tb_model
        self.lindblad_diss = lindblad_diss
        self.lindblad_diss.unit = "rad/" + self.t_unit

        # set the Hamiltonian matrix
        self.ham_matrix = q.Qobj(self.tb_ham.matrix)

        # set the initial state and iinitial density matrix
        if self.tb_ham.description == "2P":
            self.init_state = (
                self.me_kwargs.get("init_e_state"),
                self.me_kwargs.get("init_h_state"),
            )
        if self.tb_ham.description == "1P":
            if self.tb_ham.particles == ["electron"]:
                self.init_state = self.me_kwargs.get("init_e_state")
            if self.tb_ham.particles == ["hole"]:
                self.init_state = self.me_kwargs.get("init_h_state")

        self.init_matrix = self.get_init_matrix()

        # set options for the solver
        self.options = None

        # empty lists to store results
        self.reset()

        if self.verbose:
            print("Successfully initialized the ME_Solver instance.")

    def __vars__(self) -> dict:
        """Returns the instance variables as a dictionary."""
        return vars(self)

    def __repr__(self) -> str:
        """Returns a string representation of the ME_Solver instance."""
        return f"ME_Solver({self.tb_ham}, {self.lindblad_diss}, {self.me_kwargs})"

    def __eq__(self, other) -> bool:
        """Compares two ME_Solver instances for equality."""
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------

    @property
    def t_end(self):  # pylint: disable=missing-function-docstring
        return self._t_end

    @t_end.setter
    def t_end(self, new_t_end):
        old_t_end = self._t_end
        self._t_end = new_t_end

        # update the time array and reset the results
        if new_t_end != old_t_end:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    @property
    def t_steps(self):  # pylint: disable=missing-function-docstring
        return self._t_steps

    @t_steps.setter
    def t_steps(self, new_t_steps):
        old_t_steps = self._t_steps
        self._t_steps = new_t_steps

        # update the time array and reset the results
        if new_t_steps != old_t_steps:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    # --------------------------------------------------------------------

    def reset(self):
        """Resets the solver's state by clearing results and initializing dictionaries
        for populations and coherences.

        Notes
        -----
        .. note::

            - Clears the ``result`` list (for the full the all reduced density matrices).
            - Initializes `groundstate_pop`, `pop`, and `coh` dictionaries.
        """

        self.result = []
        self.groundstate_pop = {}
        self.pop = {}
        self.coh = {}
        for particle in self.tb_ham.particles:
            vars(self)["result_" + particle] = []

    def get_init_matrix(self):
        """Generate the initial state matrix for the quantum system based on the
        Hamiltonian description. The method supports two types of descriptions for the
        tight-binding Hamiltonian (tb_ham): "2P" (two- particle) and "1P" (one-
        particle). Depending on the description and the initialization parameters, the
        initial state matrix is constructed either as a delocalized state over all
        exciton states or as a localized state on a single exciton state.

        Returns
        -------
        qutip.Qobj
            The initial state matrix of the quantum system as a Qobj instance from the QuTiP library.

        Raises
        ------
        ValueError
            If the Hamiltonian description is not recognized.
        """

        init_state = None

        # 2P description
        if self.tb_ham.description == "2P":
            # new: initial state delocalized over all exciton states. This correspopnds to the initial state of the whole system.
            if self.me_kwargs["deloc_init_state"]:
                tb_basis = self.tb_ham.tb_basis
                init_states = [
                    get_eh_observable(tb_basis, "exciton", state, state)
                    for state in tb_basis
                ]
                if self.tb_ham.relaxation:
                    init_states = [
                        add_groundstate(init_state) for init_state in init_states
                    ]
                init_state = 1 / len(tb_basis) * q.Qobj(np.sum(init_states, axis=0))

            # old: initial state localized on a single exciton state
            else:
                init_state_idx = self.tb_ham.eh_basis.index(self.init_state)
                if self.tb_ham.relaxation:
                    init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx + 1)
                else:
                    init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)

        # 1P description
        elif self.tb_ham.description == "1P":
            init_state_idx = self.tb_ham.tb_basis.index(self.init_state)
            init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)

        assert init_state is not None, "Initial state is not defined."
        return init_state

    def get_result(self):
        """Calculate and return the result of the master equation solver. This method
        checks if the result has already been calculated. If not, it constructs the
        Hamiltonian matrix and solves the master equation using QuTiP's `mesolve`
        function. The result is then stored and returned.

        Returns
        -------
        list
            A list of quantum states representing the solution of the master
            equation at different time points.
        """

        # check if the result is already calculated
        if not self.result:
            # observables
            e_ops = []

            # solve the master equation
            result = q.mesolve(
                self.ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                e_ops,
                progress_bar=None,
                options=self.options,
            ).states

            # store the result
            self.result = result  # pylint: disable=attribute-defined-outside-init
        return self.result

    def get_result_particle(self, particle):
        """Retrieve the reduced density matrix for a specified particle. This method
        checks if the result has already been calculated. If not, it calculates the
        result. Then, it checks if the reduced density matrix for the specified particle
        has been calculated. If not, it calculates the reduced density matrix for the
        specified particle and stores it.

        Parameters
        ----------
        particle : str
            The particle for which the reduced density matrix is to be retrieved.

        Returns
        -------
        list
            A list of reduced density matrices for the specified particle.
        """

        # check if the result is already calculated
        if not self.result:
            self.get_result()

        # check if the result is already calculated
        if not vars(self)["result_" + particle]:
            # calculate the reduced density matrix for the specified particle
            reduced_dms = [
                get_reduced_dm(dm, particle, self.tb_model.tb_basis)
                for dm in self.result
            ]

            # store the reduced density matrix
            vars(self)["result_" + particle] = reduced_dms
        return vars(self)["result_" + particle]

    def get_pop(self):
        """Calculate and return the population of particles in the system. This method
        computes the population of particles based on the Hamiltonian description and
        the Lindblad dissipation operators. It uses the QuTiP library to solve the
        master equation and obtain the expectation values of the population operators.

        Returns
        -------
        dict
            A dictionary where the keys are particle-site combinations and the
            values are the corresponding population expectation values.

        Notes
        -----
        .. note::

            - If the population (`self.pop`) is already computed, it returns the cached result.
            - The method supports two types of Hamiltonian descriptions: "2P" and "1P".
            - For "2P" description, it uses the population operators from `self.lindblad_diss.pop_ops`.
            - For "1P" description, it constructs the population operators based on the tight-binding basis (`self.tb_ham.tb_basis`).
            - The master equation is solved using `qutip.mesolve` with the Hamiltonian matrix, initial state, time points, collapse operators, and population operators.
        """

        # check if the population is already calculated
        if not self.pop:
            # observables for the population
            e_ops = None
            if self.tb_ham.description == "2P":
                e_ops = self.lindblad_diss.pop_ops
            elif self.tb_ham.description == "1P":
                keys = [
                    self.tb_ham.particles[0] + "_" + tb_site
                    for tb_site in self.tb_ham.tb_basis
                ]
                values = [
                    q.fock_dm(self.tb_ham.matrix_dim, i)
                    for i in range(self.tb_ham.matrix_dim)
                ]
                e_ops = dict(zip(keys, values))
            assert e_ops is not None, "Population operators are not defined."

            # solve the master equation with observables
            result = q.mesolve(
                self.ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                e_ops,
                options=self.options,
            )

            # store the population values
            for particle in self.tb_ham.particles:
                for tb_site in self.tb_ham.tb_basis:
                    self.pop[particle + "_" + tb_site] = result.expect[
                        particle + "_" + tb_site
                    ]
        return self.pop

    def get_coh(self):
        """
        Calculate and return the coherence of the system.
        This method computes the coherence of the system based on the Hamiltonian
        description and the Lindblad dissipation operators. It supports two types
        of Hamiltonian descriptions: "2P" and "1P".
        For "2P" description, it uses the coherence operators from the Lindblad
        dissipation.
        For "1P" description, it constructs the coherence operators based on the
        tensor basis permutations.
        The method then solves the master equation using the QuTiP `mesolve` function
        and calculates the coherence for each particle in the system.

        Returns
        -------
        dict
            A dictionary where the keys are particle identifiers and the values are
            the computed coherence values.
        """

        # check if the coherence is already calculated
        if not self.coh:
            # observables for the coherence
            e_ops = None
            if self.tb_ham.description == "2P":
                e_ops = self.lindblad_diss.coh_ops
            if self.tb_ham.description == "1P":
                keys = [
                    self.tb_ham.particles[0] + "_" + tb_site1 + "_" + tb_site2
                    for tb_site1, tb_site2 in permutations(self.tb_ham.tb_basis, 2)
                ]
                values = [
                    q.basis(self.tb_ham.matrix_dim, i)
                    * q.basis(self.tb_ham.matrix_dim, j).dag()
                    for i, j in permutations(self.tb_ham.matrix_dim, 2)
                ]
                e_ops = dict(zip(keys, values))
            assert e_ops is not None, "Coherence operators are not defined."

            # solve the master equation with observables
            result = q.mesolve(
                self.ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                e_ops,
                options=self.options,
            )

            # store the coherence values
            for particle in self.tb_ham.particles:
                self.coh[particle] = 0
                for tb_site1, tb_site2 in permutations(self.tb_ham.tb_basis, 2):
                    self.coh[particle] += abs(
                        result.expect[particle + "_" + tb_site1 + "_" + tb_site2]
                    )

        return self.coh

    def get_groundstate_pop(self):
        """Calculate and return the ground state population. This function computes the
        ground state population of a system described by a two- particle (2P)
        Hamiltonian with relaxation. If the ground state population has already been
        computed, it returns the cached result.

        Returns
        -------
        dict
            A dictionary containing the ground state population with the key
            "groundstate".

        Raises
        ------
        AssertionError
            If the Hamiltonian description is not "2P".
        AssertionError
            If relaxation is not enabled in the Hamiltonian.
        """

        assert self.tb_ham.description == "2P", "only available for 2P description"
        assert self.tb_ham.relaxation, "only defined if relaxation is True"

        # check if the ground state population is already calculated
        if not self.groundstate_pop:
            ham_matrix = q.Qobj(self.tb_ham.matrix)

            # observables for the ground state population
            e_ops = self.lindblad_diss.groundstate_pop_ops

            # solve the master equation with observables
            result = q.mesolve(
                ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                e_ops,
                options=self.options,
            )

            # store the ground state population values
            self.groundstate_pop["groundstate"] = result.expect["groundstate"]
        return self.groundstate_pop


# --------------------------------------------------------------------------------------


def get_me_solver(upper_strand, tb_model_name, **kwargs):
    """Creates an instance of ME_Solver.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    kwargs : dict
        Additional keyword arguments for creating the ME_Solver instance.

    Returns
    -------
    ME_Solver
        An instance of ME_Solver.
    """

    dna_seq = DNA_Seq(
        upper_strand, tb_model_name, lower_strand=kwargs.get("lower_strand")
    )
    tb_ham = TB_Ham(dna_seq, **kwargs)
    lindblad_diss = Lindblad_Diss(tb_ham, **kwargs)
    me_solver = ME_Solver(tb_ham, lindblad_diss, **kwargs)
    return me_solver
