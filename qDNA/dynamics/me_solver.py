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

from ..environment import LindDiss, get_observable
from ..hamiltonian import add_groundstate
from ..io import DEFAULTS
from ..utils import check_me_solver_kwargs

from .reduced_dm import get_reduced_dm

__all__ = ["MeSolver"]

# ----------------------------------------------------------------------


class MeSolver(LindDiss):
    """
    MeSolver class for solving quantum master equations using QuTiP.

    This class extends the LindDiss class and provides functionality for solving
    quantum master equations with various configurations and parameters. It supports
    both one-particle (1P) and two-particle (2P) Hamiltonian descriptions, and includes
    methods for calculating populations, coherences, and ground state populations.

    Parameters
    ----------
    tb_sites : list
        List of tight-binding sites for the quantum system.

    Attributes
    ----------
    kwargs : dict
        Configuration parameters for the solver.
    times : numpy.ndarray
        Array of time points for the simulation.
    t_unit : str
        Unit of time for the simulation.
    ham_matrix : qutip.Qobj
        Hamiltonian matrix of the quantum system.
    init_states : list
        Initial states of the quantum system.
    init_matrix : qutip.Qobj
        Initial density matrix of the quantum system.
    result : list or None
        List of quantum states resulting from the master equation solver.
    groundstate_pop : dict or None
        Ground state population of the system.
    pop : dict or None
        Population of particles in the system.
    coh : dict or None
        Coherence of the system.
    options : dict
        Options for the QuTiP solver.
    solver_kwargs : dict
        Keyword arguments for the QuTiP solver.
    qutip_version : str
        Version of the QuTiP library being used.

    Notes
    -----
    .. note::

        - The class supports both QuTiP version 4 and version 5.
        - The Hamiltonian description can be either "1P" (one-particle) or "2P" (two-particle).
        - The solver uses QuTiP's `mesolve` function to solve the master equation.
    """

    def __init__(self, tb_sites, **kwargs):

        # Check kwargs
        self.kwargs = copy.copy(kwargs)
        self.kwargs.update(DEFAULTS["me_solver_kwargs_default"])
        self.kwargs.update(kwargs)
        check_me_solver_kwargs(**self.kwargs)

        # Initialize LindDiss
        super().__init__(tb_sites, **self.kwargs)

        # set the simulation time
        self._t_steps = int(self.kwargs.get("t_steps"))
        self._t_end = self.kwargs.get("t_end")
        self.times = np.linspace(0, self.t_end, self.t_steps)
        self.t_unit = self.kwargs.get("t_unit")
        assert self.t_steps / self.t_end > 1 / 2, (
            f"t_end {self.t_end} cannot be sufficiently resolved by t_steps {self.t_steps}. "
            "Please increase the number of steps or reduce the timespan. "
            "Alternative: change the unit of time from fs to ps (the mesolver does not know about the unit, but you do ;) )"
        )

        # TODO: ensure Hamiltonian and LindDiss match t_unit

        # set the Hamiltonian matrix (this is not unnecessary)
        self.ham_matrix = q.Qobj(self.get_matrix())

        # get initial state and initial density matrix
        self.init_states = self._get_init_states()
        self.init_matrix = self._get_init_matrix()

        # empty lists to store results
        self.result = None
        self.groundstate_pop = None
        self.pop = None
        self.coh = None
        for particle in self.particles:
            vars(self)["result_" + particle] = None

        # set options and qutip version for the solver
        self.options = {}
        self.solver_kwargs = {
            "H": self.ham_matrix,
            "rho0": self.init_matrix,
            "tlist": self.times,
            "c_ops": self.c_ops,
            "e_ops": {},
            "options": self.options,
        }
        self.qutip_version = q.__version__.split(".", maxsplit=1)[0]

    # ------------------------------------------------------------------

    def __vars__(self) -> dict:
        """Returns the instance variables as a dictionary."""
        return vars(self)

    def __repr__(self) -> str:
        """Returns a string representation of the MeSolver instance."""
        return f"MeSolver({self.tb_sites}, {self.kwargs})"

    def __eq__(self, other) -> bool:
        """Compares two MeSolver instances for equality."""
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

    # ------------------------------------------------------------------

    def reset(self):
        """Resets the solver's state by clearing results and initializing dictionaries
        for populations and coherences.

        Notes
        -----
        .. note::

            - Clears the ``result`` list (for the full the all reduced density matrices).
            - Initializes ``groundstate_pop``, ``pop``, and ``coh`` dictionaries.
        """

        self.result = None
        self.groundstate_pop = None
        self.pop = None
        self.coh = None
        for particle in self.particles:
            vars(self)["result_" + particle] = None
        self.solver_kwargs = {
            "H": self.ham_matrix,
            "rho0": self.init_matrix,
            "tlist": self.times,
            "c_ops": self.c_ops,
            "e_ops": {},
            "options": self.options,
        }

    def _get_init_states(self):

        init_e_states = self.kwargs.get("init_e_states")
        init_h_states = self.kwargs.get("init_h_states")
        init_ex_states = self.kwargs.get("init_ex_states")

        # set the initial state and iinitial density matrix
        if self.description == "2P":

            init_states = []
            for state_e in init_e_states:
                for state_h in init_h_states:
                    init_states.append((state_e, state_h))

        if self.description == "1P":
            if self.particles == ["electron"]:
                init_states = init_e_states
            if self.particles == ["hole"]:
                init_states = init_h_states
            if self.particles == ["exciton"]:
                init_states = init_ex_states
        return init_states

    def _get_init_matrix(self):
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

        init_matrix = 0

        if self.description == "2P":
            for init_state in self.init_states:
                state_matrix = get_observable(self.eh_basis, init_state, init_state)
                if self.relaxation:
                    state_matrix = add_groundstate(state_matrix)
                init_matrix += q.Qobj(state_matrix)

        if self.description == "1P":
            for init_state in self.init_states:
                init_state_idx = self.tb_basis.index(init_state)
                if self.relaxation:
                    state_matrix = q.fock_dm(self.matrix_dim, init_state_idx + 1)
                else:
                    state_matrix = q.fock_dm(self.matrix_dim, init_state_idx)
                init_matrix += q.Qobj(state_matrix)

        return init_matrix / len(self.init_states)

    # ------------------------------------------------------------------

    def _run_mesolve(self, **kwargs):
        """
        Run the mesolve function with the given arguments.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments to be passed to the mesolve function. The keys and values
            depend on the version of qutip being used.

        Returns
        -------
        list
            List of states resulting from the mesolve function.
        """

        if self.qutip_version == "5":
            kwargs["H"] = kwargs["H"].to(data_type="CSR")
            kwargs["rho0"] = kwargs["rho0"].to(data_type="CSR")
            kwargs["c_ops"] = [c_op.to(data_type="CSR") for c_op in kwargs["c_ops"]]
            kwargs["e_ops"] = {
                key: e_op.to(data_type="CSR") for key, e_op in kwargs["e_ops"].items()
            }
            kwargs["options"]["normalize_output"] = False
            kwargs["options"]["progress_bar"] = False

        if self.qutip_version == "4":
            kwargs["options"] = None

        if kwargs["e_ops"] == {}:
            return q.mesolve(**kwargs).states
        return q.mesolve(**kwargs)

    def get_result(self):
        """Calculate and return the result of the master equation solver. This method
        checks if the result has already been calculated. If not, it constructs the
        Hamiltonian matrix and solves the master equation using QuTiP's ``qutip.mesolve``
        function. The result is then stored and returned.

        Returns
        -------
        list
            A list of quantum states representing the solution of the master
            equation at different time points.
        """

        # check if the result is already calculated
        if self.result is not None:
            return self.result

        solver_kwargs = self.solver_kwargs.copy()
        solver_kwargs["e_ops"] = {}

        # store the result
        self.result = self._run_mesolve(
            **solver_kwargs
        )  # pylint: disable=attribute-defined-outside-init
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

        if vars(self)["result_" + particle] is not None:
            return vars(self)["result_" + particle]

        if self.result is None:
            self.get_result()

        reduced_dms = []
        for dm in self.result:
            reduced_dm = get_reduced_dm(dm, particle, self.tb_basis)
            reduced_dms.append(reduced_dm)
        vars(self)["result_" + particle] = reduced_dms

        return vars(self)["result_" + particle]

    # ------------------------------------------------------------------

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

            - If the population ``self.pop`` is already computed, it returns the cached result.
            - The method supports two types of Hamiltonian descriptions: "2P" and "1P".
            - The master equation is solved using ``qutip.mesolve`` with the Hamiltonian matrix, initial state, time points, collapse operators, and population operators.
        """

        if self.pop is not None:
            return self.pop
        self.pop = {}

        # solve the master equation with population observables
        solver_kwargs = self.solver_kwargs.copy()
        solver_kwargs["e_ops"] = self.pop_ops
        result = self._run_mesolve(**solver_kwargs)

        # store the population values
        for particle in self.particles:
            for tb_site in self.tb_basis:
                value = 0
                if self.qutip_version == "5":
                    value = result.e_data[particle + "_" + tb_site]
                if self.qutip_version == "4":
                    value = result.expect[particle + "_" + tb_site]
                self.pop[particle + "_" + tb_site] = value
        return self.pop

    def get_coh(self):
        """Calculate and return the coherence of the system.
        This method computes the coherence of the system based on the Hamiltonian
        description and the Lindblad dissipation operators. It supports two types
        of Hamiltonian descriptions: "2P" and "1P".
        For "2P" description, it uses the coherence operators from the Lindblad
        dissipation.
        For "1P" description, it constructs the coherence operators based on the
        tensor basis permutations.
        The method then solves the master equation using the QuTiP ``qutip.mesolve`` function
        and calculates the coherence for each particle in the system.

        Returns
        -------
        dict
            A dictionary where the keys are particle identifiers and the values are
            the computed coherence values.
        """

        if self.coh is not None:
            return self.coh
        self.coh = {}

        # solve the master equation with coherence observables
        solver_kwargs = self.solver_kwargs.copy()
        solver_kwargs["e_ops"] = self.coh_ops
        result = self._run_mesolve(**solver_kwargs)

        # store the coherence values
        for particle in self.particles:
            self.coh[particle] = 0
            for tb_site1, tb_site2 in permutations(self.tb_basis, 2):
                key = particle + "_" + tb_site1 + "_" + tb_site2
                if self.qutip_version == "5":
                    value = result.e_data[key]
                elif self.qutip_version == "4":
                    value = result.expect[key]
                else:
                    value = 0
                self.coh[particle] += np.abs(value)
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

        assert self.description == "2P", "only available for 2P description"
        assert self.relaxation, "only defined if relaxation is True"

        # check if the ground state population is already calculated
        if self.groundstate_pop is not None:
            return self.groundstate_pop
        self.groundstate_pop = {}

        # get observables for the ground state population
        solver_kwargs = self.solver_kwargs.copy()
        solver_kwargs["e_ops"] = self.groundstate_pop_ops
        result = self._run_mesolve(**solver_kwargs)

        # store the ground state population values
        value = 0
        if self.qutip_version == "5":
            value = result.e_data["groundstate"]
        if self.qutip_version == "4":
            value = result.expect["groundstate"]
        self.groundstate_pop["groundstate"] = value
        return self.groundstate_pop


# ----------------------------------------------------------------------
