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
    MeSolver class for solving master equations in quantum dynamics.
    This class extends LindDiss and provides methods for simulating quantum systems
    using master equations. It supports initialization of Hamiltonians, density matrices,
    and various observables, as well as solving the equations for populations, coherences,
    and ground state populations.

    Attributes
    ----------
    kwargs : dict
        Configuration parameters for the solver.
    times : ndarray
        Array of time points for the simulation.
    t_unit : str
        Unit of time used in the simulation.
    ham_matrix : qutip.Qobj
        Hamiltonian matrix of the system.
    init_states : list
        Initial states of the system.
    init_matrix : qutip.Qobj
        Initial density matrix of the system.
    result : list or None
        Results of the simulation.
    groundstate_pop : dict or None
        Ground state population values.
    pop : dict or None
        Population values for the system.
    coh : dict or None
        Coherence values for the system.
    options : dict
        Solver options.
    solver_kwargs : dict
        Arguments for the solver.
    qutip_version : str
        Version of QuTiP library used.
    t_end : float
    t_steps : int

    Methods
    -------
    reset()
        Resets the solver state and clears results.
    get_result()
        Solves the master equation and returns the state evolution.
    get_result_particle(particle)
        Returns the reduced density matrix for a specific particle.
    get_pop()
        Computes and returns population values.
    get_coh()
        Computes and returns coherence values.
    get_groundstate_pop()
        Computes and returns ground state population values.
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
    def t_end(self):
        """Returns the end time of the simulation."""
        return self._t_end

    @t_end.setter
    def t_end(self, new_t_end):  # pylint: disable=missing-function-docstring
        old_t_end = self._t_end
        self._t_end = new_t_end

        # update the time array and reset the results
        if new_t_end != old_t_end:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    @property
    def t_steps(self):
        """Returns the number of time steps in the simulation."""
        return self._t_steps

    @t_steps.setter
    def t_steps(self, new_t_steps):  # pylint: disable=missing-function-docstring
        old_t_steps = self._t_steps
        self._t_steps = new_t_steps

        # update the time array and reset the results
        if new_t_steps != old_t_steps:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    # ------------------------------------------------------------------

    def reset(self):

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
