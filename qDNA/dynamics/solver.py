"""
Module for solving master equations using the ME_Solver class.
"""

from itertools import permutations
import numpy as np
import qutip as q

from qDNA.tools import get_config, check_me_kwargs
from .reduced_dm import get_reduced_dm
from qDNA import DNA_Seq
from qDNA.model import TB_Ham
from qDNA.environment import Lindblad_Diss

# Shortcuts
# me: master equation
# diss: dissipator
# t: time
# init: initial
# pop: population
# coh: coherence

__all__ = ["ME_Solver", "get_me_solver"]


class ME_Solver:
    """
    A class used to solve master equations using the tight-binding Hamiltonian and Lindblad dissipator.

    Parameters
    ----------
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    lindblad_diss : LindbladDissType
        The Lindblad dissipator.
    me_kwargs : dict
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
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    tb_model : TBModelType
        The tight-binding model.
    lindblad_diss : LindbladDissType
        The Lindblad dissipator.
    init_state : tuple or str
        Initial state depending on the Hamiltonian description.
    init_matrix : q.Qobj
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
        Gets the initial density matrix.

    get_result()
        Solves the master equation and returns the result.

    get_result_particle(particle)
        Gets the result for a specific particle.

    get_pop()
        Gets the population of states.

    get_coh()
        Gets the coherence of states.

    get_groundstate_pop()
        Gets the ground state population.
    """

    def __init__(self, tb_ham, lindblad_diss, **me_kwargs):
        assert isinstance(
            tb_ham, TB_Ham
        ), "tb_ham must be an instance of the class TB_Ham"
        assert isinstance(
            lindblad_diss, Lindblad_Diss
        ), "lindblad_diss must be an instance of the class Lindblad_Diss"
        self.me_kwargs = get_config()["me_kwargs_default"]
        self.me_kwargs.update(me_kwargs)
        check_me_kwargs(**self.me_kwargs)
        self.verbose = get_config()["verbose"]
        if self.verbose:
            print("Successfully checked all inputs for the ME_Solver instance.")

        self._t_steps = int(self.me_kwargs.get("t_steps"))
        self._t_end = int(self.me_kwargs.get("t_end"))
        self.times = np.linspace(0, self.t_end, self.t_steps)
        self.t_unit = self.me_kwargs.get("t_unit")
        assert self.t_steps / self.t_end > 1 / 2, (
            f"t_end {self.t_end} cannot be sufficiently resolved by t_steps {self.t_steps}. "
            "Please increase the number of steps or reduce the timespan. "
            "Alternative: change the unit of time from fs to ps (the mesolver does not know about the unit, but you do ;) )"
        )

        self.tb_ham = tb_ham
        self.tb_ham.unit = "rad/" + self.t_unit
        self.tb_model = self.tb_ham.tb_model
        self.lindblad_diss = lindblad_diss
        self.lindblad_diss.unit = "rad/" + self.t_unit

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
        self.options = q.Options(method=self.me_kwargs.get("solver_method"))
        self.reset()

        if self.verbose:
            print("Successfully initialized the ME_Solver instance.")

    def __vars__(self) -> dict:
        """
        Returns the instance variables as a dictionary.
        """
        return vars(self)

    def __repr__(self) -> str:
        """
        Returns a string representation of the ME_Solver instance.
        """
        return f"ME_Solver({self.tb_ham}, {self.lindblad_diss}, {self.me_kwargs})"

    def __eq__(self, other) -> bool:
        """
        Compares two ME_Solver instances for equality.
        """
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------

    @property
    def t_end(self):
        return self._t_end

    @t_end.setter
    def t_end(self, new_t_end):
        old_t_end = self._t_end
        self._t_end = new_t_end
        if new_t_end != old_t_end:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    @property
    def t_steps(self):
        return self._t_steps

    @t_steps.setter
    def t_steps(self, new_t_steps):
        old_t_steps = self._t_steps
        self._t_steps = new_t_steps
        if new_t_steps != old_t_steps:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    # --------------------------------------------------------------------

    def reset(self):
        """
        Resets calculated results (for example after parameters were changed).
        """
        self.result = []
        if self.tb_ham.description == "2P":
            self.groundstate_pop = {}
            self.pop = {}
            self.coh = {}
            for particle in self.tb_ham.particles:
                vars(self)["result_" + particle] = []

    def get_init_matrix(self):
        """
        Gets the initial density matrix.

        Returns
        -------
        q.Qobj
            The initial density matrix.
        """
        if self.tb_ham.description == "2P":
            init_state_idx = self.tb_ham.eh_basis.index(self.init_state)
            if self.tb_ham.relaxation:
                init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx + 1)
            else:
                init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)
        if self.tb_ham.description == "1P":
            init_state_idx = self.tb_ham.tb_basis.index(self.init_state)
            init_state = q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)
        return init_state

    def get_result(self):
        """
        Solves the master equation and returns the result.

        Returns
        -------
        List[np.ndarray]
            The result of solving the master equation.
        """
        if not self.result:
            ham_matrix = q.Qobj(self.tb_ham.matrix)
            result = q.mesolve(
                ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                [],
                progress_bar=None,
                options=self.options,
            ).states
            self.result = result
        return self.result

    def get_result_particle(self, particle):
        """
        Gets the result for a specific particle.

        Parameters
        ----------
        particle : str
            The type of particle ('electron' or 'hole').

        Returns
        -------
        List[np.ndarray]
            The result for the specified particle.
        """
        if not self.result:
            self.get_result()
        if not vars(self)["result_" + particle]:
            vars(self)["result_" + particle] = [
                get_reduced_dm(dm, particle, self.tb_model.tb_basis)
                for dm in self.result
            ]
        return vars(self)["result_" + particle]

    def get_pop(self):
        """
        Gets the population of states.

        Returns
        -------
        Dict[str, float]
            The population of states.

        Raises
        ------
        AssertionError
            If the Hamiltonian description is not '2P'.
        """
        assert self.tb_ham.description == "2P", "only available for 2P description"
        if not self.pop:
            ham_matrix = q.Qobj(self.tb_ham.matrix)
            result = q.mesolve(
                ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                self.lindblad_diss.pop_ops,
                options=self.options,
            )
            for particle in self.tb_ham.particles:
                for tb_site in self.tb_ham.tb_basis:
                    self.pop[particle + "_" + tb_site] = result.expect[
                        particle + "_" + tb_site
                    ]
        return self.pop

    def get_coh(self):
        """
        Gets the coherence of states.

        Returns
        -------
        Dict[str, float]
            The coherence of states.

        Raises
        ------
        AssertionError
            If the Hamiltonian description is not '2P'.
        """
        assert self.tb_ham.description == "2P", "only available for 2P description"
        if not self.coh:
            ham_matrix = q.Qobj(self.tb_ham.matrix)
            result = q.mesolve(
                ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                self.lindblad_diss.coh_ops,
                options=self.options,
            )
            for particle in self.tb_ham.particles:
                self.coh[particle] = 0
                for tb_site1, tb_site2 in permutations(self.tb_ham.tb_basis, 2):
                    self.coh[particle] += abs(
                        result.expect[particle + "_" + tb_site1 + "_" + tb_site2]
                    )
        return self.coh

    def get_groundstate_pop(self):
        """
        Gets the ground state population.

        Returns
        -------
        Dict[str, float]
            The ground state population.

        Raises
        ------
        AssertionError
            If the Hamiltonian description is not '2P' or if relaxation is not enabled.
        """
        assert self.tb_ham.description == "2P", "only available for 2P description"
        if not self.groundstate_pop:
            assert (
                self.tb_ham.relaxation
            ), "function only defined if relaxation is True, otherwise the groundstate is not populated"
            ham_matrix = q.Qobj(self.tb_ham.matrix)
            result = q.mesolve(
                ham_matrix,
                self.init_matrix,
                self.times,
                self.lindblad_diss.c_ops,
                self.lindblad_diss.groundstate_pop_ops,
                options=self.options,
            )
            self.groundstate_pop["groundstate"] = result.expect["groundstate"]
        return self.groundstate_pop


# --------------------------------------------------------------------------------------


def get_me_solver(upper_strand, tb_model_name, **kwargs):
    """
    Creates an instance of ME_Solver.

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
    dna_seq = DNA_Seq(upper_strand, tb_model_name)
    tb_ham = TB_Ham(dna_seq, **kwargs)
    lindblad_diss = Lindblad_Diss(tb_ham, **kwargs)
    me_solver = ME_Solver(tb_ham, lindblad_diss, **kwargs)
    return me_solver
