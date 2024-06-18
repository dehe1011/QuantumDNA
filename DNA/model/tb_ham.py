import ast
import numpy as np
from itertools import chain
from DNA.model import get_eh_basis, basis_converter, wrap_load_tb_params, TBModelType
from utils import get_conversion_dict, get_conversion, get_config
from typing import List, Dict, Optional, Type, Tuple, Union

__all__ = ['PARTICLES', 'TB_Ham', 'TBHamType', 'set_matrix_element', 'particle_tb_Hamiltonian', 
           'eh_tb_Hamiltonian', 'add_basis_dimension', 'delete_basis_dimension', 'get_distance', 'add_interaction']

PARTICLES = ["electron", "hole", "exciton"]

# ------------------------------------------------------------------------------------------------------------

class TB_Ham:

    def __init__(self, tb_sites: List, tb_model: TBModelType, Ham_kwargs: Dict = {}):
        self.tb_model = tb_model
        self.tb_basis = self.tb_model.tb_basis
        self.tb_sites = list( chain(*tb_sites) )
        if len(self.tb_sites) != self.tb_model.num_sites:
            raise ValueError(f"The number of given tight-binding sites {len(self.tb_sites)} does not match "
            f"the dimension of the provided tight-binding model {self.tb_model.num_sites}")
        self.tb_basis_sites_dict = dict( zip( self.tb_basis, self.tb_sites) )
        self.tb_sites = np.array(self.tb_sites).reshape(self.tb_model.num_strands, self.tb_model.num_sites_per_strand)

        self.Ham_kwargs = get_config()["Ham_kwargs_default"]
        self.Ham_kwargs.update(Ham_kwargs)
        self.particle: str = self.Ham_kwargs.get("particle")
        if self.particle not in PARTICLES:
            raise ValueError(f"Unknown particle: {self.particle}. Predefined particles: {PARTICLES}")
        self._relaxation = False
        if self.particle == "exciton":
            self._interaction_param = self.Ham_kwargs.get("interaction_param")
            self._relaxation = self.Ham_kwargs.get("relaxation")
            self.eh_basis = get_eh_basis(self.tb_model.tb_dims)
            self._nearest_neighbor_cutoff = self.Ham_kwargs.get("nearest_neighbor_cutoff")

        self._unit = self.Ham_kwargs.get("unit") 
        self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()

        self._matrix = self.get_matrix()
        self.matrix_dim = self._matrix.shape[0]
        self.eigv, self.eigs = self.get_eigensystem()

    def __vars__(self) -> dict:
        return vars(self)

    @property
    def interaction_param(self):
        return self._interaction_param

    @interaction_param.setter
    def interaction_param(self, new_interaction_param):
        old_interaction_param = self._interaction_param
        self._interaction_param = new_interaction_param
        if old_interaction_param != new_interaction_param:
            self._matrix = self.get_matrix()
            self.eigv, self.eigs = self.get_eigensystem()

    @property
    def relaxation(self):
        return self._relaxation

    @relaxation.setter
    def relaxation(self, new_relaxation: bool):
        old_relaxation = self._relaxation
        self._relaxation = new_relaxation
        if new_relaxation != old_relaxation:
            if new_relaxation:
                self._matrix = add_basis_dimension(self.matrix)
            if not new_relaxation:
                self._matrix = delete_basis_dimension(self.matrix)
            self.matrix_dim = self._matrix.shape[0]
            self.eigv, self.eigs = self.get_eigensystem()

    @property
    def nearest_neighbor_cutoff(self):
        return self._nearest_neighbor_cutoff

    @nearest_neighbor_cutoff.setter
    def nearest_neighbor_cutoff(self, new_nearest_neighbor_cutoff):
        old_nearest_neighbor_cutoff = self._nearest_neighbor_cutoff 
        self._nearest_neighbor_cutoff  = new_nearest_neighbor_cutoff
        if old_nearest_neighbor_cutoff != new_nearest_neighbor_cutoff:
            self._matrix = self.get_matrix()
            self.eigv, self.eigs = self.get_eigensystem()

    @property
    def matrix(self):
        return self._matrix

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, new_unit):
        old_unit = self._unit
        self._unit = new_unit  
        if new_unit != old_unit:
            self._matrix *= get_conversion(old_unit, new_unit)
            self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()
            self.eigv *= get_conversion(old_unit, new_unit)

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, new_source):
        old_source = self._source
        self._source = new_source 
        if new_source != old_source:
            self.tb_params_electron, self.tb_params_hole = self.get_param_dicts()
            self._matrix = self.get_matrix()
            self.eigv, self.eigs = self.get_eigensystem()
            
    def get_param_dicts(self):
        self._source = self.Ham_kwargs.get("source")
        tb_params_electron, metadata = wrap_load_tb_params(
            self._source, "electron", self.tb_model.tb_model_name, load_metadata=True
        )
        tb_params_hole, metadata = wrap_load_tb_params(
            self._source, "hole", self.tb_model.tb_model_name, load_metadata=True
        )
        if self._unit != metadata["unit"]:
            tb_params_electron = get_conversion_dict(tb_params_electron, metadata["unit"], self._unit)
            tb_params_hole = get_conversion_dict(tb_params_hole, metadata["unit"], self._unit)
        return tb_params_electron, tb_params_hole

    def get_eigensystem(self) -> Tuple[np.ndarray, np.ndarray]:
            matrix = self._matrix
            if self.particle == 'exciton':
                if self.relaxation:
                    matrix = delete_basis_dimension(matrix)
            return np.linalg.eigh(matrix)

    def get_matrix(self) -> np.ndarray:
        """
        Computes and returns the Hamiltonian matrix.

        Returns:
            np.ndarray: Hamiltonian matrix.
        """
        if self.particle == "exciton":
            matrix = eh_tb_Hamiltonian(self.tb_model, self.tb_params_electron, self.tb_params_hole, self.tb_basis_sites_dict)
            if self.interaction_param:
                matrix = add_interaction(
                    matrix,
                    self.eh_basis,
                    self._interaction_param,
                    nearest_neighbor_cutoff=self._nearest_neighbor_cutoff,
                )
            if self._relaxation:
                matrix = add_basis_dimension(matrix)
            return matrix
        if self.particle == "electron":
            return particle_tb_Hamiltonian(self.tb_model, self.tb_params_electron, self.tb_basis_sites_dict)
        if self.particle == "hole":
            return particle_tb_Hamiltonian(self.tb_model, self.tb_params_hole, self.tb_basis_sites_dict)


TBHamType = Type[TB_Ham]

# ------------------------------------------------------------------------------------------------------------

def set_matrix_element(
    matrix: np.ndarray,
    tb_value: float,
    new_state: str,
    old_state: str,
    basis: List[str],
) -> np.ndarray:
    """
    Sets the matrix element for the Hamiltonian matrix ensuring hermiticity.

    Args:
        matrix (np.ndarray): The Hamiltonian matrix.
        tb_value (float): The tight-binding value to set.
        new_state (str): The new state in the basis.
        old_state (str): The old state in the basis.
        basis (List[str]): The list of basis states.

    Returns:
        np.ndarray: The updated Hamiltonian matrix.

    Example:
        >>> set_matrix_element(np.zeros((2, 2)), 1, '(0, 0)', '(1, 0)', ['(0, 0)', '(1, 0)'])
        array([[0., 1.],
               [1., 0.]])
    """
    old_state_num = basis_converter(old_state, basis)
    new_state_num = basis_converter(new_state, basis)
    matrix[new_state_num][old_state_num] += tb_value
    if old_state != new_state:
        matrix[old_state_num][new_state_num] += tb_value
    return matrix


def particle_tb_Hamiltonian(
    tb_model: TBModelType,
    tb_param_dict: Dict[str, float],
    tb_basis_sites_dict: Dict[str, str],
) -> np.ndarray:
    """
    Constructs the particle tight-binding Hamiltonian matrix.

    Args:
        tb_sites (List[str]): List of tight-binding sites.
        tb_param_dict (Dict[str, float]): Dictionary of tight-binding parameters.
        tb_config (List[Tuple[str, str, str]]): Tight-binding configuration.
        tb_site_basis (List[str]): Basis for the tight-binding sites.
        unit (str): Unit of the tight-binding parameter.
        num_strands (int): Number of strands in the system.

    Returns:
        np.ndarray: The tight-binding Hamiltonian matrix.
        
    """
    matrix = np.zeros((tb_model.num_sites, tb_model.num_sites))

    for tb_str, new_state, old_state in tb_model.tb_config:
        if tb_str == "E":
            tb_str = f"E_{tb_basis_sites_dict[old_state]}"
        else:
            tb_str = f"{tb_str}_{tb_basis_sites_dict[old_state]}{tb_basis_sites_dict[new_state]}"

        if tb_str not in tb_param_dict:
            raise ValueError(
                f"Tight-binding parameter '{tb_str}' not found in the parameter dictionary."
            )

        tb_val = tb_param_dict[tb_str]
        matrix = set_matrix_element(
            matrix, tb_val, new_state, old_state, tb_model.tb_basis
        )
    return matrix


def eh_tb_Hamiltonian(
    tb_model: TBModelType,
    tb_param_dict_electron: Dict[str, float],
    tb_param_dict_hole: Dict[str, float],
    tb_basis_sites_dict: Dict[str, str],
) -> np.ndarray:
    """
    Constructs the electron-hole tight-binding Hamiltonian matrix.

    Args:
        tb_sites (List[str]): List of tight-binding sites.
        tb_param_dict_electron (Dict[str, float]): Electron tight-binding parameters.
        tb_param_dict_hole (Dict[str, float]): Hole tight-binding parameters.
        tb_config (List[Tuple[str, str, str]]): Tight-binding configuration.
        tb_site_basis (List[str]): Basis for the tight-binding sites.
        unit (str): Unit of the tight-binding parameter.
        num_strands (int): Number of strands in the system.

    Returns:
        np.ndarray: The electron-hole tight-binding Hamiltonian matrix.

    Example:
        >>> tb_params_electron = {'E_G': 1.0, 'C_GC': 0.5}
        >>> tb_params_hole = {'E_G': 0.8, 'C_GC': 0.4}
        >>> tb_config = [('E', '(0, 0)', '(1, 0)'), ('C', '(1, 0)', '(0, 0)')]
        >>> eh_tb_Hamiltonian(['G','C'], tb_params_electron, tb_params_hole, tb_config, ['(0, 0)', '(1, 0)'], '100meV', 1)
        array([[0.5, 1.0, 0.4, 0.8],
               [1.0, 0.5, 0.8, 0.4],
               [0.4, 0.8, 0.5, 1.0],
               [0.8, 0.4, 1.0, 0.5]])
    """
    matrix_electron = particle_tb_Hamiltonian(tb_model, tb_param_dict_electron, tb_basis_sites_dict)
    matrix_hole = particle_tb_Hamiltonian(tb_model, tb_param_dict_hole, tb_basis_sites_dict)
    matrix = np.kron(np.eye(tb_model.num_sites), matrix_hole) + np.kron(matrix_electron, np.eye(tb_model.num_sites))
    return matrix

# ---------------------------------------------- relaxation ---------------------------------------


def add_basis_dimension(matrix: np.ndarray) -> np.ndarray:
    """
    Adds a dimension to the matrix to include the ground state.

    Args:
        matrix (np.ndarray): Input matrix.

    Returns:
        np.ndarray: Matrix with an additional dimension for the ground state.

    Example:
        >>> add_basis_dimension(np.array([[1, 2], [3, 4]]))
        array([[0., 0., 0.],
               [0., 1., 2.],
               [0., 3., 4.]])
    """
    N = matrix.shape[0]
    matrix = np.r_[np.zeros((1, N)), matrix]
    matrix = np.c_[np.zeros((N + 1, 1)), matrix]
    return matrix


def delete_basis_dimension(matrix: np.ndarray) -> np.ndarray:
    """
    Removes the ground state dimension from the matrix.

    Args:
        matrix (np.ndarray): Input matrix with ground state dimension.

    Returns:
        np.ndarray: Matrix without the ground state dimension.

    Example:
        >>> delete_basis_dimension(np.array([[0., 0., 0.], [0., 1., 2.], [0., 3., 4.]]))
        array([[1., 2.],
               [3., 4.]])
    """
    return matrix[1:, 1:]

# ------------------------------------------ interaction -----------------------------------------------

def get_distance(eh_basis: List[Tuple[str, str]]) -> np.ndarray:
    """
    Calculates the distance between electron and hole for each state in the basis.

    Args:
        eh_basis (List[Tuple[str, str]]): List of electron and hole positions as tuples of strings.

    Returns:
        np.ndarray: Array of distances between electron and hole for each state.

    Example:
        >>> get_distance([('(0, 0)', '(1, 1)'), ('(1, 0)', '(0, 0)')])
        array([1.41421356, 1.        ])
    """
    distance_list = []
    for position_electron, position_hole in eh_basis:
        position_electron = ast.literal_eval(position_electron)
        position_hole = ast.literal_eval(position_hole)
        distance = np.sqrt(sum( (idx_electron - idx_hole) ** 2 for idx_electron, idx_hole in zip(position_electron, position_hole)))
        distance_list.append(distance)
    return np.array(distance_list)

def add_interaction(
    matrix: np.ndarray,
    eh_basis: List[Tuple[str, str]],
    interaction_param: float,
    nearest_neighbor_cutoff: bool = False,
) -> np.ndarray:
    """
    Adds interaction terms to the Hamiltonian based on the distance between electron and hole.

    Args:
        Hamiltonian (np.ndarray): The initial Hamiltonian matrix.
        eh_basis (List[Tuple[str, str]]): List of electron and hole positions as tuples of strings.
        interaction_param (float): The interaction parameter.
        unit (str): Unit of the interaction parameter.
        nearest_neighbor_cutoff (bool, optional): If True, only nearest neighbor interactions are considered.

    Returns:
        np.ndarray: Hamiltonian matrix with interaction terms added.

    Note: 
        This works only for a hamiltonian without the additional basis element accounting for relaxation. Therefore the interaction
        should always be added before the relaxation. 

    Example:
        >>> Hamiltonian = np.array([[0, 1], [1, 0]])
        >>> eh_basis = [('(0, 0)', '(1, 1)'), ('(1, 0)', '(0, 0)')]
        >>> add_interaction(Hamiltonian, eh_basis, 1.0, '100meV')
        array([[0.        , 1.17639077],
               [1.17639077, 0.        ]])
    """

    distance_list = get_distance(eh_basis) 
    interaction_strength_list = (
        interaction_param / (1 + 3.4 * distance_list) 
    )
    if nearest_neighbor_cutoff:
        for i in np.where(distance_list > 1):
            interaction_strength_list[i] = 0

    for eh_basis_state, interaction_strength in enumerate(interaction_strength_list):
        eh_basis_state = basis_converter(eh_basis_state, eh_basis)
        matrix = set_matrix_element(
            matrix, interaction_strength, eh_basis_state, eh_basis_state, eh_basis
        )
    return matrix
