import numpy as np
from ..model.tb_basis import get_eh_distance

__all__ = [
    "set_matrix_element",
    "tb_ham_1P",
    "tb_ham_2P",
    "add_groundstate",
    "delete_groundstate",
    "add_interaction",
]

# ------------------------------------------------------------


def set_matrix_element(
    matrix,
    tb_value,
    new_state,
    old_state,
    basis,
):
    """
    Sets the matrix element for the Hamiltonian matrix ensuring hermiticity.

    Parameters
    ----------
    matrix : np.ndarray
        The Hamiltonian matrix.
    tb_value : float
        The tight-binding value to set.
    new_state : str
        The new state in the basis.
    old_state : str
        The old state in the basis.
    basis : List[str]
        The list of basis states.

    Returns
    -------
    np.ndarray
        The updated Hamiltonian matrix.

    Example
    -------
    >>> set_matrix_element(np.zeros((2, 2)), 1, '(0, 0)', '(1, 0)', ['(0, 0)', '(1, 0)'])
    array([[0., 1.],
           [1., 0.]])
    """
    old_state_idx = basis.index(old_state)
    new_state_idx = basis.index(new_state)
    matrix[new_state_idx][old_state_idx] += tb_value

    # ensure hermiticity
    if old_state != new_state:
        matrix[old_state_idx][new_state_idx] += tb_value
    return matrix


def tb_ham_1P(
    tb_model,
    tb_param_dict,
    tb_basis_sites_dict,
):
    """
    Constructs the particle tight-binding Hamiltonian matrix.

    Parameters
    ----------
    tb_model : TBModelType
        The tight-binding model.
    tb_param_dict : Dict[str, float]
        Dictionary of tight-binding parameters.
    tb_basis_sites_dict : Dict[str, str]
        Dictionary mapping the TB basis to the TB sites.

    Returns
    -------
    np.ndarray
        The tight-binding Hamiltonian matrix.

    Example
    -------
    >>> tb_model = TB_Model('model_name', (2, 2))
    >>> tb_param_dict = {'E_G': 1.0, 'C_GC': 0.5}
    >>> tb_basis_sites_dict = {'(0, 0)': 'G', '(1, 0)': 'C'}
    >>> tb_ham_1P(tb_model, tb_param_dict, tb_basis_sites_dict)
    array([[1. , 0.5],
           [0.5, 1. ]])
    """
    matrix = np.zeros((tb_model.num_sites, tb_model.num_sites))

    for tb_str, new_state, old_state in tb_model.tb_config:
        if tb_str == "E":
            tb_str = f"E_{tb_basis_sites_dict[old_state]}"
        else:
            tb_str = f"{tb_str}_{tb_basis_sites_dict[old_state]}{tb_basis_sites_dict[new_state]}"
        if tb_str[0] == "h" and tb_str not in tb_param_dict:
            tb_str = f"{tb_str[0]}_{tb_basis_sites_dict[new_state]}{tb_basis_sites_dict[old_state]}"

        if tb_str not in tb_param_dict:
            raise ValueError(
                f"Tight-binding parameter '{tb_str}' not found in the parameter dictionary."
            )

        tb_val = tb_param_dict[tb_str]
        matrix = set_matrix_element(
            matrix, tb_val, new_state, old_state, tb_model.tb_basis
        )
    return matrix


def tb_ham_2P(
    tb_model,
    tb_param_dict_electron,
    tb_param_dict_hole,
    tb_basis_sites_dict,
):
    """
    Constructs the electron-hole tight-binding Hamiltonian matrix.

    Parameters
    ----------
    tb_model : TBModelType
        The tight-binding model.
    tb_param_dict_electron : Dict[str, float]
        Electron tight-binding parameters.
    tb_param_dict_hole : Dict[str, float]
        Hole tight-binding parameters.
    tb_basis_sites_dict : Dict[str, str]
        Dictionary mapping the TB basis to the TB sites.

    Returns
    -------
    np.ndarray
        The electron-hole tight-binding Hamiltonian matrix.

    Example
    -------
    >>> tb_model = TB_Model('model_name', (2, 2))
    >>> tb_param_dict_electron = {'E_G': 1.0, 'C_GC': 0.5}
    >>> tb_param_dict_hole = {'E_G': 0.8, 'C_GC': 0.4}
    >>> tb_basis_sites_dict = {'(0, 0)': 'G', '(1, 0)': 'C'}
    >>> tb_ham_2P(tb_model, tb_param_dict_electron, tb_param_dict_hole, tb_basis_sites_dict)
    array([[1.5, 0.5, 0.4, 0. ],
           [0.5, 1. , 0. , 0. ],
           [0.4, 0. , 1.3, 0.5],
           [0. , 0. , 0.5, 0.8]])
    """
    matrix_electron = tb_ham_1P(tb_model, tb_param_dict_electron, tb_basis_sites_dict)
    matrix_hole = tb_ham_1P(tb_model, tb_param_dict_hole, tb_basis_sites_dict)
    matrix = np.kron(np.eye(tb_model.num_sites), matrix_hole) + np.kron(
        matrix_electron, np.eye(tb_model.num_sites)
    )
    return matrix


def add_groundstate(matrix):
    """
    Adds a dimension to the matrix to include the ground state.

    Parameters
    ----------
    matrix : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Matrix with an additional dimension for the ground state.

    Example
    -------
    >>> add_groundstate(np.array([[1, 2], [3, 4]]))
    array([[0., 0., 0.],
           [0., 1., 2.],
           [0., 3., 4.]])
    """

    N = matrix.shape[0]
    matrix = np.r_[np.zeros((1, N)), matrix]
    matrix = np.c_[np.zeros((N + 1, 1)), matrix]
    return matrix


def delete_groundstate(matrix):
    """
    Removes the ground state dimension from the matrix.

    Parameters
    ----------
    matrix : np.ndarray
        Input matrix with ground state dimension.

    Returns
    -------
    np.ndarray
        Matrix without the ground state dimension.

    Example
    -------
    >>> delete_groundstate(np.array([[0., 0., 0.], [0., 1., 2.], [0., 3., 4.]]))
    array([[1., 2.],
           [3., 4.]])
    """

    return matrix[1:, 1:]


def add_interaction(
    matrix,
    eh_basis,
    interaction_param,
    nn_cutoff=False,
):
    """
    Adds interaction terms to the Hamiltonian based on the distance between electron and hole.

    Parameters
    ----------
    matrix : np.ndarray
        The initial Hamiltonian matrix.
    eh_basis : List[Tuple[str, str]]
        List of electron and hole positions as tuples of strings.
    interaction_param : float
        The interaction parameter.
    nn_cutoff : bool, optional
        If True, only nearest neighbor interactions are considered.

    Returns
    -------
    np.ndarray
        Hamiltonian matrix with interaction terms added.

    Note
    ----
    This works only for a Hamiltonian without the additional basis element accounting for relaxation.
    Therefore the interaction should always be added before the relaxation.

    Example
    -------
    >>> Hamiltonian = np.array([[0, 1], [1, 0]])
    >>> eh_basis = [('(0, 0)', '(1, 1)'), ('(1, 0)', '(0, 0)')]
    >>> add_interaction(Hamiltonian, eh_basis, 1.0, True)
    array([[0.        , 1.17639077],
           [1.17639077, 0.        ]])
    """

    distance_list = get_eh_distance(eh_basis)
    # as used by Bittner
    interaction_strength_list = interaction_param / (1 + 3.4 * distance_list)

    # nearest neighbor cutoff
    if nn_cutoff:
        for i in np.where(distance_list > 1):
            interaction_strength_list[i] = 0

    # add interaction terms on the diagonal for exciton states
    for eh_basis_state_idx, interaction_strength in enumerate(
        interaction_strength_list
    ):
        eh_basis_state = eh_basis[eh_basis_state_idx]
        matrix = set_matrix_element(
            matrix, interaction_strength, eh_basis_state, eh_basis_state, eh_basis
        )
    return matrix
