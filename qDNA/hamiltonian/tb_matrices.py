from itertools import product

import numpy as np

from ..model import get_eh_distance

# ----------------------------------------------------------------------


def set_matrix_element(
    matrix,
    tb_value,
    new_state,
    old_state,
    basis,
):
    """Sets the matrix element for the Hamiltonian matrix ensuring hermiticity.

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

    Examples
    --------
    >>> set_matrix_element(np.zeros((2, 2)), 1, "(0, 0)", "(1, 0)", ["(0, 0)", "(1, 0)"])
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
    tb_dims,
    tb_config,
    tb_basis,
    tb_param_dict,
    tb_basis_sites_dict,
):
    """Constructs the particle tight-binding Hamiltonian matrix.

    Parameters
    ----------
    tb_dims : Tuple[int, int]
        Dimensions of the tight-binding model grid.
    tb_config : List[Tuple[str, str, str]]
        Configuration of tight-binding connections.
    tb_basis : List[str]
        List of basis states.
    tb_param_dict : Dict[str, float]
        Dictionary of tight-binding parameters.
    tb_basis_sites_dict : Dict[str, str]
        Dictionary mapping the TB basis to the TB sites.

    Returns
    -------
    np.ndarray
        The tight-binding Hamiltonian matrix.

    Examples
    --------
    >>> tb_dims = (2, 2)
    >>> tb_config = [("C", "(0, 0)", "(1, 0)"), ("E", "(0, 0)", "(0, 0)")]
    >>> tb_basis = ["(0, 0)", "(1, 0)", "(0, 1)", "(1, 1)"]
    >>> tb_param_dict = {"E_G": 1.0, "C_G_C": 0.5}
    >>> tb_basis_sites_dict = {"(0, 0)": "G", "(1, 0)": "C"}
    >>> tb_ham_1P(tb_dims, tb_config, tb_basis, tb_param_dict, tb_basis_sites_dict)
    array([[1. , 0.5, 0. , 0. ],
           [0.5, 0. , 0. , 0. ],
           [0. , 0. , 0. , 0. ],
           [0. , 0. , 0. , 0. ]])

    """

    n = tb_dims[0] * tb_dims[1]
    matrix = np.zeros((n, n))
    if tb_param_dict == {}:  # empty dictionary
        return matrix

    for tb_str, old_state, new_state in tb_config:
        old_state, new_state = str(old_state), str(new_state)
        if tb_str == "E":
            tb_str = f"E_{tb_basis_sites_dict[old_state]}"
        else:
            tb_str = f"{tb_str}_{tb_basis_sites_dict[old_state]}_{tb_basis_sites_dict[new_state]}"

        # for interstrand hopping the direction is not imporant
        if tb_str[0] in ["h", "r"] and tb_str not in tb_param_dict:
            tb_str = f"{tb_str.split('_')[0]}_{tb_str.split('_')[2]}_{tb_str.split('_')[1]}"

        if tb_str not in tb_param_dict:
            raise ValueError(
                f"Tight-binding parameter '{tb_str}' not found in the parameter dictionary."
            )

        tb_val = tb_param_dict[tb_str]
        matrix = set_matrix_element(matrix, tb_val, new_state, old_state, tb_basis)
    return matrix


def tb_ham_2P(
    tb_dims,
    tb_config,
    tb_basis,
    tb_params,
    tb_basis_sites_dict,
):
    """Constructs the electron-hole tight-binding Hamiltonian matrix.

    Parameters
    ----------
    tb_dims : Tuple[int, int]
        Dimensions of the tight-binding model grid.
    tb_config : List[Tuple[str, str, str]]
        Configuration of tight-binding connections.
    tb_basis : List[str]
        List of basis states.
    tb_params : Dict[str, Dict[str, float]]
        Dictionary containing electron, hole, and optionally exciton tight-binding parameters.
    tb_basis_sites_dict : Dict[str, str]
        Dictionary mapping the TB basis to the TB sites.

    Returns
    -------
    np.ndarray
        The electron-hole tight-binding Hamiltonian matrix.

    """

    matrix_electron = tb_ham_1P(
        tb_dims, tb_config, tb_basis, tb_params["electron"], tb_basis_sites_dict
    )
    matrix_hole = tb_ham_1P(tb_dims, tb_config, tb_basis, tb_params["hole"], tb_basis_sites_dict)

    dim = matrix_hole.shape[0]
    matrix = np.zeros((dim**2, dim**2))

    # exciton matrix only if tb_params["exciton"] exists.
    if "exciton" in tb_params:
        matrix_exciton = tb_ham_1P(
            tb_dims, tb_config, tb_basis, tb_params["exciton"], tb_basis_sites_dict
        )

        # exciton matrix
        if not np.allclose(matrix_exciton, np.zeros((dim, dim))):
            for i, j in product(range(dim), repeat=2):
                basis_matrix = np.zeros((dim, dim))
                basis_matrix[i, j] = 1
                matrix += matrix_exciton[i, j] * np.kron(basis_matrix, basis_matrix)

    num_sites = tb_dims[0] * tb_dims[1]
    matrix += np.kron(np.eye(num_sites), matrix_hole)
    matrix += np.kron(matrix_electron, np.eye(num_sites))
    return matrix


def add_groundstate(matrix):
    """Adds a dimension to the matrix to include the ground state.

    Parameters
    ----------
    matrix : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Matrix with an additional dimension for the ground state.

    Examples
    --------
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
    """Removes the ground state dimension from the matrix.

    Parameters
    ----------
    matrix : np.ndarray
        Input matrix with ground state dimension.

    Returns
    -------
    np.ndarray
        Matrix without the ground state dimension.

    Examples
    --------
    >>> delete_groundstate(np.array([[0., 0., 0.], [0., 1., 2.], [0., 3., 4.]]))
    array([[1., 2.],
           [3., 4.]])
    """

    return matrix[1:, 1:]


def add_interaction(
    matrix,
    eh_basis,
    interaction_param,
    interaction_type,
    nn_cutoff=False,
):
    """Adds interaction terms to the Hamiltonian based on the distance between electron
    and hole.

    Parameters
    ----------
    matrix : np.ndarray
        The initial Hamiltonian matrix.
    eh_basis : List[Tuple[str, str]]
        List of electron and hole positions as tuples of strings.
    interaction_param : float
        The interaction parameter.
    interaction_type : str
        The type of interaction. Either 'Coulomb' or 'Exchange'.
    nn_cutoff : bool, optional
        If True, only nearest neighbor interactions are considered.

    Returns
    -------
    np.ndarray
        Hamiltonian matrix with interaction terms added.

    Notes
    -----
    .. note::

        This works only for a Hamiltonian without the additional basis element accounting for relaxation.
        Therefore the interaction should always be added before the relaxation.

    Examples
    --------
    >>> Hamiltonian = np.array([[0, 1], [1, 0]])
    >>> eh_basis = [("(0, 0)", "(1, 1)"), ("(1, 0)", "(0, 0)")]
    >>> add_interaction(Hamiltonian, eh_basis, 1.0, "Coulomb", True)
    array([[0.        , 1.17639077],
           [1.17639077, 0.        ]])
    """

    distance_list = get_eh_distance(eh_basis)
    assert interaction_type in [
        "Coulomb",
        "Exchange",
    ], "Interaction type not supported."

    # as used by Bittner
    interaction_strength_list = []
    if interaction_type == "Coulomb":
        interaction_strength_list = interaction_param / (1 + 3.4 / 1 * distance_list)
    elif interaction_type == "Exchange":
        interaction_strength_list = interaction_param * np.exp(-3.4 / 0.5 * distance_list)

    # nearest neighbor cutoff
    if nn_cutoff:
        for i in np.where(distance_list > 1):
            interaction_strength_list[i] = 0

    # add interaction terms on the diagonal for exciton states
    for eh_basis_state_idx, interaction_strength in enumerate(interaction_strength_list):
        eh_basis_state = eh_basis[eh_basis_state_idx]
        matrix = set_matrix_element(
            matrix, interaction_strength, eh_basis_state, eh_basis_state, eh_basis
        )
    return matrix


# ----------------------------------------------------------------------
