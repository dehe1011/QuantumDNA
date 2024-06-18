import numpy as np
from typing import Union, List, Tuple, Any
from itertools import product
from ast import literal_eval

# -------------------------------------------- Tight-binding basis ----------------------------------------------

def get_tb_basis(tb_dims: Tuple[int, int]) -> List[str]:
    """
    Generates a symbolic string representation of the indices '(i, j)' where i indicates the strand number and j the site number inside the strand.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.

    Returns:
        List of string representations of site indices.

    Examples:
        >>> get_tb_basis((2, 2))
        ['(0, 0)', '(0, 1)', '(1, 0)', '(1, 1)']
    """
    num_sites = tb_dims[0] * tb_dims[1]
    tb_basis = [int_to_str(tb_dims, site_num) for site_num in range(num_sites)]
    return tb_basis

def str_to_int(tb_dims: Tuple[int, int], tb_basis_str: str) -> int:
    """
    Converts a string representation of a site index to an integer representation.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        tb_basis_str: String representation of the site index.

    Returns:
        Integer representation of the site index.

    Example:
        >>> str_to_int( (1,3), '(0, 2)')
        2
    """
    return tuple_to_int(tb_dims, str_to_tuple(tb_basis_str))

def int_to_str(tb_dims: Tuple[int, int], tb_basis_int: int) -> str:
    """
    Converts an integer representation of a site index to a string representation.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        tb_basis_num: Integer representation of the site index.

    Returns:
        String representation of the site index.

    Example:
        >>> int_to_str( (1,3), 2 )
        '(0, 2)'
    """
    return tuple_to_str(int_to_tuple(tb_dims, tb_basis_int))

def str_to_tuple(tb_basis_str: str) -> Tuple[int, int]:
    """
    Converts a string representation of a site index to a tuple index.

    Args:
        tb_basis_str: String representation of the site index.

    Returns:
        Tuple representation of the site index.
    """
    return literal_eval(tb_basis_str)

def tuple_to_str(tb_basis_tuple: Tuple[int, int]) -> str:
    """
    Converts a tuple representation of a site index to a string representation.

    Args:
        tb_basis_idx: Tuple representation of the site index.

    Returns:
        String representation of the site index.
    """
    return str(tb_basis_tuple)

def int_to_tuple(tb_dims: Tuple[int, int], tb_basis_int: int) -> Tuple[int, int]:
    """
    Converts an integer representation of a site index to a tuple index.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        tb_basis_num: Integer representation of the site index.

    Returns:
        Tuple representation of the site index.
    """
    num_strands, num_sites_per_strand = tb_dims
    return (tb_basis_int // num_sites_per_strand, tb_basis_int % num_sites_per_strand)

def tuple_to_int(tb_dims: Tuple[int, int], tb_basis_tuple: Tuple[int, int]) -> int:
    """
    Converts a tuple representation of a site index to an integer representation.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        tb_basis_idx: Tuple representation of the site index.

    Returns:
        Integer representation of the site index.
    """
    num_strands, num_sites_per_strand = tb_dims
    strand_num, site_in_strand_num = tb_basis_tuple
    return strand_num * num_sites_per_strand + site_in_strand_num

# --------------------------------------- Electron-hole basis --------------------------------------------------

def get_eh_basis(tb_dims: Tuple[int, int]) -> List[Tuple[str, str]]:
    """
    Generates a symbolic string representation of electron-hole basis indices.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.

    Returns:
        List of tuples representing electron-hole basis indices.

    Examples:
        >>> get_eh_basis((2, 2))
        [('(0, 0)', '(0, 0)'), ('(0, 0)', '(0, 1)'), ('(0, 0)', '(1, 0)'), ('(0, 0)', '(1, 1)'),
         ('(0, 1)', '(0, 0)'), ('(0, 1)', '(0, 1)'), ('(0, 1)', '(1, 0)'), ('(0, 1)', '(1, 1)'),
         ('(1, 0)', '(0, 0)'), ('(1, 0)', '(0, 1)'), ('(1, 0)', '(1, 0)'), ('(1, 0)', '(1, 1)'),
         ('(1, 1)', '(0, 0)'), ('(1, 1)', '(0, 1)'), ('(1, 1)', '(1, 0)'), ('(1, 1)', '(1, 1)')]
    """
    tb_basis = get_tb_basis(tb_dims)
    eh_basis = list(product(tb_basis, tb_basis))
    return eh_basis

def get_particle_eh_states(particle: str, tb_basis_element: str, tb_basis: List[str]) -> List[Tuple[str, str]]:
    """
    Generates a list of electron-hole states for a given particle and site basis element.

    Args:
        particle: The type of particle ('electron', 'hole', or 'exciton').
        tb_basis_element: The site basis element for which to generate the states.
        tb_basis: The list of all site basis elements.

    Returns:
        List of tuples representing the electron-hole states.

    Examples:
        >>> get_particle_eh_states('electron', '(0, 2)', get_tb_basis((1,3)))
        [('(0, 2)', '(0, 0)'), ('(0, 2)', '(0, 1)'), ('(0, 2)', '(0, 2)')]
    """
    if tb_basis_element not in tb_basis:
        raise ValueError(f"Given basis element {tb_basis_element} not in tight-binding basis {tb_basis}")
    
    if particle == 'electron':
        return list(product([tb_basis_element], tb_basis))
    elif particle == 'hole':
        return list(product(tb_basis, [tb_basis_element]))
    elif particle == 'exciton':
        return list(product([tb_basis_element], [tb_basis_element]))
    else:
        raise ValueError(f"Invalid particle type: {particle}")

# --------------------------------------- Basis conversion -----------------------------------------

def basis_converter(basis_element: Any, basis: List[Any]) -> Any:
    """
    Converts between basis index representations (integer to string or vice versa).

    Args:
        basis_element: The basis element to convert (can be an integer, string, or tuple).
        basis: The list of basis elements for reference.

    Returns:
        Corresponding basis element in the other representation.

    Example: 
        >>> basis_converter( 'e2', ['e1', 'e2', 'e3'])
        1
    """
    if isinstance(basis_element, int):
        return basis[basis_element]
    else:
        return basis.index(basis_element)

def tb_basis_converter(tb_dims: Tuple[int, int], tb_basis_element: Union[str, int]) -> Union[int, str]:
    """
    Converts between site basis index representations.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        tb_basis_element: The basis element to convert (can be an integer or string).

    Returns:
        Corresponding basis element in the other representation.

    Examples:
        >>> tb_basis_converter((2, 3), '(1, 2)')
        5
        >>> tb_basis_converter((2, 3), 5)
        '(1, 2)'
    """
    return basis_converter(tb_basis_element, get_tb_basis(tb_dims))

def eh_basis_converter(tb_dims: Tuple[int, int], eh_basis_element: Union[int, Tuple[str, str]]) -> Union[Tuple[str, str], int]:
    """
    Converts between electron-hole basis index representations.

    Args:
        tb_dims: A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
        eh_basis_element: The basis element to convert (can be an integer or tuple of strings).

    Returns:
        Corresponding basis element in the other representation.

    Examples:
        >>> eh_basis_converter((2, 3), ('(0, 2)', '(1, 1)'))
        16
        >>> eh_basis_converter((2, 3), 16)
        ('(0, 2)', '(1, 1)')
    """
    return basis_converter(eh_basis_element, get_eh_basis(tb_dims))

# -------------------------- Basis change from local to global basis (eigenbasis) -------------------------------

def basis_change(
    matrix: np.ndarray, states: np.ndarray, liouville: bool = False
) -> np.ndarray:
    """
    Performs a basis change of the given matrix.

    Args:
        matrix: The matrix to change basis.
        states: The old basis expressed as a vector in the new basis.
        liouville: Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns:
        Matrix in the new basis.
    """
    if liouville:
        states = np.kron(states, states.conjugate())
    return np.matmul(states, np.matmul(matrix, states.conjugate().T))

def global_to_local(
    matrix: np.ndarray, eigs: np.ndarray, liouville: bool = False
) -> np.ndarray:
    """
    Performs a basis change from the eigenbasis to the site basis.

    Args:
        matrix: The matrix to change basis.
        eigs: The eigenbasis (old) expressed in the site basis (new).
        liouville: Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns:
        Matrix in the site basis.
    """
    return basis_change(matrix, eigs, liouville=liouville)

def local_to_global(
    matrix: np.ndarray, eigs: np.ndarray, liouville: bool = False
) -> np.ndarray:
    """
    Performs a basis change from the site basis to the eigenbasis.

    Args:
        matrix: The matrix to change basis.
        eigs: The site basis (old) expressed in the eigenbasis (new).
        liouville: Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns:
        Matrix in the eigenbasis.
    """
    return basis_change(matrix, eigs.conjugate().T, liouville=liouville)
