from itertools import combinations

import numpy as np
import scipy.constants as c

# ----------------------------------------------------------------------


def get_prefactor(lcao_param, vector, orbital_atoms):
    """
    Calculate the prefactor for Slater-Koster matrix elements.

    Parameters
    ----------
    lcao_param : dict
        Dictionary containing Slater-Koster parameters including 'd0', 'd0H', 'cutoff_radius',
        and 'b' (hydrogen correction factor).
    vector : np.array
        Vector connecting the two atoms in Cartesian coordinates (x, y, z).
    orbital_atoms : list or str
        List or string of atom symbols involved in the orbital connection.

    Returns
    -------
    float
        The calculated prefactor for the Slater-Koster matrix element, including hydrogen
        correction factor if needed.

    Notes
    -----
    .. note::
        - For "interbase" connections, an exponential decay factor is applied. Based on
        hydrogen involvement, either d0H or d0 is used as the characteristic distance.
        - For "intrabase" connections, a simple 1/d^2 Harrison-type dependence is used.

    Examples
    --------
    >>> lcao_param = {'d0': 2.0, 'd0H': 1.5, 'b': 0.8, 'cutoff_radius': 10.0}
    >>> vector = np.array([1.0, 0.0, 0.0])
    >>> atoms = ['C', 'C']
    >>> prefactor = get_prefactor(lcao_param, vector, atoms)
    >>> print(prefactor)
    """

    distance = np.linalg.norm(vector)

    # Constants for unit conversion
    unit_factor = c.hbar**2 / (c.m_e * c.angstrom**2 * c.e)

    # Determine d0 and decay constant depending on hydrogen involvement
    is_hydrogen = (orbital_atoms[0] == "H") or (orbital_atoms[1] == "H")
    d0 = lcao_param["d0"] if not is_hydrogen else lcao_param["d0H"]
    exponent = np.exp(-2 / d0 * (distance - d0))

    if distance >= d0:
        prefactor = unit_factor / (d0**2) * exponent
    else:
        prefactor = unit_factor / (distance**2)

    # Apply hydrogen correction factor if needed
    correction = 1.0
    if "H" in orbital_atoms:
        correction *= lcao_param["b"] ** orbital_atoms.count("H")

    return prefactor * correction


# pylint: disable=too-many-branches, too-many-statements
def get_overlap(lcao_param, vector, orbital_types):
    """
    Calculate the interaction integral between two orbitals using Slater-Koster rules.
    This function computes the interaction between two atomic orbitals based on their
    orbital types (s, px, py, pz) and the vector connecting the two atoms using the
    provided  Slater-Koster parameterization. Directional cosines are given as the
    components of the normalized vector.

    Parameters
    ----------
    lcao_param : dict
        Dictionary containing Slater-Koster parameters.
    vector : np.array
        Vector connecting the two atoms in Cartesian coordinates (x, y, z).
    orbital_types : list of str
        List containing two orbital type designators, e.g., ["s", "px"].
        Valid orbital types are: "s", "px", "py", "pz".

    Returns
    -------
    float
        The overlap integral between the two orbitals.

    Examples
    --------
    >>> lcao_param = {'chi_sssigma': 0.5, 'chi_spsigma': 0.3, 'chi_ppsigma': 0.2, 'chi_pppi': 0.1}
    >>> vector = np.array([1.0, 0.0, 0.0])
    >>> overlap = get_overlap(lcao_param, vector, ["s", "s"])
    >>> print(overlap)
    """

    distance = np.linalg.norm(vector)
    vector_norm = vector / distance

    V_sssigma = lcao_param["chi_sssigma"]
    V_spsigma = lcao_param["chi_spsigma"]
    V_ppsigma = lcao_param["chi_ppsigma"]
    V_pppi = lcao_param["chi_pppi"]

    if orbital_types == ["s", "s"]:
        overlap = V_sssigma

    elif orbital_types == ["s", "px"]:
        eta_1 = -vector_norm[0]
        overlap = V_spsigma * eta_1

    elif orbital_types == ["px", "s"]:
        eta_1 = vector_norm[0]
        overlap = V_spsigma * eta_1

    elif orbital_types == ["s", "py"]:
        eta_1 = -vector_norm[1]
        overlap = V_spsigma * eta_1

    elif orbital_types == ["py", "s"]:
        eta_1 = vector_norm[1]
        overlap = V_spsigma * eta_1

    elif orbital_types == ["s", "pz"]:
        eta_1 = -vector_norm[2]
        overlap = V_spsigma * eta_1

    elif orbital_types == ["pz", "s"]:
        eta_1 = vector_norm[2]
        overlap = V_spsigma * eta_1

    elif orbital_types in [["px", "py"], ["py", "px"]]:
        eta_1 = vector_norm[0]
        eta_2 = vector_norm[1]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    elif orbital_types in [["py", "pz"], ["pz", "py"]]:
        eta_1 = vector_norm[1]
        eta_2 = vector_norm[2]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    elif orbital_types in [["px", "pz"], ["pz", "px"]]:
        eta_1 = vector_norm[0]
        eta_2 = vector_norm[2]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    elif orbital_types == ["px", "px"]:
        eta_1 = vector_norm[0]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    elif orbital_types == ["py", "py"]:
        eta_1 = vector_norm[1]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    elif orbital_types == ["pz", "pz"]:
        eta_1 = vector_norm[2]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    else:
        overlap = None
        raise ValueError(f"Unknown orbital types: {orbital_types}")

    return overlap


def calc_orbital_interaction(lcao_param, orbitals, orbitals_coordinates):
    """
    Calculate the interaction between two orbitals based on LCAO parameters.

    Parameters
    ----------
    lcao_param : dict
        Dictionary containing the LCAO parametrization.
    orbitals : list of str
        List of the two orbital identifiers in the format "<atom>_<orbital_type>", e.g., ["C_s", "N_px"].
    orbitals_coordinates : ndarray
        Array of shape (2, 3) representing the coordinates of the two orbitals.

    Returns
    -------
    float
        The calculated orbital interaction value. Returns 0 if the distance between orbitals is zero.

    Examples
    --------
    >>> lcao_param = {'d0': 2.0, 'd0H': 1.5, 'b': 0.8, 'cutoff_radius': 10.0, 'chi_sssigma': 0.5, 'chi_spsigma': 0.3, 'chi_ppsigma': 0.2, 'chi_pppi': 0.1}
    >>> orbitals = ["C_s", "N_px"]
    >>> orbitals_coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    >>> interaction = calc_orbital_interaction(lcao_param, orbitals, orbitals_coordinates)
    >>> print(interaction)
    """

    orbital_atoms = [orbital.split("_")[0] for orbital in orbitals]
    orbital_types = [orbital.split("_")[1] for orbital in orbitals]

    vector = orbitals_coordinates[1] - orbitals_coordinates[0]
    distance = np.linalg.norm(vector)
    if distance == 0:
        return 0

    prefactor = get_prefactor(lcao_param, vector, orbital_atoms)
    overlap = get_overlap(lcao_param, vector, orbital_types)
    return prefactor * overlap


def calc_orbital_energy(lcao_param, orbital):
    """
    Calculate the energy of a specified orbital using LCAO parameters.

    Parameters
    ----------
    lcao_param : dict
        Dictionary containing LCAO parameters.
        Expected keys are in the format "E_<atom><orbital_type>".
    orbital : str
        Orbital identifier in the format "<atom>_<orbital_type>".

    Returns
    -------
    float
        Energy of the specified orbital.

    Examples
    --------
    >>> lcao_param = {'E_Cs': -10.0, 'E_Cpx': -5.0, 'E_Ns': -12.0, 'E_Npx': -6.0}
    >>> orbital = "C_s"
    >>> energy = calc_orbital_energy(lcao_param, orbital)
    >>> print(energy)
    """

    orbital_atom, orbital_type = orbital.split("_")
    return lcao_param["E_" + orbital_atom + orbital_type[0]]


def calc_H_intra(lcao_param, comp):
    """
    Calculate the intra-base Hamiltonian matrix for a given component.

    Parameters
    ----------
    lcao_param : dict
        Parameters for the Linear Combination of Atomic Orbitals (LCAO) model.
    comp : object
        Component containing orbital information and coordinates.

    Returns
    -------
    np.ndarray
        Intra-base Hamiltonian matrix of shape (n, n), where `n` is the number
        of orbitals in the component.
    """

    n = comp.num_orbitals
    H_intra = np.zeros((n, n))

    # orbital overlap
    for i, j in combinations(range(n), r=2):
        orbitals = [comp.orbitals[i], comp.orbitals[j]]
        coords = [comp.orbitals_coordinates[i], comp.orbitals_coordinates[j]]
        value = calc_orbital_interaction(lcao_param, orbitals, coords)
        H_intra[i, j] = value
        H_intra[j, i] = value

    # orbital energy
    for i in range(n):
        orbital = comp.orbitals[i]
        value = calc_orbital_energy(lcao_param, orbital)
        H_intra[i, i] = value
    return H_intra


def calc_H_inter(lcao_param, comp1, comp2):
    """
    Calculate the interbase Hamiltonian matrix for two components.

    Parameters
    ----------
    lcao_param : dict
        Parameters for the Linear Combination of Atomic Orbitals (LCAO) method.
    comp1 : object
        First component containing orbital information and coordinates.
    comp2 : object
        Second component containing orbital information and coordinates.

    Returns
    -------
    np.ndarray
        Inter-base Hamiltonian matrix of shape (n, m), where `n` and `m` are the
        numbers of orbitals in the components.
    """

    n1, n2 = comp1.num_orbitals, comp2.num_orbitals
    H_inter = np.zeros((n1, n2), dtype=float)

    for i in range(n1):
        for j in range(n2):
            orbitals = [comp1.orbitals[i], comp2.orbitals[j]]
            coords = [comp1.orbitals_coordinates[i], comp2.orbitals_coordinates[j]]
            value = calc_orbital_interaction(lcao_param, orbitals, coords)
            H_inter[i, j] = value
    return H_inter


# ----------------------------------------------------------------------
