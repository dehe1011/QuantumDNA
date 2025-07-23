from itertools import combinations

import numpy as np
import scipy.constants as c

# ----------------------------------------------------------------------


def _get_prefactor(lcao_param, vector, orbital_atoms, connection_type):
    distance = np.linalg.norm(vector)

    # Constants for unit conversion
    unit_factor = c.hbar**2 / (c.m_e * c.angstrom**2 * c.e)

    if connection_type == "interbase":
        # Determine d0 and decay constant depending on hydrogen involvement
        is_hydrogen = (orbital_atoms[0] == "H") ^ (orbital_atoms[1] == "H")
        d0 = lcao_param["d0H"] if is_hydrogen else lcao_param["d0"]
        exponent = np.exp(-2 / d0 * (distance - d0))
        prefactor = unit_factor / (d0**2) * exponent

    elif connection_type == "intrabase":
        # Check if orbitals are too far apart to contribute
        if distance >= lcao_param["cutoff_radius"]:
            return 0
        prefactor = unit_factor / (distance**2)

    else:
        raise ValueError(f"Unknown connection_type: {connection_type}")

    # Apply hydrogen correction factor if needed
    correction = 1.0
    if "H" in orbital_atoms:
        correction *= lcao_param["b"] ** orbital_atoms.count("H")

    return prefactor * correction


# pylint: disable=too-many-branches, too-many-statements
def _get_overlap(lcao_param, vector, orbital_types):
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


def calc_orbital_interaction(
    lcao_param, orbitals, orbitals_coordinates, connection_type
):
    """
    Calculate the interaction between two orbitals based on LCAO parameters.

    Parameters
    ----------
    lcao_param : dict
        Dictionary containing Linear Combination of Atomic Orbitals (LCAO) parameters.
    orbitals : list of str
        List of orbital identifiers in the format "atom_orbital".
    orbitals_coordinates : ndarray
        Array of shape (2, 3) representing the coordinates of the two orbitals.
    connection_type : str
        Type of connection between the orbitals (e.g., covalent, ionic).

    Returns
    -------
    float
        The calculated orbital interaction value. Returns 0 if the distance between orbitals is zero.
    """

    orbital_atoms = [orbital.split("_")[0] for orbital in orbitals]
    orbital_types = [orbital.split("_")[1] for orbital in orbitals]

    vector = orbitals_coordinates[1] - orbitals_coordinates[0]
    distance = np.linalg.norm(vector)
    if distance == 0:
        return 0

    prefactor = _get_prefactor(lcao_param, vector, orbital_atoms, connection_type)
    overlap = _get_overlap(lcao_param, vector, orbital_types)
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
        Component containing orbital information.

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
        value = calc_orbital_interaction(lcao_param, orbitals, coords, "intrabase")
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
        A 2D array representing the interbase Hamiltonian matrix.
    """

    n1, n2 = comp1.num_orbitals, comp2.num_orbitals
    H_inter = np.zeros((n1, n2), dtype=float)

    for i in range(n1):
        for j in range(n2):
            orbitals = [comp1.orbitals[i], comp2.orbitals[j]]
            coords = [comp1.orbitals_coordinates[i], comp2.orbitals_coordinates[j]]
            value = calc_orbital_interaction(lcao_param, orbitals, coords, "interbase")
            H_inter[i, j] = value
    return H_inter


# ----------------------------------------------------------------------
