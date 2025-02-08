import numpy as np
import scipy.constants as c

from .save_load import PARAMETRIZATION as PARAM

# ----------------------------------------


def calc_orbital_energy(orbital):
    """
    Calculate the energy of a given orbital.

    Parameters
    ----------
    orbital : Orbital
        An instance of the Orbital class containing information about the atom and orbital type.

    Returns
    -------
    float
        The energy of the specified orbital.
    """

    return PARAM["E_" + orbital.atom + orbital.orbital_type[0]]


def calc_orbital_overlap(
    orbital1, orbital2, connection
):  # pylint: disable=too-many-branches, too-many-statements
    """
    Calculate the overlap between two orbitals using the Harrison expression and the Slater-Koster two-center interactions.

    Parameters
    ----------
    orbital1 : Orbital
        The first orbital object containing its coordinates, atom type, orbital type, and atom identifier.
    orbital2 : Orbital
        The second orbital object containing its coordinates, atom type, orbital type, and atom identifier.
    connection : str
        The type of connection between the orbitals, either "interbase" or "intrabase".

    Returns
    -------
    float
        The calculated overlap between the two orbitals. Returns 0 if the distance between orbitals is zero,
        if the orbitals belong to the same atom in an intrabase connection, or if the distance exceeds the cutoff radius.

    Notes
    -----
    .. note::

        The directional cosines are used as projections on the coordinate axis rather than actual cosine values.
        The function includes a hydrogen correction factor obtained by optimization.
    """

    vector = np.array(orbital2.coordinates) - np.array(orbital1.coordinates)
    distance = np.linalg.norm(vector)

    if distance == 0:
        return 0

    vector_norm = vector / distance

    orbital_types = [orbital1.orbital_type, orbital2.orbital_type]

    # ----------------------------------------

    prefactor = None

    # interbase overlap, exponential decay term
    if connection == "interbase":
        if (orbital1.atom == "H") ^ (orbital2.atom == "H"):
            exponent = np.e ** (-2 / PARAM["d0H"] * (distance - PARAM["d0H"]))
            prefactor = (
                c.hbar**2 / (c.m_e * c.angstrom**2 * c.e * PARAM["d0H"] ** 2) * exponent
            )

        else:
            exponent = np.e ** (-2 / PARAM["d0"] * (distance - PARAM["d0"]))
            prefactor = (
                c.hbar**2 / (c.m_e * c.angstrom**2 * c.e * PARAM["d0"] ** 2) * exponent
            )

    # intrabase overlap, cutoff radius
    if connection == "intrabase":
        # exclude overlap between orbitals of the same atom
        if orbital1.atom_identifier == orbital2.atom_identifier:
            return 0

        # exclude overlap between orbitals that are too far apart
        if distance >= PARAM["cutoff_radius"]:
            return 0

        prefactor = c.hbar**2 / (c.m_e * c.angstrom**2 * c.e * distance**2)

    assert prefactor is not None, "Invalid connection type."

    # ----------------------------------------

    # hydrogen correction factor obtained by optimization
    correction = 1
    if orbital1.atom == "H":
        correction *= PARAM["b"]
    if orbital2.atom == "H":
        correction *= PARAM["b"]

    # ----------------------------------------

    overlap = None

    V_sssigma = PARAM["chi_sssigma"]
    V_spsigma = PARAM["chi_spsigma"]
    V_ppsigma = PARAM["chi_ppsigma"]
    V_pppi = PARAM["chi_pppi"]

    if orbital_types == ["s", "s"]:
        overlap = V_sssigma

    if orbital_types == ["s", "px"]:
        eta_1 = -vector_norm[0]
        overlap = V_spsigma * eta_1

    if orbital_types == ["px", "s"]:
        eta_1 = vector_norm[0]
        overlap = V_spsigma * eta_1

    if orbital_types == ["s", "py"]:
        eta_1 = -vector_norm[1]
        overlap = V_spsigma * eta_1

    if orbital_types == ["py", "s"]:
        eta_1 = vector_norm[1]
        overlap = V_spsigma * eta_1

    if orbital_types == ["s", "pz"]:
        eta_1 = -vector_norm[2]
        overlap = V_spsigma * eta_1

    if orbital_types == ["pz", "s"]:
        eta_1 = vector_norm[2]
        overlap = V_spsigma * eta_1

    # --------------------------------------

    if orbital_types in [["px", "py"], ["py", "px"]]:
        eta_1 = vector_norm[0]
        eta_2 = vector_norm[1]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    if orbital_types in [["py", "pz"], ["pz", "py"]]:
        eta_1 = vector_norm[1]
        eta_2 = vector_norm[2]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    if orbital_types in [["px", "pz"], ["pz", "px"]]:
        eta_1 = vector_norm[0]
        eta_2 = vector_norm[2]
        overlap = eta_1 * eta_2 * (V_ppsigma - V_pppi)

    if orbital_types == ["px", "px"]:
        eta_1 = vector_norm[0]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    if orbital_types == ["py", "py"]:
        eta_1 = vector_norm[1]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    if orbital_types == ["pz", "pz"]:
        eta_1 = vector_norm[2]
        overlap = eta_1**2 * V_ppsigma + (1 - eta_1**2) * V_pppi

    assert overlap is not None, "Invalid orbital types."

    return correction * prefactor * overlap
