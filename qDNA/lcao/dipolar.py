import numpy as np
import scipy.constants as c

# ----------------------------------------------------------------------


def calc_dipolar_coupling(monomer1, monomer2):
    """
    Calculate the dipolar coupling between two monomers.

    Parameters
    ----------
    monomer1 : object
        An object representing the first monomer
    monomer2 : object
        An object representing the second monomer

    Returns
    -------
    float
        The dipolar coupling constant in units of eV.
    """

    R1, R2 = monomer1.center_of_mass, monomer2.center_of_mass
    if np.allclose(R1, R2):
        return 0
    mu1, mu2 = monomer1.dipole_moment, monomer2.dipole_moment
    d = np.linalg.norm(R1 - R2)
    prefactor = 1 / (4 * np.pi * c.epsilon_0) * 1 / d**3
    factor = mu1 @ mu2 - 3 / d**2 * (mu1 @ (R2 - R1)) * (mu2 @ (R2 - R1))
    return prefactor * factor / (c.angstrom * c.e)


# ----------------------------------------------------------------------
