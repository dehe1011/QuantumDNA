import math
from scipy import constants as c
import numpy as np

# --------------------------------- spectral densities ----------------------------------------
SPECTRAL_DENSITIES = ['debye', 'ohmic', 'renger-marcus']

def debye_spectral_density(omega, cutoff_freq, reorg_energy):
    return 2 * reorg_energy * omega * cutoff_freq / (omega**2 + cutoff_freq**2)

def ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent):
    # the exponent specifies a sub- or superohmic bath 
    return (np.pi * reorg_energy * omega / cutoff_freq) * np.exp(- omega / cutoff_freq)

def renger_marcus_spectral_density(omega, reorg_energy):
    # Parameters as given in the paper
    s1, s2 = 0.8, 0.5
    hbar_w1, hbar_w2 = 0.069, 0.24  # meV
    # Convert frequencies into rad ps^-1
    w1 = hbar_w1 * 1e-15 * c.e / c.hbar  # rad ps^-1
    w2 = hbar_w2 * 1e-15 * c.e / c.hbar

    spec_dens = 0
    for si, wi in zip([s1, s2], [w1, w2]):
        tmp = si * omega**3 * np.exp(- np.sqrt(omega / wi))
        tmp /= (math.factorial(7) * 2 * wi**4)
        spec_dens += tmp

    scaling = np.pi * reorg_energy * (42 * w1 * w2) / (s1 * w2 + s2 * w1)
    return scaling * spec_dens

# -------------------------------------- rates ------------------------------------------------

def bose_einstein_distrib(omega, temperature):
    return 1. / (np.exp(c.hbar * omega * 1e12 / (c.k * temperature)) - 1)

def rate_constant_redfield(omega, deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent = 1):
    # returns eq.(1.25)
    if omega == 0:
        # Redfield rate constant at zero frequency is the dephasing rate
        if deph_rate is None:
            deph_rate = dephasing_rate(cutoff_freq, reorg_energy, temperature)
        return deph_rate

    if spectral_density == 'debye':
        spec_omega_ij = debye_spectral_density(omega, cutoff_freq, reorg_energy)
        spec_omega_ji = debye_spectral_density(-omega, cutoff_freq, reorg_energy)
    elif spectral_density == 'ohmic':
        spec_omega_ij = ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent)
        spec_omega_ji = ohmic_spectral_density(-omega, cutoff_freq, reorg_energy, exponent)
    elif spectral_density == 'renger-marcus':
        spec_omega_ij = renger_marcus_spectral_density(omega, reorg_energy)
        spec_omega_ji = renger_marcus_spectral_density(-omega, reorg_energy)

    n_omega_ij = bose_einstein_distrib(omega, temperature)
    n_omega_ji = bose_einstein_distrib(-omega, temperature)
    return 2 * ((spec_omega_ij * (1 + n_omega_ij)) + (spec_omega_ji * n_omega_ji) )

def dephasing_rate(cutoff_freq, reorg_energy, temperature):
    # limit of the Redfield rate equation as the frequency approaches zero
    return (4 * reorg_energy * c.k * temperature/ (c.hbar * cutoff_freq * 1e12))

# Try the following code to see this:
# from sympy import *
# hbar, w, wc, lam, k, T = symbols('hbar w w_c \\lambda k T')
# limit(2 * ((1 + (1 / (exp(hbar * w / (k * T)) - 1)))
#                    * (2 * lam * w * wc / (w**2 + wc**2))
#                    + (1 / (exp(hbar * - w / (k * T))))
#                    * (2 * lam * -w * wc / ((-w)**2 + wc**2))),
#               w, 0, '+')

