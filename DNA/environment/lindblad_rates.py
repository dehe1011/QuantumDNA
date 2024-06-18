import numpy as np
import scipy.constants as c

__all__ = ['debye_spectral_density', 'ohmic_spectral_density', 'bose_einstein_distrib', 'rate_constant_redfield', 'dephasing_rate' ]

# --------------------------- bath spectral densities --------------------------------------

def debye_spectral_density(omega, cutoff_freq, reorg_energy):
    if omega <= 0: return 0
    return 2 * reorg_energy * omega * cutoff_freq / (omega**2 + cutoff_freq**2)

def ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent):
    # the exponent specifies a sub- or superohmic bath 
    if omega <= 0: return 0
    return (np.pi * reorg_energy * omega / cutoff_freq) * np.exp(- omega / cutoff_freq)

# ----------------------------- Lindblad rates -------------------------------------------

def bose_einstein_distrib(omega, temperature):
    return 1. / (np.exp(c.hbar * omega * 1e12 / (c.k * temperature)) - 1)

def rate_constant_redfield(omega, deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent):
    # returns eq.(1.25) in Abbott     
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
