from scipy import constants
import numpy as np

from evolution import DYNAMICS_MODELS, get_initial_dm, get_eq_state, time_evo_heom, time_evo_lindblad
from hamiltonian import system_hamiltonian, system_bath_coupling_op
from lindbladian import LINDBLAD_MODELS, hamiltonian_superop, lindbladian_superop

# LINDBLAD_MODELS 
    # Local Lindblad: deph_rate
    # Local and global thermalizing Lindblad: deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent
    # HEOM: bath_cutoff, matsubara_freqs, matsubara_coeffs, matsubara_terms

init_site_pop = [1]
alpha_beta = (0,-15.5) 
epsi_delta = (0, -31)
time_interval = 5 # fs
timesteps = 500
deph_rate = 11
temperature = 77.0
cutoff_freq = 6.024
reorg_energy = 11.87 
spectral_density = 'debye'
matsubara_terms = 2
matsubara_coeffs = None # [184.80661935-100.j  65.28935061  +0.j]
matsubara_freqs = None # [10.         63.33992887]
bath_cutoff = 20

class QuantumSystem:
    def __init__(self, sites, interaction_model, dynamics_model, **settings):
        self.sites = sites 
        self.interaction_model = interaction_model 
        self.dynamics_model = dynamics_model
        # Initial state 
        self.init_site_pop = settings.get('init_site_pop', [1]) #
        # Hamiltonian 
        self.epsi_delta = settings.get('epsi_delta', (0, -31) ) 
        # Time
        self.time_interval = settings.get('time_interval', 5)
        self.timesteps = settings.get('timesteps', 500) 
        # Bath: Spectral density
        self.deph_rate = settings.get('deph_rate', 11) 
        self.temperature = settings.get('temperature', 500) 
        self.cutoff_freq = settings.get('cutoff_freq', 6.024) 
        self.reorg_energy = settings.get('reorg_energy', 11.87) 
        self.spectral_density = settings.get('spectral_density', 'debye') 
        self.ohmic_exponent = settings.get('ohmic_exponent', 1) 
        # Bath: HEOM
        self.matsubara_terms = settings.get('matsubara_terms', 2) 
        self.matsubara_coeffs = settings.get('matsubara_coeffs',None)
        self.matsubara_freqs = settings.get('matsubara_freqs', None) 
        self.bath_cutoff = settings.get('bath_cutoff', 20) 

        self.initial_density_matrix = get_initial_dm(self.sites, self.init_site_pop)
        self.hamiltonian = system_hamiltonian(self.interaction_model, self.epsi_delta)
        self.equilibrium_state = get_eq_state(self.dynamics_model, self.sites, self.hamiltonian, self.temperature)
        self.coupling_op = system_bath_coupling_op() # so far only sigma_z coupling possible
        self.hamiltonian_superop = hamiltonian_superop(self.hamiltonian)
        self.lindbladian_superop = lindbladian_superop(self.sites, self.dynamics_model, self.hamiltonian, self.deph_rate, self.cutoff_freq, 
                                                    self.reorg_energy, self.temperature, self.spectral_density, self.ohmic_exponent)  

    @property
    def time_evolution(self):
        if self.dynamics_model in LINDBLAD_MODELS:
            superop = self.hamiltonian_superop + self.lindbladian_superop
            return time_evo_lindblad(self.initial_density_matrix, superop, self.timesteps,self.time_interval, 
                                     self.dynamics_model, self.hamiltonian, self.temperature)

        if self.dynamics_model == 'HEOM':
            temperature = (self.temperature * 1e-12 * (constants.k / constants.hbar))  # convert K to rad/ps
            tmp = time_evo_heom(self.initial_density_matrix, self.timesteps, self.time_interval * 1e-3, self.hamiltonian, self.coupling_op, 
                                    self.reorg_energy, temperature,self.bath_cutoff, self.matsubara_terms, self.cutoff_freq, self.matsubara_coeffs, 
                                    self.matsubara_freqs )
            evolution, self.matsubara_coeffs, self.matsubara_freqs = tmp
            return evolution
