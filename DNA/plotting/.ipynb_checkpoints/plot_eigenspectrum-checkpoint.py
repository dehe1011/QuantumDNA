import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from DNA.dynamics import get_reduced_dm_eigs
from DNA.model import TBHamType, PARTICLES
from utils import get_conversion, save_fig

COLORS = dict( zip(PARTICLES, sns.color_palette()[:3]) ) 

# -------------------------------------------------------------------------------------
    
def plot_eigv(ax, tb_ham: TBHamType, energy_unit: str = '100meV'):
    """
    Plots the eigenenergies.
    """
    eigv = tb_ham.eigv * get_conversion(tb_ham.unit, energy_unit)
    ax.set_title("Eigenvalues")
    ax.plot(eigv, 'o')
    ax.set_ylabel('Energy in ' + energy_unit )

def plot_eigs(ax, tb_ham: TBHamType, eigenstate_idx: int):
    """
    Plots the distribution of eigenstates over tight-binding sites. 
    """
    description = tb_ham.particle
    if description == 'exciton':
        for particle in PARTICLES:
            reduced_dm_eigs = get_reduced_dm_eigs(tb_ham, particle, eigenstate_idx)
            eigs_distribution = np.diag( reduced_dm_eigs )
            eigs_prob_distribution = np.multiply( eigs_distribution, eigs_distribution.conj() ).real
            ax.plot( tb_ham.tb_basis, eigs_prob_distribution, label=particle, color=COLORS[particle])
    if description in ('electron', 'hole'):
        particle = description
        eigs_distribution = np.diag( tb_ham.eigs )
        eigs_prob_distribution = np.multiply( eigs_distribution, eigs_distribution.conj() ).real
        ax.plot( tb_ham.tb_basis, eigs_prob_distribution, label=particle, color=COLORS[particle])
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Distribution of Eigenstate {eigenstate_idx}")
    ax.legend()
    
    
    
    