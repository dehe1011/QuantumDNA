from utils import load_TB_params, complementary_strand, add_basis_dimension, H_particle, H_exciton, H_eh_config, get_conversion, get_config
import qutip as q

__all__ = ['Hamiltonian', 'H_p_single_base', 'H_single_base', 'H_p_base_pair', 'H_base_pair']

# ------------------------------------------ class for TB Hamiltonians ----------------------------------------------
class Hamiltonian:
    def __init__(self, DNAstring, Ham_kwargs={}):
        self.DNAstring = DNAstring 

        self.default = get_config()['Ham_kwargs_default']
        self.source = Ham_kwargs.get('source', self.default['source'])
        self.J = Ham_kwargs.get('J', self.default['J'])
        self.V_F = Ham_kwargs.get('V_F', self.default['V_F'])
        self.diagonal_hopping = Ham_kwargs.get('diagonal_hopping', self.default['diagonal_hopping'])
        self.description = Ham_kwargs.get('description', self.default['description'])
        self.relaxation = Ham_kwargs.get('relaxation', self.default['relaxation'])
        self.particle = Ham_kwargs.get('particle', self.default['particle'])
        self.unit = Ham_kwargs.get('unit', self.default['unit']) # Hamiltonian given in units of a rate 
        self.verbose = Ham_kwargs.get('verbose', self.default['verbose'])

        if self.description == 'single_base': self.N = 2*len(self.DNAstring)
        elif self.description == 'base_pair': self.N = len(self.DNAstring)
        if self.particle != 'eh': self.relaxation = False
        if self.verbose: 
            print(f'Hamiltonian object with parameters: \n {vars(self)}')

    @property
    def matrix(self):
        if self.particle == 'eh' and self.description == 'single_base':
            H = H_single_base(self.DNAstring, source = self.source, diagonal_hopping = self.diagonal_hopping, 
                          relaxation = self.relaxation, J = self.J, V_F = self.V_F)
        if self.particle == 'eh' and self.description == 'base_pair':
            H = H_base_pair(self.DNAstring, source = self.source, relaxation = self.relaxation, J = self.J, V_F = self.V_F)
        if self.particle == ('electron' or 'hole') and self.description == 'single_base':
            H = H_p_single_base(self.DNAstring, self.particle, source = self.source, diagonal_hopping = self.diagonal_hopping)
        if self.particle == ('electron' or 'hole') and self.description == 'base_pair':
            H = H_p_base_pair(self.DNAstring, self.particle, source = self.source)
        return q.Qobj( (H * get_conversion('100meV', self.unit) ).full() )

    @matrix.setter
    def matrix(self, new_matrix):
        self._matrix = new_matrix
        
# --------------------------- single base Hamiltonian -----------------------------------------------------
    
def H_p_single_base(DNAstring, particle, source = 'Hawke2010', diagonal_hopping = True):
    description = 'single_base'
    N, lun = 2 * len(DNAstring), len(DNAstring)
    DNA = DNAstring + complementary_strand(DNAstring)
    TB_params = load_TB_params(source, particle, description, folder = 'stored_data/tb_params', load_metadata = False)
    # construction of the TB Hamiltonian 
    H=0
    for i in range(N):
        H += TB_params[DNA[i]]*H_particle(N,i,i)
    for i in range(lun-1):
        H += TB_params['t_'+DNA[i]+DNA[i+1]] * (H_particle(N,i,i+1)+H_particle(N,i+1,i))
    for i in range(lun+1,N):
        H += TB_params['t_'+DNA[i]+DNA[i-1]] * (H_particle(N,i-1,i)+H_particle(N,i,i-1))
    for i in range(lun):
        H += TB_params['h_'+DNA[i]+DNA[i+lun]] * (H_particle(N,i,i+lun)+H_particle(N,i+lun,i))
    if diagonal_hopping == True:
        for i in range(lun-1):
            H += TB_params['r+_'+DNA[i]+DNA[i+lun+1]] * (H_particle(N,i,i+lun+1)+H_particle(N,i+lun+1,i))
            H += TB_params['r-_'+DNA[i+1]+DNA[i+lun]] * (H_particle(N,i+1,i+lun)+H_particle(N,i+lun,i+1))
    return H # multiply by 100 meV/hbar to obtain a rate in units of rad/s 

def H_I_single_base(N,J):
    # neigherst neighbor interaction 
    lun=N//2
    H=0
    U=[J,J/(1+3.4),J/(1+3.4)]
    for i in range(N):
        H += U[0] * H_eh_config(N,i,i,i,i)
    for i in range(lun-1): 
        H += U[1] * (H_eh_config(N,i,i,i+1,i+1) + H_eh_config(N,i+1,i+1,i,i)) # next neighbor in top chain
        H += U[1] * (H_eh_config(N,i+lun,i+lun,i+lun+1,i+lun+1) + H_eh_config(N,i+lun+1,i+lun+1,i+lun,i+lun)) 
        # next neighbor in bottom chain
    for i in range(lun):
        H += U[2] * (H_eh_config(N,i,i,i+lun,i+lun) + H_eh_config(N,i+lun,i+lun,i,i)) # next neighbor between the chains
    return H 

def H_F_single_base(N,V_F):
    # Förster resonant energy transfer (FRET)
    # this feature is deprecated since we are not sure about the physical correctness of this implementation
    lun=N//2
    H=0
    for i in range(lun-1):
        H += V_F * (H_exciton(N,i,i+1)+H_exciton(N,i+1,i))
    for i in range(lun+1,N):
        H += V_F * (H_exciton(N,i-1,i)+H_exciton(N,i,i-1))
    for i in range(lun):
        H += V_F * (H_exciton(N,i,i+lun)+H_exciton(N,i+lun,i))
    return H
    
def H_single_base(DNAstring, source = 'Hawke2010', diagonal_hopping = True, relaxation = True, J = 0, V_F = 0):
    N=2*len(DNAstring)
    DNA = DNAstring+complementary_strand(DNAstring)
    # construction of the total TB Hamiltonian 
    H = 0
    H += q.tensor( q.qeye(N), H_p_single_base(DNAstring, 'hole', source = source, diagonal_hopping = diagonal_hopping) )
    H += q.tensor( H_p_single_base(DNAstring, 'electron', source = source, diagonal_hopping = diagonal_hopping), q.qeye(N) )
    if J!=0: H += H_I_single_base(N,J)
    if V_F!=0: H += H_F_single_base(N,V_F)
    if relaxation == True: H = add_basis_dimension(N,H)
    return H

#--------------------------------------- base pair Hamiltonian -------------------------------------------------

def H_p_base_pair(DNAstring, particle, source = 'Hawke2010'):
    description = 'base_pair'
    N=len(DNAstring)
    DNA = DNAstring
    TB_params = load_TB_params(source, particle, description, folder = 'stored_data/tb_params', load_metadata = False)
    # construction of the TB Hamiltonian 
    H=0
    for i in range(N):
        H += TB_params[DNA[i]] * H_particle(N,i,i)
    for i in range(N-1):
        H += TB_params['t_'+DNA[i]+DNA[i+1]] * (H_particle(N,i,i+1)+H_particle(N,i+1,i))
    return H 

def H_I_base_pair(N,J):
    H=0
    U=[J,J/(1+3.4)]
    for i in range(N):
        H += U[0] * H_eh_config(N,i,i,i,i)
    for i in range(N-1): 
        H += U[1] * ( H_eh_config(N,i,i,i+1,i+1) + H_eh_config(N,i+1,i+1,i,i) )
    return H 

def H_F_base_pair(N,V_F):
    H=0
    for i in range(N-1):
        H += V_F * ( H_exciton(N,i,i+1) + H_exciton(N,i+1,i) )
    return H
        
def H_base_pair(DNAstring, source = 'Hawke2010', relaxation = True, J = 0, V_F = 0):
    N=len(DNAstring)
    # construction of the total TB Hamiltonian 
    H = 0
    H += q.tensor( q.qeye(N), H_p_base_pair(DNAstring, 'hole', source = source) )
    H += q.tensor( H_p_base_pair(DNAstring, 'electron', source = source), q.qeye(N) )
    if J!=0: H += H_I_base_pair(N,J)
    if V_F!=0: H += H_F_base_pair(N,V_F)
    if relaxation == True: H = add_basis_dimension(N,H)
    return H

