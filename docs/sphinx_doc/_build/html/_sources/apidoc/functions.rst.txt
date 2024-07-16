Functions
=========

This is based on :cite:`Herb2024` .

Tight-binding basis
-------------------

.. automodule:: qDNA.model
   :members: get_tb_basis, get_eh_basis, get_eh_distance, get_particle_eh_states, basis_change, local_to_global, global_to_local
   :show-inheritance: False
   
Tight-binding configuration
---------------------------

.. automodule:: qDNA.model
   :members: get_tb_config
   :show-inheritance: False
   :noindex:

Tight-binding Hamiltonian 
-------------------------

.. automodule:: qDNA.model
   :members: set_matrix_element, tb_ham_1P, tb_ham_2P, add_groundstate, delete_groundstate, add_interaction
   :show-inheritance: False
   :noindex:

Save and load tight-binding parameters
--------------------------------------

.. automodule:: qDNA.model
   :members: save_tb_params, load_tb_params, wrap_save_tb_params, wrap_load_tb_params
   :show-inheritance: False
   :noindex:
   
   
Reduced density matrix
----------------------

.. automodule:: qDNA.dynamics
   :members: get_reduced_dm, get_reduced_dm_eigs
   :show-inheritance: False

Master equation solver
----------------------

.. automodule:: qDNA.dynamics
   :members: get_me_solver
   :show-inheritance: False
   :noindex:
   
   
Exciton Lifetime 
----------------

.. automodule:: qDNA.evaluation
   :members: calc_lifetime, calc_lifetime_dict
   :show-inheritance: False
   
Average charge separation
-------------------------

.. automodule:: qDNA.evaluation
   :members: calc_dipole, calc_dipole_wrapper, calc_dipole_dict
   :show-inheritance: False
   :noindex:
   
Equilibrium states
------------------

.. automodule:: qDNA.evaluation
   :members: get_therm_eq_state, get_deph_eq_state
   :show-inheritance: False
   :noindex:


DNA sequences
-------------

.. automodule:: qDNA
   :members: create_upper_strands
   :show-inheritance: False
   
   
Observables
-----------

.. automodule:: qDNA.environment
   :members: get_tb_observable, get_eh_observable, get_pop_particle, get_coh_particle
   :show-inheritance: False

Lindblad rates
--------------

.. automodule:: qDNA.environment
   :members: debye_spectral_density, ohmic_spectral_density, bose_einstein_distrib, dephasing_rate
   :show-inheritance: False
   :noindex:
   
Lindblad operators
------------------

.. automodule:: qDNA.environment
   :members: get_relax_op, get_glob_therm_op, get_glob_therm_ops, get_loc_therm_op, get_loc_therm_ops, get_loc_deph_ops, get_glob_deph_ops, get_loc_deph_p_ops, get_glob_deph_p_ops
   :show-inheritance: False
   :noindex:
   
   

   
   
   
   
   
   
   
   
   
   
   
   
   
   


