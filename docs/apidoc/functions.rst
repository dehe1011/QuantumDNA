Functions
=========

Calculate Tight-Binding Parameters
----------------------------------

.. automodule:: qDNA.lcao
   :members: calc_orbital_energy, calc_orbital_overlap, load_xyz, convert_json_to_xyz, convert_pdb_to_xyz, calc_tb_energies_monomers, calc_tb_params_dimer, calc_tb_params
   :show-inheritance: False

Save and Load Tight-Binding Parameters
--------------------------------------

.. automodule:: qDNA.hamiltonian
   :members: save_tb_params, load_tb_params, wrap_save_tb_params, wrap_load_tb_params
   :show-inheritance: False

Tight-Binding Basis
-------------------

There are three relevant bases:

- The `tight-binding basis` or one-particle basis, which is the basis in which the tight-binding Hamiltonian is defined.
- The `electron-hole basis` or two-particle basis, which is the basis in which the electron-hole states are defined.
- The `eigenstate basis`, which is the basis in which the Hamiltonian is diagonal.

The first two are refered to as `local` bases, the third as the `global` basis. The functions below define these bases and allow to change between them.

.. automodule:: qDNA.model
   :members: get_tb_basis, get_eh_basis, get_eh_distance, get_particle_eh_states, basis_change, local_to_global, global_to_local
   :show-inheritance: False

Tight-Binding Configuration
---------------------------

.. automodule:: qDNA.model
   :members: get_tb_config
   :show-inheritance: False

DNA sequences
-------------

.. automodule:: qDNA
   :members: create_upper_strands
   :show-inheritance: False

Tight-Binding Hamiltonian
-------------------------

.. automodule:: qDNA.hamiltonian
   :members: set_matrix_element, tb_ham_1P, tb_ham_2P, add_groundstate, delete_groundstate, add_interaction
   :show-inheritance: False

Lindblad rates
--------------

These functions is adapted from the quantum_HEOM GitHub repository :cite:`Abbott2020`.

.. automodule:: qDNA.environment
   :members: debye_spectral_density, ohmic_spectral_density, bose_einstein_distrib, dephasing_rate
   :show-inheritance: False

Lindblad operators
------------------

.. automodule:: qDNA.environment
   :members: get_relax_op, get_glob_therm_op, get_glob_therm_ops, get_loc_therm_op, get_loc_therm_ops, get_loc_deph_ops, get_glob_deph_ops, get_loc_deph_p_ops, get_glob_deph_p_ops
   :show-inheritance: False

Master Equation Solver
----------------------

.. automodule:: qDNA.dynamics
   :members: get_me_solver
   :show-inheritance: False

Reduced Density Matrix
----------------------

.. automodule:: qDNA.dynamics
   :members: get_reduced_dm, get_reduced_dm_eigs
   :show-inheritance: False

Exciton Observables
-------------------

.. automodule:: qDNA.evaluation
   :members: calc_lifetime, calc_lifetime_dict
   :show-inheritance: False

.. automodule:: qDNA.evaluation
   :members: calc_dipole, calc_dipole_wrapper, calc_dipole_dict
   :show-inheritance: False

Observables
-----------

.. automodule:: qDNA.environment
   :members: get_tb_observable, get_eh_observable, get_pop_particle, get_coh_particle
   :show-inheritance: False

Equilibrium states
------------------

.. automodule:: qDNA.evaluation
   :members: get_therm_eq_state, get_deph_eq_state
   :show-inheritance: False
