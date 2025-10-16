Functions
=========

Calculate Tight-Binding Parameters
----------------------------------

.. autofunction:: qDNA.lcao.calc_orbital_energy
.. autofunction:: qDNA.lcao.calc_orbital_overlap
.. autofunction:: qDNA.lcao.load_xyz
.. autofunction:: qDNA.lcao.convert_json_to_xyz
.. autofunction:: qDNA.lcao.convert_pdb_to_xyz
.. autofunction:: qDNA.lcao.calc_tb_energies_monomers
.. autofunction:: qDNA.lcao.calc_tb_params_dimer
.. autofunction:: qDNA.lcao.calc_tb_params

Save and Load Tight-Binding Parameters
--------------------------------------

.. autofunction:: qDNA.hamiltonian.save_tb_params
.. autofunction:: qDNA.hamiltonian.load_tb_params
.. autofunction:: qDNA.hamiltonian.wrap_save_tb_params
.. autofunction:: qDNA.hamiltonian.wrap_load_tb_params

Tight-Binding Basis
-------------------

There are three relevant bases:

- The `tight-binding basis` or one-particle basis, which is the basis in which the tight-binding Hamiltonian is defined.
- The `electron-hole basis` or two-particle basis, which is the basis in which the electron-hole states are defined.
- The `eigenstate basis`, which is the basis in which the Hamiltonian is diagonal.

The first two are refered to as `local` bases, the third as the `global` basis. The functions below define these bases and allow to change between them.

.. autofunction:: qDNA.model.get_tb_basis
.. autofunction:: qDNA.model.get_eh_basis
.. autofunction:: qDNA.model.get_eh_distance
.. autofunction:: qDNA.model.get_particle_eh_states
.. autofunction:: qDNA.model.basis_change
.. autofunction:: qDNA.model.local_to_global
.. autofunction:: qDNA.model.global_to_local


Tight-Binding Configuration
---------------------------

.. autofunction:: qDNA.model.get_tb_config


DNA sequences
-------------

.. autofunction:: qDNA.create_upper_strands


Tight-Binding Hamiltonian
-------------------------

.. autofunction:: qDNA.hamiltonian.set_matrix_element
.. autofunction:: qDNA.hamiltonian.tb_ham_1P
.. autofunction:: qDNA.hamiltonian.tb_ham_2P
.. autofunction:: qDNA.hamiltonian.add_groundstate
.. autofunction:: qDNA.hamiltonian.delete_groundstate
.. autofunction:: qDNA.hamiltonian.add_interaction

Lindblad rates
--------------

These functions is adapted from the quantum_HEOM GitHub repository :cite:`Abbott2020`.

.. autofunction:: qDNA.environment.debye_spectral_density
.. autofunction:: qDNA.environment.ohmic_spectral_density
.. autofunction:: qDNA.environment.bose_einstein_distrib
.. autofunction:: qDNA.environment.dephasing_rate

Lindblad operators
------------------

.. autofunction:: qDNA.environment.get_relax_op
.. autofunction:: qDNA.environment.get_glob_therm_op
.. autofunction:: qDNA.environment.get_glob_therm_ops
.. autofunction:: qDNA.environment.get_loc_therm_op
.. autofunction:: qDNA.environment.get_loc_therm_ops
.. autofunction:: qDNA.environment.get_loc_deph_ops
.. autofunction:: qDNA.environment.get_glob_deph_ops
.. autofunction:: qDNA.environment.get_loc_deph_p_ops
.. autofunction:: qDNA.environment.get_glob_deph_p_ops


Master Equation Solver
----------------------

.. autofunction:: qDNA.dynamics.get_me_solver


Reduced Density Matrix
----------------------

.. autofunction:: qDNA.dynamics.get_reduced_dm
.. autofunction:: qDNA.dynamics.get_reduced_dm_eigs


Exciton Observables
-------------------

.. autofunction:: qDNA.evaluation.calc_lifetime
.. autofunction:: qDNA.evaluation.calc_lifetime_dict


.. autofunction:: qDNA.evaluation.calc_dipole
.. autofunction:: qDNA.evaluation.calc_dipole_wrapper
.. autofunction:: qDNA.evaluation.calc_dipole_dict


Observables
-----------

.. autofunction:: qDNA.environment.get_tb_observable
.. autofunction:: qDNA.environment.get_eh_observable
.. autofunction:: qDNA.environment.get_pop_particle
.. autofunction:: qDNA.environment.get_coh_particle

Equilibrium states
------------------

.. autofunction:: qDNA.evaluation.get_therm_eq_state
.. autofunction:: qDNA.evaluation.get_deph_eq_state
