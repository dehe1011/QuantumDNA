Functions
=========

Calculate Tight-Binding Parameters
----------------------------------

.. autofunction:: qDNA.lcao.calc_orbital_interaction
.. autofunction:: qDNA.lcao.calc_orbital_energy
.. autofunction:: qDNA.lcao.calc_H_intra
.. autofunction:: qDNA.lcao.calc_H_inter

Save and Load Tight-Binding Parameters
--------------------------------------

.. autofunction:: qDNA.io.save_tb_params
.. autofunction:: qDNA.io.load_tb_params
.. autofunction:: qDNA.io.delete_tb_params

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

.. autofunction:: qDNA.model.get_tb_couplings
.. autofunction:: qDNA.model.get_tb_energies
.. autofunction:: qDNA.model.get_tb_config


Tight-Binding Hamiltonian
-------------------------

.. autofunction:: qDNA.hamiltonian.set_matrix_element
.. autofunction:: qDNA.hamiltonian.tb_ham_1P
.. autofunction:: qDNA.hamiltonian.tb_ham_2P
.. autofunction:: qDNA.hamiltonian.add_groundstate
.. autofunction:: qDNA.hamiltonian.delete_groundstate
.. autofunction:: qDNA.hamiltonian.add_interaction

Lindblad rates--------------

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


Reduced Density Matrix
----------------------

.. autofunction:: qDNA.dynamics.get_reduced_dm
.. autofunction:: qDNA.dynamics.get_reduced_dm_eigs

Observables
-----------

.. autofunction:: qDNA.environment.get_tb_observable
.. autofunction:: qDNA.environment.get_eh_observable
.. autofunction:: qDNA.environment.get_pop_observable
.. autofunction:: qDNA.environment.get_coh_observable

Equilibrium states
------------------

.. autofunction:: qDNA.evaluation.get_therm_eq_state
.. autofunction:: qDNA.evaluation.get_deph_eq_state

Unit Conversion
---------------

.. autofunction:: qDNA.utils.get_conversion
.. autofunction:: qDNA.utils.get_all_conversions
.. autofunction:: qDNA.utils.get_conversion_dict


Hamiltonian Analysis
--------------------

.. autofunction:: qDNA.hamiltonian.calc_average_pop
.. autofunction:: qDNA.hamiltonian.calc_amplitudes
.. autofunction:: qDNA.hamiltonian.calc_frequencies
.. autofunction:: qDNA.hamiltonian.get_pop_fourier
.. autofunction:: qDNA.hamiltonian.calc_ipr_hamiltonian


Density Matrix Analysis
-----------------------

.. autofunction:: qDNA.dynamics.calc_trace_distance
.. autofunction:: qDNA.dynamics.calc_purity
.. autofunction:: qDNA.dynamics.calc_coherence
.. autofunction:: qDNA.dynamics.calc_ipr_dm

Evaluation
----------

.. autofunction:: qDNA.evaluation.evaluate_pdb
