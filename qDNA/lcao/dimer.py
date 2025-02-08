import os
import json

import numpy as np

from .slater_koster import calc_orbital_overlap

# ----------------------------------------------


class Dimer:
    """
    A class to represent a dimer consisting of two molecules.

    Parameters
    ----------
    molecule1 : Molecule
        The first molecule in the dimer.
    molecule2 : Molecule
        The second molecule in the dimer.
    identifier : str, optional
        An identifier for the dimer. If not provided, it is generated from the identifiers of the two molecules.

    Attributes
    ----------
    molecule1 : Molecule
        The first molecule in the dimer.
    molecule2 : Molecule
        The second molecule in the dimer.
    identifier : str, optional
        An identifier for the dimer. If not provided, it is generated from the identifiers of the two molecules.
    H_int : numpy.ndarray
        The interaction Hamiltonian matrix between the two molecules.
    t_HOMO : float
        The hopping parameter for the highest occupied molecular orbital (HOMO).
    t_HOMO_1 : float
        The hopping parameter for the first HOMO-1.
    t_HOMO_2 : float
        The hopping parameter for the second HOMO-2.
    t_LUMO : float
        The hopping parameter for the lowest unoccupied molecular orbital (LUMO).
    t_LUMO_1 : float
        The hopping parameter for the first LUMO-1.
    t_LUMO_2 : float
        The hopping parameter for the second LUMO-2.

    Methods
    -------
    save_results(directory="results"):
        Saves the results of the calculation to a JSON file.
    """

    def __init__(self, molecule1, molecule2, identifier=None):
        self.molecule1 = molecule1
        self.molecule2 = molecule2

        # Identifier
        self.identifier = identifier
        if self.identifier is None:
            self.identifier = molecule1.identifier + molecule2.identifier

        # Interaction Hamiltonian (not the full LCAO Hamiltonian)
        self.H_int = self._calc_H_int()

        # Hopping parameters for electron and hole
        self.t_HOMO = self.molecule1.HOMO @ self.H_int @ self.molecule2.HOMO
        self.t_HOMO_1 = self.molecule1.HOMO_1 @ self.H_int @ self.molecule2.HOMO_1
        self.t_HOMO_2 = self.molecule1.HOMO_2 @ self.H_int @ self.molecule2.HOMO_2
        self.t_LUMO = self.molecule1.LUMO @ self.H_int @ self.molecule2.LUMO
        self.t_LUMO_1 = self.molecule1.LUMO_1 @ self.H_int @ self.molecule2.LUMO_1
        self.t_LUMO_2 = self.molecule1.LUMO_2 @ self.H_int @ self.molecule2.LUMO_2

    def __repr__(self):
        return f'BasePair({self.molecule1}, {self.molecule2}, identifier = "{self.identifier}")'

    def save_results(self, directory="results"):
        """
        Save the results of the calculation to a JSON file.

        Parameters
        ----------
        directory : str, optional
            The directory where the results file will be saved (default is "results").

        Notes
        -----
        .. note::

            The results include the HOMO and LUMO energies of two molecules and their coupling terms.
            The results are saved in a JSON file named after the identifier of the object.
        """

        dictionary = {
            "E1_HOMO": round(float(self.molecule1.E_HOMO), 4),
            "E2_HOMO": round(float(self.molecule2.E_HOMO), 4),
            "t_HOMO": round(float(self.t_HOMO), 4),
            "E1_LUMO": round(float(self.molecule1.E_LUMO), 4),
            "E2_LUMO": round(float(self.molecule2.E_LUMO), 4),
            "t_LUMO": round(float(self.t_LUMO), 4),
        }

        filename = self.identifier + ".json"
        filepath = os.path.join(directory, filename)

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(dictionary, f, indent=2)
        print("Results saved at" + filepath)

    def _calc_H_int(self):
        H_int = np.zeros((self.molecule1.num_orbitals, self.molecule2.num_orbitals))
        for i, orbital1 in enumerate(self.molecule1.orbitals):
            for j, orbital2 in enumerate(self.molecule2.orbitals):
                value = calc_orbital_overlap(orbital1, orbital2, "interbase")
                H_int[i, j] = value
        return H_int
