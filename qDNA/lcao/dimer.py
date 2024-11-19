import os
import json

import numpy as np

from .slater_koster import calc_orbital_overlap

# ----------------------------------------------


class Dimer:
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
        dict = {
            "E1_HOMO": round(float(self.molecule1.E_HOMO), 4),
            "E2_HOMO": round(float(self.molecule2.E_HOMO), 4),
            "t_HOMO": round(float(self.t_HOMO), 4),
            "E1_LUMO": round(float(self.molecule1.E_LUMO), 4),
            "E2_LUMO": round(float(self.molecule2.E_LUMO), 4),
            "t_LUMO": round(float(self.t_LUMO), 4),
        }

        filename = self.identifier + ".json"
        filepath = os.path.join(directory, filename)

        with open(filepath, "w") as f:
            json.dump(dict, f, indent=2)
        print("Results saved at" + filepath)

    def _calc_H_int(self):
        H_int = np.zeros((self.molecule1.num_orbitals, self.molecule2.num_orbitals))
        for i, orbital1 in enumerate(self.molecule1.orbitals):
            for j, orbital2 in enumerate(self.molecule2.orbitals):
                value = calc_orbital_overlap(orbital1, orbital2, "interbase")
                H_int[i, j] = value
        return H_int
