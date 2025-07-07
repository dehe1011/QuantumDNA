# pylint: skip-file

import customtkinter as ctk
from ..io import OPTIONS, DEFAULTS

# ----------------------------------------------------------------------


class TBHamFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        super().__init__(master)
        self.pack(fill="both", expand=True)

        # Store configurations
        self.options = OPTIONS
        self.defaults = DEFAULTS["tb_ham_kwargs_default"]

        # --------------------------------------------------------------

        # TB parameters label
        self.tb_params_label = ctk.CTkLabel(
            self, text="TB Parameters", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.tb_params_label.grid(row=0, column=0, pady=10, padx=10)

        # Source label and combo box
        self.source_label = ctk.CTkLabel(self, text="Source:")
        self.source_label.grid(row=1, column=0, padx=10, pady=10)
        self.source_combo = ctk.CTkComboBox(self, values=self.options["sources"])
        self.source_combo.set(self.defaults["source"])
        self.source_combo.grid(row=1, column=1, padx=10, pady=10)

        # Description label and combo box
        self.description_label = ctk.CTkLabel(self, text="Description:")
        self.description_label.grid(row=2, column=0, padx=10, pady=10)
        self.description_combo = ctk.CTkComboBox(
            self, values=self.options["descriptions"]
        )
        self.description_combo.set(self.defaults["description"])
        self.description_combo.grid(row=2, column=1, padx=10, pady=10)

        # Particles label and checkboxes
        self.particles_label = ctk.CTkLabel(self, text="Particles:")
        self.particles_label.grid(row=3, column=0, padx=10, pady=10)
        self.particles = self.options["particles"]
        self.selected_particles = {
            particle: ctk.BooleanVar(value=True) for particle in self.particles
        }
        for idx, particle in enumerate(self.particles):
            checkbox = ctk.CTkCheckBox(
                self, text=particle, variable=self.selected_particles[particle]
            )
            checkbox.grid(row=3 + idx, column=1, padx=10, pady=5, sticky="w")

        # Unit label and combo box
        self.unit_label = ctk.CTkLabel(self, text="Unit:")
        self.unit_label.grid(row=6, column=0, padx=10, pady=10)
        self.unit_combo = ctk.CTkComboBox(self, values=self.options["units"])
        self.unit_combo.set(self.defaults["unit"])
        self.unit_combo.grid(row=6, column=1, padx=10, pady=10)

        # --------------------------------------------------------------

        # Exciton parameters label
        self.tb_params_label = ctk.CTkLabel(
            self,
            text="Exciton \nParameters",
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self.tb_params_label.grid(row=7, column=0, pady=10, padx=10)

        # Coulomb interaction label and entry
        self.coulomb_param_label = ctk.CTkLabel(self, text="Coulomb Interaction:")
        self.coulomb_param_label.grid(row=8, column=0, padx=10, pady=10)
        self.coulomb_param_entry = ctk.CTkEntry(self)
        self.coulomb_param_entry.insert(0, str(self.defaults["coulomb_param"]))
        self.coulomb_param_entry.grid(row=8, column=1, padx=10, pady=10)

        # Exchange interaction label and entry
        self.exchange_param_label = ctk.CTkLabel(self, text="Exchange Interaction:")
        self.exchange_param_label.grid(row=9, column=0, padx=10, pady=10)
        self.exchange_param_entry = ctk.CTkEntry(self)
        self.exchange_param_entry.insert(0, str(self.defaults["exchange_param"]))
        self.exchange_param_entry.grid(row=9, column=1, padx=10, pady=10)

        # NN Cutoff checkbox
        self.nn_cutoff_var = ctk.BooleanVar(value=self.defaults["nn_cutoff"])
        self.nn_cutoff_check = ctk.CTkCheckBox(
            self, text="Nearest-Neighbor \nCutoff", variable=self.nn_cutoff_var
        )
        self.nn_cutoff_check.grid(row=10, column=1, padx=10, pady=10, sticky="w")

        # relaxation label
        self.relaxation_label = ctk.CTkLabel(
            self, text="DNA Relaxation", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.relaxation_label.grid(row=11, column=0, pady=10, padx=10)

        # relaxation checkbox
        self.relaxation_var = ctk.BooleanVar(value=self.defaults["relaxation"])
        self.relaxation_check = ctk.CTkCheckBox(
            self, text="Groundstate", variable=self.relaxation_var
        )
        self.relaxation_check.grid(row=12, column=1, padx=10, pady=10, sticky="w")

        # --------------------------------------------------------------

    def get_tb_ham_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        tb_ham_kwargs = {
            "source": self.source_combo.get(),
            "description": self.description_combo.get(),
            "particles": [
                particle
                for particle, var in self.selected_particles.items()
                if var.get()
            ],
            "unit": self.unit_combo.get(),
            "coulomb_param": float(self.coulomb_param_entry.get()),
            "exchange_param": float(self.exchange_param_entry.get()),
            "relaxation": self.relaxation_var.get(),
            "nn_cutoff": self.nn_cutoff_var.get(),
        }
        return tb_ham_kwargs


# ----------------------------------------------------------------------
