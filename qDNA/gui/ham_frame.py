# pylint: skip-file

import customtkinter as ctk

# --------------------------------------------------


class HamFrame(ctk.CTkFrame):
    def __init__(self, master, configs, **kwargs):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        super().__init__(master)
        self.pack(fill="both", expand=True)

        # Store configurations
        self.configs = configs
        self.ham_kwargs_default = self.configs["ham_kwargs_default"]

        self.tb_params_label = ctk.CTkLabel(
            self, text="TB Parameters", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.tb_params_label.grid(row=0, column=0, pady=10, padx=10)

        # Source widgets
        self.ham_source_label = ctk.CTkLabel(self, text="Source:")
        self.ham_source_label.grid(row=1, column=0, padx=10, pady=10)
        self.ham_source_combo = ctk.CTkComboBox(self, values=self.configs["SOURCES"])
        self.ham_source_combo.set(self.ham_kwargs_default["source"])
        self.ham_source_combo.grid(row=1, column=1, padx=10, pady=10)

        # Description widgets
        self.ham_description_label = ctk.CTkLabel(self, text="Description:")
        self.ham_description_label.grid(row=2, column=0, padx=10, pady=10)
        self.ham_description_combo = ctk.CTkComboBox(
            self, values=self.configs["DESCRIPTIONS"]
        )
        self.ham_description_combo.set(self.ham_kwargs_default["description"])
        self.ham_description_combo.grid(row=2, column=1, padx=10, pady=10)

        # Particles widgets
        self.ham_particles_label = ctk.CTkLabel(self, text="Particles:")
        self.ham_particles_label.grid(row=3, column=0, padx=10, pady=10)
        self.particles = self.configs["PARTICLES"]
        self.selected_particles = {
            particle: ctk.BooleanVar(value=True) for particle in self.particles
        }
        for idx, particle in enumerate(self.particles):
            checkbox = ctk.CTkCheckBox(
                self, text=particle, variable=self.selected_particles[particle]
            )
            checkbox.grid(row=3 + idx, column=1, padx=10, pady=5, sticky="w")

        # Unit widgets
        self.ham_unit_label = ctk.CTkLabel(self, text="Unit:")
        self.ham_unit_label.grid(row=6, column=0, padx=10, pady=10)
        self.ham_unit_combo = ctk.CTkComboBox(self, values=self.configs["UNITS"])
        self.ham_unit_combo.set(self.ham_kwargs_default["unit"])
        self.ham_unit_combo.grid(row=6, column=1, padx=10, pady=10)

        self.tb_params_label = ctk.CTkLabel(
            self,
            text="Excitonic \nInteractions",
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self.tb_params_label.grid(row=7, column=0, pady=10, padx=10)

        # Coulomb interaction widgets
        self.ham_coulomb_param_label = ctk.CTkLabel(self, text="Coulomb Interaction:")
        self.ham_coulomb_param_label.grid(row=8, column=0, padx=10, pady=10)
        self.ham_coulomb_param_entry = ctk.CTkEntry(self)
        self.ham_coulomb_param_entry.insert(
            0, str(self.ham_kwargs_default["coulomb_param"])
        )
        self.ham_coulomb_param_entry.grid(row=8, column=1, padx=10, pady=10)

        # Exchange interaction widgets
        self.ham_exchange_param_label = ctk.CTkLabel(self, text="Exchange Interaction:")
        self.ham_exchange_param_label.grid(row=9, column=0, padx=10, pady=10)
        self.ham_exchange_param_entry = ctk.CTkEntry(self)
        self.ham_exchange_param_entry.insert(
            0, str(self.ham_kwargs_default["exchange_param"])
        )
        self.ham_exchange_param_entry.grid(row=9, column=1, padx=10, pady=10)

        # NN Cutoff checkbox
        self.ham_nn_cutoff_var = ctk.BooleanVar(
            value=self.ham_kwargs_default["nn_cutoff"]
        )
        self.ham_nn_cutoff_check = ctk.CTkCheckBox(
            self, text="Nearest-Neighbor \nCutoff", variable=self.ham_nn_cutoff_var
        )
        self.ham_nn_cutoff_check.grid(row=10, column=1, padx=10, pady=10, sticky="w")

        self.relaxation_label = ctk.CTkLabel(
            self, text="DNA Relaxation", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.relaxation_label.grid(row=11, column=0, pady=10, padx=10)

        # Relaxation checkbox
        self.ham_relaxation_var = ctk.BooleanVar(
            value=self.ham_kwargs_default["relaxation"]
        )
        self.ham_relaxation_check = ctk.CTkCheckBox(
            self, text="Groundstate", variable=self.ham_relaxation_var
        )
        self.ham_relaxation_check.grid(row=12, column=1, padx=10, pady=10, sticky="w")

    def get_ham_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        ham_kwargs = {
            "source": self.ham_source_combo.get(),
            "description": self.ham_description_combo.get(),
            "particles": [
                particle
                for particle, var in self.selected_particles.items()
                if var.get()
            ],
            "unit": self.ham_unit_combo.get(),
            "coulomb_param": float(self.ham_coulomb_param_entry.get()),
            "exchange_param": float(self.ham_exchange_param_entry.get()),
            "relaxation": self.ham_relaxation_var.get(),
            "nn_cutoff": self.ham_nn_cutoff_var.get(),
        }
        return ham_kwargs
