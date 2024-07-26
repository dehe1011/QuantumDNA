import customtkinter as ctk

from qDNA import wrap_save_tb_params
from qDNA.tools import get_config


class HamFrame(ctk.CTkFrame):
    def __init__(self, master, configs, **kwargs):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.configs = configs
        self.ham_kwargs_default = self.configs["ham_kwargs_default"]

        # widgets
        self.ham_source_label = ctk.CTkLabel(self, text="Source:")
        self.ham_source_label.grid(row=0, column=0, padx=10, pady=10)

        self.ham_source_combo = ctk.CTkComboBox(self, values=self.configs["SOURCES"])
        self.ham_source_combo.set(self.ham_kwargs_default["source"])
        self.ham_source_combo.grid(row=0, column=1, padx=10, pady=10)

        self.ham_description_label = ctk.CTkLabel(self, text="Description:")
        self.ham_description_label.grid(row=1, column=0, padx=10, pady=10)

        self.ham_description_combo = ctk.CTkComboBox(
            self, values=self.configs["DESCRIPTIONS"]
        )
        self.ham_description_combo.set(self.ham_kwargs_default["description"])
        self.ham_description_combo.grid(row=1, column=1, padx=10, pady=10)

        self.ham_particles_label = ctk.CTkLabel(self, text="Particles:")
        self.ham_particles_label.grid(row=2, column=0, padx=10, pady=10)

        self.particles = self.configs["PARTICLES"]
        self.selected_particles = {
            particle: ctk.BooleanVar(value=True) for particle in self.particles
        }
        for idx, particle in enumerate(self.particles):
            checkbox = ctk.CTkCheckBox(
                self, text=particle, variable=self.selected_particles[particle]
            )
            checkbox.grid(row=3 + idx, column=0, padx=20, pady=5, columnspan=2)

        self.ham_unit_label = ctk.CTkLabel(self, text="Unit:")
        self.ham_unit_label.grid(row=6, column=0, padx=10, pady=10)

        self.ham_unit_combo = ctk.CTkComboBox(self, values=self.configs["UNITS"])
        self.ham_unit_combo.set(self.ham_kwargs_default["unit"])
        self.ham_unit_combo.grid(row=6, column=1, padx=10, pady=10)

        self.ham_interaction_param_label = ctk.CTkLabel(self, text="Interaction Param:")
        self.ham_interaction_param_label.grid(row=7, column=0, padx=10, pady=10)

        self.ham_interaction_param_entry = ctk.CTkEntry(self)
        self.ham_interaction_param_entry.insert(
            0, str(self.ham_kwargs_default["interaction_param"])
        )
        self.ham_interaction_param_entry.grid(row=7, column=1, padx=10, pady=10)

        self.ham_relaxation_var = ctk.BooleanVar(
            value=self.ham_kwargs_default["relaxation"]
        )
        self.ham_relaxation_check = ctk.CTkCheckBox(
            self, text="Relaxation", variable=self.ham_relaxation_var
        )
        self.ham_relaxation_check.grid(row=12, column=0, padx=10, pady=10, columnspan=2)

        self.ham_nn_cutoff_var = ctk.BooleanVar(
            value=self.ham_kwargs_default["nn_cutoff"]
        )
        self.ham_nn_cutoff_check = ctk.CTkCheckBox(
            self, text="NN Cutoff", variable=self.ham_nn_cutoff_var
        )
        self.ham_nn_cutoff_check.grid(row=13, column=0, padx=10, pady=10, columnspan=2)

    def get_ham_kwargs(self):
        """
        Returns the values of widgets with get() method in dictionary format.
        """

        ham_kwargs = {
            "source": self.ham_source_combo.get(),
            "description": self.ham_description_combo.get(),
            "particles": [
                particle
                for particle, var in self.selected_particles.items()
                if var.get()
            ],
            "unit": self.ham_unit_combo.get(),
            "interaction_param": float(self.ham_interaction_param_entry.get()),
            "relaxation": self.ham_relaxation_var.get(),
            "nn_cutoff": self.ham_nn_cutoff_var.get(),
        }
        return ham_kwargs
