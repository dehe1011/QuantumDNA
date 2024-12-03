# pylint: skip-file

import customtkinter as ctk

# --------------------------------------------------


class DissFrame(ctk.CTkFrame):
    def __init__(self, master, configs, **kwargs):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.configs = configs
        self.diss_kwargs_default = self.configs["diss_kwargs_default"]

        # widgets
        self.dephasing_label = ctk.CTkLabel(
            self, text="Dephasing", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.dephasing_label.grid(row=0, column=0, pady=10, padx=10)

        self.diss_loc_deph_rate_label = ctk.CTkLabel(self, text="Local Dephasing Rate:")
        self.diss_loc_deph_rate_label.grid(row=1, column=0, padx=10, pady=10)

        self.diss_loc_deph_rate_entry = ctk.CTkEntry(self)
        self.diss_loc_deph_rate_entry.insert(
            0, str(self.diss_kwargs_default["loc_deph_rate"])
        )
        self.diss_loc_deph_rate_entry.grid(row=1, column=1, padx=10, pady=10)

        self.diss_glob_deph_rate_label = ctk.CTkLabel(
            self, text="Global Dephasing Rate:"
        )
        self.diss_glob_deph_rate_label.grid(row=2, column=0, padx=10, pady=10)

        self.diss_glob_deph_rate_entry = ctk.CTkEntry(self)
        self.diss_glob_deph_rate_entry.insert(
            0, str(self.diss_kwargs_default["glob_deph_rate"])
        )
        self.diss_glob_deph_rate_entry.grid(row=2, column=1, padx=10, pady=10)

        # -------------------------------

        self.relaxation_label = ctk.CTkLabel(
            self, text="Relaxation", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.relaxation_label.grid(row=3, column=0, pady=10, padx=10)

        self.diss_relax_rate_label = ctk.CTkLabel(self, text="Relaxation Rate:")
        self.diss_relax_rate_label.grid(row=4, column=0, padx=10, pady=10)

        self.diss_relax_rate_entry = ctk.CTkEntry(self)
        self.diss_relax_rate_entry.insert(
            0, str(self.diss_kwargs_default["relax_rate"])
        )
        self.diss_relax_rate_entry.grid(row=4, column=1, padx=10, pady=10)

        # -----------------------------------

        self.thermalizing_label = ctk.CTkLabel(
            self, text="Thermalizing", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.thermalizing_label.grid(row=5, column=0, pady=10, padx=10)

        self.diss_loc_therm_var = ctk.BooleanVar(
            value=self.diss_kwargs_default["loc_therm"]
        )
        self.diss_loc_therm_check = ctk.CTkCheckBox(
            self, text="Local", variable=self.diss_loc_therm_var
        )
        self.diss_loc_therm_check.grid(row=6, column=0, padx=10, pady=10)

        self.diss_glob_therm_var = ctk.BooleanVar(
            value=self.diss_kwargs_default["glob_therm"]
        )
        self.diss_glob_therm_check = ctk.CTkCheckBox(
            self, text="Global", variable=self.diss_glob_therm_var
        )
        self.diss_glob_therm_check.grid(row=6, column=1, padx=10, pady=10)

        self.diss_deph_rate_label = ctk.CTkLabel(self, text="Dephasing Rate:")
        self.diss_deph_rate_label.grid(row=7, column=0, padx=10, pady=10)

        self.diss_deph_rate_entry = ctk.CTkEntry(self)
        self.diss_deph_rate_entry.insert(0, str(self.diss_kwargs_default["deph_rate"]))
        self.diss_deph_rate_entry.grid(row=7, column=1, padx=10, pady=10)

        self.diss_cutoff_freq_label = ctk.CTkLabel(self, text="Cutoff Frequency:")
        self.diss_cutoff_freq_label.grid(row=8, column=0, padx=10, pady=10)

        self.diss_cutoff_freq_entry = ctk.CTkEntry(self)
        self.diss_cutoff_freq_entry.insert(
            0, str(self.diss_kwargs_default["cutoff_freq"])
        )
        self.diss_cutoff_freq_entry.grid(row=8, column=1, padx=10, pady=10)

        self.diss_reorg_energy_label = ctk.CTkLabel(self, text="Reorganization Energy:")
        self.diss_reorg_energy_label.grid(row=9, column=0, padx=10, pady=10)

        self.diss_reorg_energy_entry = ctk.CTkEntry(self)
        self.diss_reorg_energy_entry.insert(
            0, str(self.diss_kwargs_default["reorg_energy"])
        )
        self.diss_reorg_energy_entry.grid(row=9, column=1, padx=10, pady=10)

        self.diss_temperature_label = ctk.CTkLabel(self, text="Temperature (K):")
        self.diss_temperature_label.grid(row=10, column=0, padx=10, pady=10)

        self.diss_temperature_entry = ctk.CTkEntry(self)
        self.diss_temperature_entry.insert(
            0, str(self.diss_kwargs_default["temperature"])
        )
        self.diss_temperature_entry.grid(row=10, column=1, padx=10, pady=10)

        self.diss_spectral_density_label = ctk.CTkLabel(self, text="Spectral Density:")
        self.diss_spectral_density_label.grid(row=11, column=0, padx=10, pady=10)

        self.diss_spectral_density_combo = ctk.CTkComboBox(
            self, values=self.configs["SPECTRAL_DENSITIES"]
        )
        self.diss_spectral_density_combo.set(
            self.diss_kwargs_default["spectral_density"]
        )
        self.diss_spectral_density_combo.grid(row=11, column=1, padx=10, pady=10)

        self.diss_exponent_label = ctk.CTkLabel(self, text="Exponent:")
        self.diss_exponent_label.grid(row=12, column=0, padx=10, pady=10)

        self.diss_exponent_entry = ctk.CTkEntry(self)
        self.diss_exponent_entry.insert(0, str(self.diss_kwargs_default["exponent"]))
        self.diss_exponent_entry.grid(row=12, column=1, padx=10, pady=10)

    def set_uniform_relax(self):
        if self.diss_uniform_relaxation_var.get():
            self.diss_relax_rate_entry.configure(state="normal")
            self.diss_relax_rates_entry.configure(state="disabled")
        else:
            self.diss_relax_rate_entry.configure(state="disabled")
            self.diss_relax_rates_entry.configure(state="normal")

    def get_diss_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        diss_kwargs = {
            "loc_deph_rate": float(self.diss_loc_deph_rate_entry.get()),
            "glob_deph_rate": float(self.diss_glob_deph_rate_entry.get()),
            "uniform_relaxation": True,
            "relax_rate": float(self.diss_relax_rate_entry.get()),
            "loc_therm": self.diss_loc_therm_var.get(),
            "glob_therm": self.diss_glob_therm_var.get(),
            "deph_rate": float(self.diss_deph_rate_entry.get()),
            "cutoff_freq": float(self.diss_cutoff_freq_entry.get()),
            "reorg_energy": float(self.diss_reorg_energy_entry.get()),
            "temperature": float(self.diss_temperature_entry.get()),
            "spectral_density": self.diss_spectral_density_combo.get(),
            "exponent": float(self.diss_exponent_entry.get()),
        }
        return diss_kwargs
