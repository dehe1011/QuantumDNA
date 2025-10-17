# pylint: skip-file

import customtkinter as ctk
from ..io import OPTIONS, DEFAULTS

# ----------------------------------------------------------------------


class LindDissFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)

        # Store configurations
        self.options = OPTIONS
        self.defaults = DEFAULTS["lind_diss_kwargs_default"]

        # --------------------------------------------------------------

        # dephasing label
        self.dephasing_label = ctk.CTkLabel(
            self, text="Dephasing", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.dephasing_label.grid(row=0, column=0, pady=10, padx=10)

        # local dephasing rate label and entry
        self.loc_deph_rate_label = ctk.CTkLabel(self, text="Local Dephasing Rate:")
        self.loc_deph_rate_label.grid(row=1, column=0, padx=10, pady=10)

        self.loc_deph_rate_entry = ctk.CTkEntry(self)
        self.loc_deph_rate_entry.insert(0, str(self.defaults["loc_deph_rate"]))
        self.loc_deph_rate_entry.grid(row=1, column=1, padx=10, pady=10)

        # global dephasing rate label and entry
        self.glob_deph_rate_label = ctk.CTkLabel(self, text="Global Dephasing Rate:")
        self.glob_deph_rate_label.grid(row=2, column=0, padx=10, pady=10)

        self.glob_deph_rate_entry = ctk.CTkEntry(self)
        self.glob_deph_rate_entry.insert(0, str(self.defaults["glob_deph_rate"]))
        self.glob_deph_rate_entry.grid(row=2, column=1, padx=10, pady=10)

        # --------------------------------------------------------------

        # relaxation label
        self.relaxation_label = ctk.CTkLabel(
            self, text="Relaxation", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.relaxation_label.grid(row=3, column=0, pady=10, padx=10)

        # relaxation rate label and entry
        self.relax_rate_label = ctk.CTkLabel(self, text="Relaxation Rate:")
        self.relax_rate_label.grid(row=4, column=0, padx=10, pady=10)

        self.relax_rate_entry = ctk.CTkEntry(self)
        # self.relax_rate_entry.insert(
        #     0, str(self.defaults["relax_rate"])
        # )
        self.relax_rate_entry.insert(0, 3.0)
        self.relax_rate_entry.grid(row=4, column=1, padx=10, pady=10)

        # -----------------------------------

        # thermalization label
        self.therm_label = ctk.CTkLabel(
            self, text="Thermalization", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.therm_label.grid(row=5, column=0, pady=10, padx=10)

        # local thermalization checkbox
        self.loc_therm_var = ctk.BooleanVar(value=self.defaults["loc_therm"])
        self.loc_therm_check = ctk.CTkCheckBox(self, text="Local", variable=self.loc_therm_var)
        self.loc_therm_check.grid(row=6, column=0, padx=10, pady=10)

        # global thermalization checkbox
        self.glob_therm_var = ctk.BooleanVar(value=self.defaults["glob_therm"])
        self.glob_therm_check = ctk.CTkCheckBox(self, text="Global", variable=self.glob_therm_var)
        self.glob_therm_check.grid(row=6, column=1, padx=10, pady=10)

        # --------------------------------------------------------------

        # dephasing rate label and entry
        self.deph_rate_label = ctk.CTkLabel(self, text="Dephasing Rate:")
        self.deph_rate_label.grid(row=7, column=0, padx=10, pady=10)

        self.deph_rate_entry = ctk.CTkEntry(self)
        self.deph_rate_entry.insert(0, str(self.defaults["deph_rate"]))
        self.deph_rate_entry.grid(row=7, column=1, padx=10, pady=10)

        # cutoff frequency label and entry
        self.cutoff_freq_label = ctk.CTkLabel(self, text="Cutoff Frequency:")
        self.cutoff_freq_label.grid(row=8, column=0, padx=10, pady=10)

        self.cutoff_freq_entry = ctk.CTkEntry(self)
        self.cutoff_freq_entry.insert(0, str(self.defaults["cutoff_freq"]))
        self.cutoff_freq_entry.grid(row=8, column=1, padx=10, pady=10)

        # reorganization energy label and entry
        self.reorg_energy_label = ctk.CTkLabel(self, text="Reorganization Energy:")
        self.reorg_energy_label.grid(row=9, column=0, padx=10, pady=10)

        self.reorg_energy_entry = ctk.CTkEntry(self)
        self.reorg_energy_entry.insert(0, str(self.defaults["reorg_energy"]))
        self.reorg_energy_entry.grid(row=9, column=1, padx=10, pady=10)

        # temperature label and entry
        self.temperature_label = ctk.CTkLabel(self, text="Temperature (K):")
        self.temperature_label.grid(row=10, column=0, padx=10, pady=10)

        self.temperature_entry = ctk.CTkEntry(self)
        self.temperature_entry.insert(0, str(self.defaults["temperature"]))
        self.temperature_entry.grid(row=10, column=1, padx=10, pady=10)

        # spectral density label and combo box
        self.spectral_density_label = ctk.CTkLabel(self, text="Spectral Density:")
        self.spectral_density_label.grid(row=11, column=0, padx=10, pady=10)

        self.spectral_density_combo = ctk.CTkComboBox(
            self, values=self.options["spectral_densities"]
        )
        self.spectral_density_combo.set(self.defaults["spectral_density"])
        self.spectral_density_combo.grid(row=11, column=1, padx=10, pady=10)

        # exponent label and entry
        self.exponent_label = ctk.CTkLabel(self, text="Exponent:")
        self.exponent_label.grid(row=12, column=0, padx=10, pady=10)

        self.exponent_entry = ctk.CTkEntry(self)
        self.exponent_entry.insert(0, str(self.defaults["exponent"]))
        self.exponent_entry.grid(row=12, column=1, padx=10, pady=10)

        # --------------------------------------------------------------

    def set_uniform_relax(self):
        if self.uniform_relaxation_var.get():
            self.relax_rate_entry.configure(state="normal")
            self.relax_rates_entry.configure(state="disabled")
        else:
            self.relax_rate_entry.configure(state="disabled")
            self.relax_rates_entry.configure(state="normal")

    def get_lind_diss_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        lind_diss_kwargs = {
            "loc_deph_rate": float(self.loc_deph_rate_entry.get()),
            "glob_deph_rate": float(self.glob_deph_rate_entry.get()),
            "uniform_relaxation": True,
            "relax_rate": float(self.relax_rate_entry.get()),
            "loc_therm": self.loc_therm_var.get(),
            "glob_therm": self.glob_therm_var.get(),
            "deph_rate": float(self.deph_rate_entry.get()),
            "cutoff_freq": float(self.cutoff_freq_entry.get()),
            "reorg_energy": float(self.reorg_energy_entry.get()),
            "temperature": float(self.temperature_entry.get()),
            "spectral_density": self.spectral_density_combo.get(),
            "exponent": float(self.exponent_entry.get()),
        }
        return lind_diss_kwargs


# ----------------------------------------------------------------------
