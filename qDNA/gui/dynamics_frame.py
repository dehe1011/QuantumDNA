import customtkinter as ctk

from qDNA.tools import get_config


class DynamicsFrame(ctk.CTkFrame):
    def __init__(self, master, configs, **kwargs):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.configs = configs
        self.me_kwargs_default = self.configs["me_kwargs_default"]

        # widgets
        self.me_t_steps_label = ctk.CTkLabel(self, text="Time Steps:")
        self.me_t_steps_label.grid(row=0, column=0, padx=10, pady=10)

        self.me_t_steps_entry = ctk.CTkEntry(self)
        self.me_t_steps_entry.insert(0, str(self.me_kwargs_default["t_steps"]))
        self.me_t_steps_entry.grid(row=0, column=1, padx=10, pady=10)

        self.me_t_end_label = ctk.CTkLabel(self, text="Time End:")
        self.me_t_end_label.grid(row=1, column=0, padx=10, pady=10)

        self.me_t_end_entry = ctk.CTkEntry(self)
        self.me_t_end_entry.insert(0, str(self.me_kwargs_default["t_end"]))
        self.me_t_end_entry.grid(row=1, column=1, padx=10, pady=10)

        self.me_t_unit_label = ctk.CTkLabel(self, text="Time Unit:")
        self.me_t_unit_label.grid(row=2, column=0, padx=10, pady=10)

        self.me_t_unit_combo = ctk.CTkComboBox(self, values=self.configs["T_UNITS"])
        self.me_t_unit_combo.set(self.me_kwargs_default["t_unit"])
        self.me_t_unit_combo.grid(row=2, column=1, padx=10, pady=10)

        self.me_init_e_state_label = ctk.CTkLabel(self, text="Initial Electron State:")
        self.me_init_e_state_label.grid(row=3, column=0, padx=10, pady=10)

        self.me_init_e_state_combo = ctk.CTkComboBox(
            self, values=self.me_kwargs_default["init_e_state"]
        )
        self.me_init_e_state_combo.grid(row=3, column=1, padx=10, pady=10)
        self.me_init_e_state_combo.set(self.me_kwargs_default["init_e_state"])

        self.me_init_h_state_label = ctk.CTkLabel(self, text="Initial Hole State:")
        self.me_init_h_state_label.grid(row=4, column=0, padx=10, pady=10)

        self.me_init_h_state_combo = ctk.CTkComboBox(
            self, values=self.me_kwargs_default["init_h_state"]
        )
        self.me_init_h_state_combo.grid(row=4, column=1, padx=10, pady=10)
        self.me_init_h_state_combo.set(self.me_kwargs_default["init_h_state"])

    def get_me_kwargs(self):
        """
        Returns the values of widgets with get() method in dictionary format.
        """

        me_kwargs = {
            "t_steps": float(self.me_t_steps_entry.get()),
            "t_end": float(self.me_t_end_entry.get()),
            "t_unit": self.me_t_unit_combo.get(),
            "init_e_state": self.me_init_e_state_combo.get(),
            "init_h_state": self.me_init_h_state_combo.get(),
        }
        return me_kwargs
