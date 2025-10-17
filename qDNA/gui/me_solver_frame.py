# pylint: skip-file

import customtkinter as ctk
from ..io import OPTIONS, DEFAULTS

# ----------------------------------------------------------------------


class MeSolverFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the options tab. This means that the options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        tb_basis = ["(0, 0)"]

        # Store configurations
        self.options = OPTIONS
        self.defaults = DEFAULTS["me_solver_kwargs_default"]

        # --------------------------------------------------------------

        # simulation time label
        self.time_label = ctk.CTkLabel(
            self, text="Simulation Time", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.time_label.grid(row=0, column=0, pady=10, padx=10)

        # t steps label and entry
        self.t_steps_label = ctk.CTkLabel(self, text="Time Steps:")
        self.t_steps_label.grid(row=1, column=0, padx=10, pady=10)

        self.t_steps_entry = ctk.CTkEntry(self)
        self.t_steps_entry.insert(0, str(self.defaults["t_steps"]))
        self.t_steps_entry.grid(row=1, column=1, padx=10, pady=10)

        # t end label and entry
        self.t_end_label = ctk.CTkLabel(self, text="Time End:")
        self.t_end_label.grid(row=2, column=0, padx=10, pady=10)

        self.t_end_entry = ctk.CTkEntry(self)
        self.t_end_entry.insert(0, str(self.defaults["t_end"]))
        self.t_end_entry.grid(row=2, column=1, padx=10, pady=10)

        # t unit label and combo box
        self.t_unit_label = ctk.CTkLabel(self, text="Time Unit:")
        self.t_unit_label.grid(row=3, column=0, padx=10, pady=10)

        self.t_unit_combo = ctk.CTkComboBox(self, values=self.options["t_units"])
        self.t_unit_combo.set(self.defaults["t_unit"])
        self.t_unit_combo.grid(row=3, column=1, padx=10, pady=10)

        # --------------------------------------------------------------

        # initial state label
        self.initial_state_label = ctk.CTkLabel(
            self, text="Initial State", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.initial_state_label.grid(row=4, column=0, pady=10, padx=10)

        # init e state label and combo box
        self.init_e_state_label = ctk.CTkLabel(self, text="Initial Electron State:")
        self.init_e_state_label.grid(row=5, column=0, padx=10, pady=10)

        self.init_e_state_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.init_e_state_combo.grid(row=5, column=1, padx=10, pady=10)
        self.init_e_state_combo.set(self.defaults["init_e_states"][0])

        # init h state label and combo box
        self.init_h_state_label = ctk.CTkLabel(self, text="Initial Hole State:")
        self.init_h_state_label.grid(row=6, column=0, padx=10, pady=10)

        self.init_h_state_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.init_h_state_combo.grid(row=6, column=1, padx=10, pady=10)
        self.init_h_state_combo.set(self.defaults["init_h_states"][0])

        # --------------------------------------------------------------

    def get_me_solver_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        me_solver_kwargs = {
            "t_steps": float(self.t_steps_entry.get()),
            "t_end": float(self.t_end_entry.get()),
            "t_unit": self.t_unit_combo.get(),
            "init_e_states": [self.init_e_state_combo.get()],
            "init_h_states": [self.init_h_state_combo.get()],
        }
        return me_solver_kwargs


# ----------------------------------------------------------------------
