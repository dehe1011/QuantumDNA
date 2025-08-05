# pylint: skip-file

import customtkinter as ctk

from .tb_ham_frame import TBHamFrame
from .lind_diss_frame import LindDissFrame
from .me_solver_frame import MeSolverFrame
from .gui_utils import change_state_all_widgets

# ----------------------------------------------------------------------


class OptionsTab(ctk.CTkTabview):
    def __init__(self, master):
        """
        Notes:
            This tab is located in the options frame. This means that the options frame is its master.
            The tab itself is the master for three further frames.
        """

        # initialization of the ctk.CTkTabview class
        super().__init__(master)

        self.tb_ham_tab = self.add("Hamiltonian")
        self.lind_diss_tab = self.add("Dissipator")
        self.me_solver_tab = self.add("Dynamics")
        # self.set("Hamiltonian")

        self.tb_ham_frame = TBHamFrame(self.tb_ham_tab)
        self.lind_diss_frame = LindDissFrame(self.lind_diss_tab)
        self.me_solver_frame = MeSolverFrame(self.me_solver_tab)


class OptionsFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the main window. This means that the main window is its master.
            The frame itself is the master for the tab options_tab.
            The frame uses the press_second_confirm(), enable_initial_frame() from the master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # options label
        self.label = ctk.CTkLabel(self, text="Options", font=ctk.CTkFont(size=20, weight="bold"))
        self.label.grid(row=0, column=0, pady=10, padx=10, columnspan=2)

        # options tab
        self.options_tab = OptionsTab(self)
        self.options_tab.grid(row=1, column=0, columnspan=2, pady=10, padx=10)

        # back button
        self.back_button = ctk.CTkButton(self, text="Back", command=master.enable_initial_frame)
        self.back_button.grid(row=2, column=0, pady=10, padx=10)

        # second confirm button
        self.second_confirm_button = ctk.CTkButton(
            self, text="Confirm", command=master.press_second_confirm
        )
        self.second_confirm_button.grid(row=2, column=1, pady=10, padx=10)

    def change_state(self, state):
        """Changes the state of certain widgets (between 'normal' and 'disabled')."""
        change_state_all_widgets(self.options_tab.tb_ham_frame, state=state)
        change_state_all_widgets(self.options_tab.lind_diss_frame, state=state)
        change_state_all_widgets(self.options_tab.me_solver_frame, state=state)
        self.second_confirm_button.configure(state=state)
        self.back_button.configure(state=state)


# ----------------------------------------------------------------------
