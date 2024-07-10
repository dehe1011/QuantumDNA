import customtkinter as ctk
from .ham_frame import HamFrame
from .diss_frame import DissFrame
from .dynamics_frame import DynamicsFrame
from .user_interface_utils import change_state_all_widgets

class OptionsTab(ctk.CTkTabview):
    def __init__(self, master, configs, **kwargs):
        super().__init__(master)
       
        self.ham_tab = self.add("Hamiltonian")
        self.diss_tab = self.add("Dissipator")
        self.dynamics_tab = self.add("Dynamics")
        self.set("Hamiltonian")

        self.ham_frame = HamFrame(self.ham_tab, configs, **kwargs)
        self.diss_frame = DissFrame(self.diss_tab, configs, **kwargs)
        self.dynamics_frame = DynamicsFrame(self.dynamics_tab, configs, **kwargs)

class OptionsFrame(ctk.CTkFrame):
    def __init__(self, master, configs, **kwargs):
        super().__init__(master)
        
        self.label = ctk.CTkLabel(self, text="Options", font=ctk.CTkFont(size=20, weight="bold"))
        self.label.grid(row=0, column=0, pady=10, padx=10, columnspan=2)

        self.options_tab = OptionsTab(self, configs, **kwargs)
        self.options_tab.grid(row=1, column=0, columnspan=2, pady=10, padx=10)

        self.back_button = ctk.CTkButton(self, text="Back", command=master.enable_initial_frame)
        self.back_button.grid(row=2, column=0, pady=20, padx=10)

        self.second_confirm_button = ctk.CTkButton(self, text="Confirm", command=master.press_second_confirm)
        self.second_confirm_button.grid(row=2, column=1, pady=10, padx=10)

    def change_state(self, state):
        change_state_all_widgets(self.options_tab.ham_frame, state=state)
        change_state_all_widgets(self.options_tab.diss_frame, state=state)
        change_state_all_widgets(self.options_tab.dynamics_frame, state=state)
        self.second_confirm_button.configure(state=state)
        self.back_button.configure(state=state)
        