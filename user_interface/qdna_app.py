import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import customtkinter as ctk
from utils import get_config
from DNA import get_me_solver, TB_Ham, DNA_Seq
from user_interface import InitialFrame, ConfigFrame, OptionsFrame, ScrollableConsoleFrame, PlotOptionsFrame, PlottingWindow

class qDNA_App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("QuantumDNA")
        self.configs = get_config()
        self.tb_basis = ['(0, 0)']
        
        self.initial_frame = InitialFrame(self)
        self.initial_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        
        self.config_frame = ConfigFrame(self)
        self.config_frame.grid(row=1, column=0, padx=20, pady=20, sticky="nsew")

        self.options_frame = OptionsFrame(self)
        self.options_frame.grid(row=0, column=1, rowspan = 2, padx=20, pady=20, sticky="nsew")
        
        self.scrollable_console_frame = ScrollableConsoleFrame(self)
        self.scrollable_console_frame.grid(row=1, column=2, padx=20, pady=20, sticky="nsew")

        self.plot_options_frame = PlotOptionsFrame(self)
        self.plot_options_frame.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")

        self.options_frame.change_state('disabled')
        self.plot_options_frame.change_state('disabled')

    def press_first_confirm(self):
        self.enable_options_frame()
        self.upper_strand = self.initial_frame.upper_strand_entry.get()
        self.tb_model_name = self.initial_frame.tb_model_combo.get()
        self.tb_ham = TB_Ham(DNA_Seq(self.upper_strand, self.tb_model_name))
        self.tb_basis = self.tb_ham.tb_basis
        self.options_frame.options_tab.dynamics_frame.me_init_e_state_combo.configure(values=self.tb_basis)
        self.options_frame.options_tab.dynamics_frame.me_init_h_state_combo.configure(values=self.tb_basis)
        self.plot_options_frame.plotting_options_tab.pop_frame.tb_site_combo.configure(values=self.tb_basis)
    
    def enable_initial_frame(self):
        self.initial_frame.change_state('normal')
        self.options_frame.change_state('disabled')
        self.plot_options_frame.change_state('disabled')
        
    def enable_options_frame(self):
        self.initial_frame.change_state('disabled')
        self.options_frame.change_state('normal')
        self.plot_options_frame.change_state('disabled')

    def enable_plotting_frame(self):
        self.initial_frame.change_state('disabled')
        self.options_frame.change_state('disabled')
        self.plot_options_frame.change_state('normal')

    def submit(self):
        self.kwargs = self.options_frame.get_kwargs()
        self.plotting_kwargs = self.plot_options_frame.get_plotting_kwargs()

        self.me_solver = get_me_solver(self.upper_strand, self.tb_model_name, **self.kwargs)
        self.plotting_window = PlottingWindow(self)  

if __name__ == "__main__":
    app = qDNA_App()
    app.mainloop()