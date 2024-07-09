import customtkinter as ctk
from .user_interface_utils import change_state_all_widgets

class PopFrame(ctk.CTkFrame):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master, **kwargs)
        self.pack(fill="both", expand=True)

        self.tb_site_label = ctk.CTkLabel(self, text="TB site:")
        self.tb_site_label.grid(row=0, column=0, padx=10, pady=10)
        
        self.tb_site_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.tb_site_combo.grid(row=0, column=1, padx=10, pady=10)
        self.tb_site_combo.set('(0, 0)')

    def get_pop_kwargs(self):
        pop_kwargs = {"init_tb_site": self.tb_site_combo.get(),
        }
        return pop_kwargs

class PlotOptionsTab(ctk.CTkTabview):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master, **kwargs)
        
        self.pop_tab = self.add("Population")
        self.coh_tab = self.add("Coherence")
        self.fourier_tab = self.add("Fourier")
        self.set("Population")

        self.pop_frame = PopFrame(self.pop_tab, tb_basis, **kwargs)

class PlotOptionsFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)

        self.label = ctk.CTkLabel(self, text="Plotting", font=ctk.CTkFont(size=20, weight="bold"))
        self.label.grid(row=0, column=0, columnspan=2, pady=10, padx=10)

        self.plotting_options_tab = PlotOptionsTab(self, master.tb_basis)
        self.plotting_options_tab.grid(row=1, column=0, columnspan=2, pady=10, padx=10)

        self.back_button = ctk.CTkButton(self, text="Back", command=master.enable_options_frame)
        self.back_button.grid(row=3, column=0, pady=10, padx=10)
        
        self.submit_button = ctk.CTkButton(self, text="Submit", command=master.submit)
        self.submit_button.grid(row=3, column=1, pady=10, padx=10)

    def get_plotting_kwargs(self):
        self.pop_kwargs = self.plotting_options_tab.pop_frame.get_pop_kwargs() 
        return dict(**self.pop_kwargs)
      
    def change_state(self, state):
        change_state_all_widgets(self.plotting_options_tab.pop_frame, state=state)
        self.submit_button.configure(state=state)
        self.back_button.configure(state=state)
        