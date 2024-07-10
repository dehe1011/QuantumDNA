import customtkinter as ctk
from itertools import product

from DNA import calc_lifetime, calc_dipole
from .user_interface_utils import change_state_all_widgets

class LifetimeFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.kwargs = kwargs
        
        self.lifetime_button = ctk.CTkButton(self, text="Calculate exciton lifetime", command = self.calc_lifetime)
        self.lifetime_button.grid(row=0, column=0, padx=10, pady=10)
        self.dipole_button = ctk.CTkButton(self, text="Calculate charge separation", command = self.calc_dipole)
        self.dipole_button.grid(row=1, column=0, padx=10, pady=10)
        
    def calc_lifetime(self):
        lifetime = calc_lifetime(**self.kwargs)
        print(f"Exciton Lifetime: {lifetime} fs")
    
    def calc_dipole(self):
        dipole = calc_dipole(**self.kwargs)
        print(f"Charge separation: {dipole} A") 
 
class PopFrame(ctk.CTkFrame):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master)
        self.pack(fill="both", expand=True)

        self.tb_site_label = ctk.CTkLabel(self, text="TB site:")
        self.tb_site_label.grid(row=0, column=0, padx=10, pady=10)

        self.tb_site_combo = ctk.CTkComboBox(self, values=['All']+tb_basis)
        self.tb_site_combo.grid(row=0, column=1, padx=10, pady=10)
        self.tb_site_combo.set('All')

    def get_pop_kwargs(self):
        pop_kwargs = {"init_tb_site": self.tb_site_combo.get(),
        }
        return pop_kwargs
    
class CohFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master)
        self.pack(fill="both", expand=True)
        
    def get_coh_kwargs(self):
        coh_kwargs = {}
        return coh_kwargs
    
class FourierFrame(ctk.CTkFrame):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master)
        self.pack(fill="both", expand=True)

        self.tb_site_label = ctk.CTkLabel(self, text="TB site:")
        self.tb_site_label.grid(row=0, column=0, padx=10, pady=10)
        
        self.tb_site_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.tb_site_combo.grid(row=0, column=1, padx=10, pady=10)
        
        self.init_e_site_label = ctk.CTkLabel(self, text="Initial electron \n TB site:")
        self.init_e_site_label.grid(row=1, column=0, padx=10, pady=10)
        
        self.init_e_site_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.init_e_site_combo.grid(row=1, column=1, padx=10, pady=10)
        
        self.init_h_site_label = ctk.CTkLabel(self, text="Initial hole \n TB site:")
        self.init_h_site_label.grid(row=2, column=0, padx=10, pady=10)
        
        self.init_h_site_combo = ctk.CTkComboBox(self, values=tb_basis)
        self.init_h_site_combo.grid(row=2, column=1, padx=10, pady=10)
        
        self.x_axis_combo = ctk.CTkComboBox(self, values=['period', 'frequency'])
        self.x_axis_combo.grid(row=3, column=0, padx=10, pady=10)
        self.x_axis_combo.set('period')
              
        self.average_pop_var = ctk.BooleanVar(value=False)
        self.average_pop_check = ctk.CTkCheckBox(self, text="Calculate average population", variable=self.average_pop_var)
        self.average_pop_check.grid(row=4, column=0, columnspan=2, padx=10, pady=10)
             
    def get_fourier_kwargs(self):
        fourier_kwargs = {"end_state": self.tb_site_combo.get(),
                      "init_state": ( self.init_e_site_combo.get(), self.init_h_site_combo.get() ),
                      "x_axis": self.x_axis_combo.get(),
        }
        return fourier_kwargs
    
# ----------------------------------------------
    
class PlotOptionsTab(ctk.CTkTabview):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master)
        
        self.pop_tab = self.add("Population")
        self.coh_tab = self.add("Coherence")
        self.fourier_tab = self.add("Fourier")
        self.lifetime_tab = self.add("Lifetime")
        self.set("Population")

        self.pop_frame = PopFrame(self.pop_tab, tb_basis, **kwargs)
        self.coh_frame = CohFrame(self.coh_tab, **kwargs)
        self.fourier_frame = FourierFrame(self.fourier_tab, tb_basis, **kwargs)
        self.lifetime_frame = LifetimeFrame(self.lifetime_tab, **kwargs)

class PlotOptionsFrame(ctk.CTkFrame):
    def __init__(self, master, tb_basis, **kwargs):
        super().__init__(master)

        self.label = ctk.CTkLabel(self, text="Plotting", font=ctk.CTkFont(size=20, weight="bold"))
        self.label.grid(row=0, column=0, columnspan=2, pady=10, padx=10)

        self.plotting_options_tab = PlotOptionsTab(self, tb_basis, **kwargs)
        self.plotting_options_tab.grid(row=1, column=0, columnspan=2, pady=10, padx=10)

        self.back_button = ctk.CTkButton(self, text="Back", command=master.enable_options_frame)
        self.back_button.grid(row=3, column=0, pady=10, padx=10)
        
        self.submit_button = ctk.CTkButton(self, text="Submit", command=master.submit)
        self.submit_button.grid(row=3, column=1, pady=10, padx=10)
      
    def change_state(self, state):
        change_state_all_widgets(self.plotting_options_tab.pop_frame, state=state)
        change_state_all_widgets(self.plotting_options_tab.coh_frame, state=state)
        change_state_all_widgets(self.plotting_options_tab.fourier_frame, state=state)
        change_state_all_widgets(self.plotting_options_tab.lifetime_frame, state=state)
        self.submit_button.configure(state=state)
        self.back_button.configure(state=state)
        