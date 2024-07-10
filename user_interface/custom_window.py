import customtkinter as ctk
import webbrowser
from DNA import wrap_save_tb_params

class CustomFrame(ctk.CTkFrame):
    def __init__(self, master, configs):
        super().__init__(master)
        self.configs = configs
        
        self.tb_params_label = ctk.CTkLabel(self, text="TB parameters:")
        self.tb_params_label.grid(row=0, column=0, columnspan=2, pady=10, padx=10)
    
        self.tb_params_entry = ctk.CTkEntry(self)
        self.tb_params_entry.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

        self.source_label = ctk.CTkLabel(self, text="Source:")
        self.source_label.grid(row=2, column=0, columnspan=2, padx=10, pady=10)
    
        self.source_entry = ctk.CTkEntry(self)
        self.source_entry.grid(row=3, column=0, columnspan=2, padx=10, pady=10)
        
        self.particle_label = ctk.CTkLabel(self, text="Particle:")
        self.particle_label.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

        self.particle_combobox = ctk.CTkComboBox(self, values=self.configs['PARTICLES'] )
        self.particle_combobox.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

        self.tb_model_label = ctk.CTkLabel(self, text="TB Model:")
        self.tb_model_label.grid(row=6, column=0, columnspan=2, padx=10, pady=10)

        self.tb_model_combobox = ctk.CTkComboBox(self, values=self.configs['TB_MODELS'] )
        self.tb_model_combobox.grid(row=7, column=0, columnspan=2, padx=10, pady=10)

        self.unit_label = ctk.CTkLabel(self, text="Unit:")
        self.unit_label.grid(row=8, column=0, columnspan=2, padx=10, pady=10)

        self.unit_combobox = ctk.CTkComboBox(self, values=self.configs['UNITS'] )
        self.unit_combobox.grid(row=9, column=0, columnspan=2, padx=10, pady=10)

        check_var = ctk.BooleanVar(value=True)
        self.methylation_checkbox = ctk.CTkCheckBox(self, text="Methylation", variable=check_var, onvalue=True, offvalue=False)
        self.methylation_checkbox.grid(row=10, column=0, columnspan=2, padx=10, pady=10)
       
        self.save_button = ctk.CTkButton(self, text="Save", command=master.save)
        self.save_button.grid(row=11, column=0, padx=10, pady=10)
        
        self.cancel_button = ctk.CTkButton(self, text="Cancel", command=master.cancel)
        self.cancel_button.grid(row=11, column=1, padx=10, pady=10)

class CustomWindow(ctk.CTkToplevel):
    def __init__(self, master, configs):
        super().__init__(master)
        self.configs = configs
        self.title("Customize")
        
        self.custom_frame = CustomFrame(self, configs)
        self.custom_frame.grid(row=0, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

    def save(self):
        tb_param_dict = self.custom_frame.tb_params_entry.get()
        source = self.custom_frame.source_entry.get()
        particle = self.custom_frame.particle_combobox.get()
        tb_model_name = self.custom_frame.tb_model_combobox.get()
        unit = self.custom_frame.unit_combobox.get()
        notes = f"methlation {self.custom_frame.methylation_checkbox.get()}"
        wrap_save_tb_params(tb_param_dict, source, particle, tb_model_name, unit, notes=notes)
            
        self.destroy()
    
    def cancel(self):
        self.destroy()