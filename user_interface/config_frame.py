import customtkinter as ctk

class ConfigFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master)
        padx = 20
        self.configs = master.configs
        
        self.open_github_button = ctk.CTkButton(self, text="GitHub", command=master.open_github)
        self.open_github_button.grid(row=8, column=0, pady=10, padx=padx)
    
        self.open_custom_window_button = ctk.CTkButton(self, text="Customize", command=master.open_custom_window)
        self.open_custom_window_button.grid(row=9, column=0, pady=10, padx=padx)
    
        self.appearance_mode_label = ctk.CTkLabel(self, text="Appearance Mode:", anchor="w")
        self.appearance_mode_label.grid(row=10, column=0, pady=10, padx=padx)
        
        self.appearance_mode_combo = ctk.CTkComboBox(self, values=["Light", "Dark", "System"], command=ctk.set_appearance_mode)
        self.appearance_mode_combo.grid(row=11, column=0, pady=10, padx=padx)
        self.appearance_mode_combo.set('System')
