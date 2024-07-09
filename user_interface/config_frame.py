import customtkinter as ctk
import webbrowser
from .custom_window import CustomWindow

class ConfigFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        padx = 20
        self.configs = master.configs
        
        self.open_github_button = ctk.CTkButton(self, text="GitHub", command=self.open_github)
        self.open_github_button.grid(row=8, column=0, pady=10, padx=padx)
    
        self.open_custom_window_button = ctk.CTkButton(self, text="Customize", command=self.open_custom_window)
        self.open_custom_window_button.grid(row=9, column=0, pady=10, padx=padx)
    
        self.appearance_mode_label = ctk.CTkLabel(self, text="Appearance Mode:", anchor="w")
        self.appearance_mode_label.grid(row=10, column=0, pady=10, padx=padx)
        
        self.appearance_mode_combo = ctk.CTkComboBox(self, values=["Light", "Dark", "System"], command=ctk.set_appearance_mode)
        self.appearance_mode_combo.grid(row=11, column=0, pady=10, padx=padx)
        self.appearance_mode_combo.set('System')

    def open_github(self):
        webbrowser.open('https://github.com/dehe1011/quantum_DNA')

    def open_custom_window(self):
        self.custom_window = CustomWindow(self)