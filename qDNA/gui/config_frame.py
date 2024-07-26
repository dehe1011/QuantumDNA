import customtkinter as ctk


class ConfigFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        """
        Notes:
            This frame is located in the main/ menu window. This means that the main window is its master.
            The frame uses the commands open_custom_window() and open_github() from the master

        Widgets with get() method:
            appearance_mode_combo
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # widgets
        self.open_github_button = ctk.CTkButton(
            self, text="GitHub", command=master.open_github
        )
        self.open_github_button.grid(row=8, column=0, pady=10, padx=10)

        self.open_custom_window_button = ctk.CTkButton(
            self, text="Customize", command=master.open_custom_window
        )
        self.open_custom_window_button.grid(row=9, column=0, pady=10, padx=10)

        self.appearance_mode_label = ctk.CTkLabel(
            self, text="Appearance Mode:", anchor="w"
        )
        self.appearance_mode_label.grid(row=10, column=0, pady=10, padx=10)

        self.appearance_mode_combo = ctk.CTkComboBox(
            self, values=["Light", "Dark", "System"], command=ctk.set_appearance_mode
        )
        self.appearance_mode_combo.grid(row=11, column=0, pady=10, padx=10)
        self.appearance_mode_combo.set("System")
