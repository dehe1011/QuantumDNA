# pylint: skip-file

import customtkinter as ctk

# --------------------------------------------------


class ConfigFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the main/ menu window. This means that the main window is its master.
            The frame uses the commands open_pdb_window(), open_fasta_window() from the master
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # widgets
        self.advanced_label = ctk.CTkLabel(
            self, text="Advanced", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.advanced_label.grid(row=0, column=0, pady=10, padx=10)

        self.open_fasta_window_button = ctk.CTkButton(
            self, text="FASTA Input", command=master.open_fasta_window
        )
        self.open_fasta_window_button.grid(row=1, column=0, pady=10, padx=10)

        self.open_pdb_window_button = ctk.CTkButton(
            self, text="PDB Input", command=master.open_pdb_window
        )
        self.open_pdb_window_button.grid(row=2, column=0, pady=10, padx=10)


class HelpFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the main/ menu window. This means that the main window is its master.
            The frame uses the commands open_github(), open_documentation(), open_tutorials() from the master
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # widgets
        self.help_label = ctk.CTkLabel(
            self, text="Help", font=ctk.CTkFont(size=15, weight="bold")
        )
        self.help_label.grid(row=0, column=0, pady=10, padx=10)

        self.open_github_button = ctk.CTkButton(
            self, text="GitHub", command=master.open_github
        )
        self.open_github_button.grid(row=1, column=0, pady=10, padx=10)

        self.open_github_button = ctk.CTkButton(
            self, text="Documentation", command=master.open_documentation
        )
        self.open_github_button.grid(row=2, column=0, pady=10, padx=10)

        self.open_github_button = ctk.CTkButton(
            self, text="Tutorials", command=master.open_tutorials
        )
        self.open_github_button.grid(row=3, column=0, pady=10, padx=10)
