# pylint: skip-file

import customtkinter as ctk

# --------------------------------------------------


class InitialFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        """
        Notes:
            This frame is located in the main/ menu window. This means that the main window is its master.
            The frame uses the commands press_first_confirm() from the master.

        Widgets with get() method:
            upper_strand_entry, tb_model_combo
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # widgets
        self.logo_label = ctk.CTkLabel(
            self, text="QuantumDNA", font=ctk.CTkFont(size=20, weight="bold")
        )
        self.logo_label.grid(row=0, column=0, pady=10, padx=10)

        self.grid_rowconfigure(1, weight=1)

        self.upper_strand_label = ctk.CTkLabel(
            self, text="Upper DNA Strand \n(5'-3' direction):"
        )
        self.upper_strand_label.grid(row=2, column=0, pady=0, padx=10)

        self.upper_strand_entry = ctk.CTkEntry(self)
        self.upper_strand_entry.grid(row=3, column=0, pady=10, padx=10)
        self.upper_strand_entry.insert(0, "G_C_G")

        self.lower_strand_label = ctk.CTkLabel(
            self, text="Lower DNA Strand \n(3'-5' direction):"
        )
        self.lower_strand_label.grid(row=4, column=0, pady=0, padx=10)

        self.lower_strand_entry = ctk.CTkEntry(self)
        self.lower_strand_entry.grid(row=5, column=0, pady=10, padx=10)
        self.lower_strand_entry.insert(0, "auto complete")

        self.tb_model_label = ctk.CTkLabel(self, text="TB Model:")
        self.tb_model_label.grid(row=6, column=0, pady=0, padx=10)

        self.tb_model_combo = ctk.CTkComboBox(self, values=master.configs["TB_MODELS"])
        self.tb_model_combo.grid(row=7, column=0, pady=10, padx=10)
        self.tb_model_combo.set("ELM")

        self.grid_rowconfigure(8, weight=1)

        self.first_confirm_button = ctk.CTkButton(
            self, text="Confirm", command=master.press_first_confirm
        )
        self.first_confirm_button.grid(row=9, column=0, pady=10, padx=10)

    def change_state(self, state):
        """Changes the state of certain widgets (between 'normal' and 'disabled')."""
        self.upper_strand_entry.configure(state=state)
        self.tb_model_combo.configure(state=state)
        self.first_confirm_button.configure(state=state)
