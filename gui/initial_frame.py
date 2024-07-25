import customtkinter as ctk


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

        self.upper_strand_label = ctk.CTkLabel(
            self, text="DNA Sequence \n (upper strand):"
        )
        self.upper_strand_label.grid(row=1, column=0, pady=10, padx=10)

        self.upper_strand_entry = ctk.CTkEntry(self)
        self.upper_strand_entry.grid(row=2, column=0, pady=10, padx=10)
        self.upper_strand_entry.insert(0, "GCG")

        self.tb_model_label = ctk.CTkLabel(self, text="TB Model:")
        self.tb_model_label.grid(row=3, column=0, pady=10, padx=10)

        self.tb_model_combo = ctk.CTkComboBox(self, values=master.configs["TB_MODELS"])
        self.tb_model_combo.grid(row=4, column=0, pady=10, padx=10)
        self.tb_model_combo.set("WM")

        self.first_confirm_button = ctk.CTkButton(
            self, text="Confirm", command=master.press_first_confirm
        )
        self.first_confirm_button.grid(row=5, column=0, pady=10, padx=10)

    def change_state(self, state):
        """
        Changes the state of certain widgets (between 'normal' and 'disabled').
        """
        self.upper_strand_entry.configure(state=state)
        self.tb_model_combo.configure(state=state)
        self.first_confirm_button.configure(state=state)
