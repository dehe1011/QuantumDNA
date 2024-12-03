# pylint: skip-file
import os

import customtkinter as ctk
from tkinter import filedialog, messagebox

from .. import DATA_DIR
from ..lcao import convert_pdb_to_xyz, calc_tb_params
from ..hamiltonian import wrap_save_tb_params

# --------------------------------------------------


class PDBFrame(ctk.CTkFrame):
    def __init__(self, master, configs):
        """
        Notes:
            This frame is located in the custom window. This means that the custom window is its master.
            The frame uses the commands save() and cancel() from the master

        Widgets with get() method:
            tb_params_entry, source_entry, particle_combobox, tb_model_combobox, unit_combobox, methylation_checkbox
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.configs = configs

        # widgets

        # Upload button

        self.logo_label = ctk.CTkLabel(
            self, text="QuantumDNA", font=ctk.CTkFont(size=20, weight="bold")
        )
        self.logo_label.grid(row=0, column=0, pady=10, padx=10)

        self.upload_pdb_button = ctk.CTkButton(
            self, text="Upload PDB file", command=master.open_pdb_file
        )
        self.upload_pdb_button.grid(row=1, column=0, padx=10, pady=10)

        self.identifier_label = ctk.CTkLabel(self, text="Identifier:")
        self.identifier_label.grid(row=0, column=1, padx=10, pady=0)

        self.identifier_entry = ctk.CTkEntry(self)
        self.identifier_entry.grid(row=1, column=1, padx=10, pady=10)

        self.tb_model_label = ctk.CTkLabel(self, text="TB Model:")
        self.tb_model_label.grid(row=2, column=1, padx=10, pady=0)

        self.tb_model_combobox = ctk.CTkComboBox(self, values=self.configs["TB_MODELS"])
        self.tb_model_combobox.grid(row=3, column=1, padx=10, pady=10)

        self.notes_label = ctk.CTkLabel(self, text="Notes:")
        self.notes_label.grid(row=4, column=1, padx=10, pady=0)

        self.notes_entry = ctk.CTkEntry(self)
        self.notes_entry.grid(row=5, column=1, padx=10, pady=10)

        self.text_label = ctk.CTkLabel(
            self,
            text="(Tip: Restart the session to apply and use \nthe calculated parameters in simulations.)",
        )
        self.text_label.grid(row=6, column=0, columnspan=2, pady=10, padx=10)

        self.save_button = ctk.CTkButton(self, text="Save", command=master.save)
        self.save_button.grid(row=7, column=0, padx=10, pady=10)

        self.cancel_button = ctk.CTkButton(self, text="Cancel", command=master.cancel)
        self.cancel_button.grid(row=7, column=1, padx=10, pady=10)


class PDBWindow(ctk.CTkToplevel):
    def __init__(self, master, configs):
        """
        Notes:
            This window is a toplevel window of the main window. This means that the main window is its master.
            The frame uses the commands save() and cancel() from the master
        """

        # initialization of the ctk.CTkToplevel class
        super().__init__(master)
        self.title("PDB Input")

        # add a instance of CustomFrame to the window
        self.custom_frame = PDBFrame(self, configs)
        self.custom_frame.grid(
            row=0, column=0, columnspan=2, padx=10, pady=10, sticky="nsew"
        )

    def open_pdb_file(self):
        pdb_file_path = filedialog.askopenfilename(
            title="Select PDB file", filetypes=[("Text files", "*.pdb")]
        )
        if not pdb_file_path:
            messagebox.showwarning("Warning", "No PDB file selected.")
            return

        self.pdb_file_path = pdb_file_path

    def get_custom_frame_params(self):
        """Makes all parameters of custom_frame available in this window."""
        self.identifier = self.custom_frame.identifier_entry.get()
        self.tb_model_name = self.custom_frame.tb_model_combobox.get()
        self.notes = self.custom_frame.notes_entry.get()

    def save(self):
        """Saves the tight-binding parameters and closes the window."""
        self.get_custom_frame_params()

        convert_pdb_to_xyz(self.pdb_file_path)
        # calculate tight-binding parameters
        directories = [self.pdb_file_path.split(".")[0]]
        HOMO_dict, LUMO_dict = calc_tb_params(directories, self.tb_model_name)

        wrap_save_tb_params(
            HOMO_dict,
            self.identifier,
            "hole",
            self.tb_model_name,
            unit="meV",
            notes=self.notes,
        )
        wrap_save_tb_params(
            LUMO_dict,
            self.identifier,
            "electron",
            self.tb_model_name,
            unit="meV",
            notes=self.notes,
        )
        self.destroy()

    def cancel(self):
        """Closes the window."""
        self.destroy()
