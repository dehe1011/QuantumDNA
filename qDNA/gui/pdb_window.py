# pylint: skip-file
import customtkinter as ctk
from tkinter import filedialog, messagebox

from ..io import OPTIONS
from ..lcao import Oligomer
from .plotting_window import PlottingWindow2

# ----------------------------------------------------------------------


class PDBFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the custom window. This means that the custom window is its master.
            The frame uses the commands save() and cancel() from the master

        Widgets with get() method:
            tb_params_entry, source_entry, particle_combobox, tb_model_combobox, unit_combobox, methylation_checkbox
        """
        super().__init__(master)

        # initialization of the ct

        self.logo_label = ctk.CTkLabel(
            self, text="QuantumDNA", font=ctk.CTkFont(size=20, weight="bold")
        )
        self.logo_label.grid(row=0, column=0, pady=10, padx=10)

        self.upload_pdb_button = ctk.CTkButton(
            self, text="Upload PDB file", command=master.open_pdb_file
        )
        self.upload_pdb_button.grid(row=1, column=0, padx=10, pady=10)

        # param id label and combobox
        self.param_id_label = ctk.CTkLabel(self, text="LCAO Parametrization:")
        self.param_id_label.grid(row=1, column=1, padx=10, pady=0)

        self.param_id_combobox = ctk.CTkComboBox(self, values=OPTIONS["param_ids"])
        self.param_id_combobox.grid(row=2, column=1, padx=10, pady=10)
        self.param_id_combobox.set("MSF")

        # tb model name label and combobox
        self.tb_model_label = ctk.CTkLabel(self, text="TB Model:")
        self.tb_model_label.grid(row=3, column=1, padx=10, pady=0)

        self.tb_model_combobox = ctk.CTkComboBox(self, values=OPTIONS["tb_models"])
        self.tb_model_combobox.grid(row=4, column=1, padx=10, pady=10)
        self.tb_model_combobox.set("ELM")

        self.text_label = ctk.CTkLabel(
            self,
            text="(Tip: Restart the session to apply and use \nthe calculated parameters in simulations.)",
        )
        self.text_label.grid(row=6, column=0, columnspan=2, pady=10, padx=10)

        self.cancel_button = ctk.CTkButton(self, text="Cancel", command=master.cancel)
        self.cancel_button.grid(row=7, column=0, padx=10, pady=10)

        self.save_button = ctk.CTkButton(self, text="Save", command=master.save)
        self.save_button.grid(row=7, column=1, padx=10, pady=10)

        # plot label and particle combobox
        self.plot_label = ctk.CTkLabel(self, text="Plot TB parameters", font=ctk.CTkFont(size=15))
        self.plot_label.grid(row=8, column=0, columnspan=2, pady=10, padx=10)

        self.particle_combobox = ctk.CTkComboBox(
            self,
            values=OPTIONS["particles"],
        )
        self.particle_combobox.grid(row=9, column=0, padx=10, pady=10)
        self.particle_combobox.set("hole")

        self.plot_button = ctk.CTkButton(self, text="Plot", command=master.plot)
        self.plot_button.grid(row=9, column=1, padx=10, pady=10)


class PDBWindow(ctk.CTkToplevel):
    def __init__(self, master):
        """
        Notes:
            This window is a toplevel window of the main window. This means that the main window is its master.
            The frame uses the commands save() and cancel() from the master
        """

        # initialization of the ctk.CTkToplevel class
        super().__init__(master)
        self.title("PDB Input")

        # add a instance of CustomFrame to the window
        self.custom_frame = PDBFrame(self)
        self.custom_frame.grid(row=0, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

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
        self.param_id = self.custom_frame.param_id_combobox.get()
        self.tb_model_name = self.custom_frame.tb_model_combobox.get()

    def cancel(self):
        """Closes the window."""
        self.destroy()

    def plot(self):
        self.particle = self.custom_frame.particle_combobox.get()
        self.plotting_window = PlottingWindow2(self)

    def save(self):
        """Saves the tight-binding parameters and closes the window."""
        self.get_custom_frame_params()

        self.oligomer = Oligomer(
            self.pdb_file_path, param_id=self.param_id, tb_model_name=self.tb_model_name
        )
        self.oligomer.save_tb_params()
        print(
            f"TB parameters saved successfully. Access them as {self.oligomer.filename_pdb}_{self.oligomer.tb_model_name}."
        )
        # self.destroy()


# ----------------------------------------------------------------------
