import json
import multiprocessing
import time
import tkinter as tk
from tkinter import filedialog, messagebox

import customtkinter as ctk
import pandas as pd

from qDNA import (
    calc_dipole_dict,
    calc_dipole_moment_dict,
    calc_exciton_transfer_dict,
    calc_lifetime_dict,
)
from .scrollable_console_frame import ScrollableConsoleFrame


def format_time(seconds):
    if seconds >= 3600:
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        return f"{hours}h {minutes}m {seconds}s" if minutes or seconds else f"{hours}h"
    elif seconds >= 60:
        minutes = seconds // 60
        seconds = seconds % 60
        return f"{minutes}m {seconds}s" if seconds else f"{minutes}m"
    else:
        return f"{seconds}s"


class FastaWindow(ctk.CTkToplevel):
    def __init__(self, master):
        super().__init__(master)

        self.title("Fasta Input")
        self.setup_ui()

        self.upper_strands = []
        self.identifiers = []

        self.fasta_data = []
        self.kwargs = {}
        self.tb_model_name = "ELM"
        self.filename = "delete"
        self.num_cpu = multiprocessing.cpu_count() - 1

    def setup_ui(self):
        self.logo_label = ctk.CTkLabel(
            self, text="QuantumDNA", font=ctk.CTkFont(size=20, weight="bold")
        )
        self.logo_label.grid(row=0, column=0, pady=10, padx=10)

        # Upload buttons
        self.upload_fasta_button = ctk.CTkButton(
            self, text="Upload FASTA file", command=self.open_fasta_file
        )
        self.upload_fasta_button.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        self.upload_kwargs_button = ctk.CTkButton(
            self, text="Upload Options file", command=self.open_kwargs_file
        )
        self.upload_kwargs_button.grid(row=2, column=0, padx=10, pady=10, sticky="w")

        # Calculation options
        self.lifetime_var = ctk.BooleanVar(value=True)
        self.lifetime_checkbox = ctk.CTkCheckBox(
            self,
            text="Calculate Exciton Lifetime",
            state="disabled",
            variable=self.lifetime_var,
        )
        self.lifetime_checkbox.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        self.dipole_var = ctk.BooleanVar(value=True)
        self.dipole_checkbox = ctk.CTkCheckBox(
            self, text="Calculate Average Charge Separation", variable=self.dipole_var
        )
        self.dipole_checkbox.grid(row=1, column=1, padx=10, pady=10, sticky="w")

        self.dipole_moment_var = ctk.BooleanVar(value=True)
        self.dipole_moment_checkbox = ctk.CTkCheckBox(
            self, text="Calculate Dipole Moment", variable=self.dipole_moment_var
        )
        self.dipole_moment_checkbox.grid(row=2, column=1, padx=10, pady=10, sticky="w")

        self.exciton_transfer_var = ctk.BooleanVar(value=True)
        self.exciton_transfer_checkbox = ctk.CTkCheckBox(
            self,
            text="Calculate Average Exciton Population",
            variable=self.exciton_transfer_var,
        )
        self.exciton_transfer_checkbox.grid(
            row=3, column=1, padx=10, pady=10, sticky="w"
        )

        self.estimate_comp_time_button = ctk.CTkButton(
            self, text="Estimate Computation Time", command=self.estimate_comp_time
        )
        self.estimate_comp_time_button.grid(
            row=4, column=0, columnspan=2, padx=10, pady=10, sticky="ew"
        )

        # Submit button
        self.submit_button = ctk.CTkButton(
            self, text="Submit", command=self.process_files
        )
        self.submit_button.grid(
            row=5, column=0, columnspan=2, padx=10, pady=10, sticky="ew"
        )

        self.scrollable_console_frame = ScrollableConsoleFrame(self)
        self.scrollable_console_frame.grid(
            row=6, column=0, columnspan=2, padx=10, pady=10, sticky="nsew"
        )

    def open_fasta_file(self):
        fasta_file_path = filedialog.askopenfilename(
            title="Select FASTA file", filetypes=[("Text files", "*.txt")]
        )
        if not fasta_file_path:
            messagebox.showwarning("Warning", "No FASTA file selected.")
            return

        try:
            with open(fasta_file_path, "r") as file:
                self.fasta_data = file.readlines()
                self.parse_fasta_data()
            if not self.fasta_data:
                raise ValueError("FASTA file is empty.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open FASTA file: {e}")
            self.fasta_data = []

    def open_kwargs_file(self):
        kwargs_file_path = filedialog.askopenfilename(
            title="Select Options file", filetypes=[("JSON files", "*.json")]
        )
        if not kwargs_file_path:
            messagebox.showwarning("Warning", "No Options file selected.")
            return

        try:
            with open(kwargs_file_path, "r") as json_file:
                self.kwargs = json.load(json_file)
            if not self.kwargs:
                raise ValueError("Options file is empty.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open Options file: {e}")
            self.kwargs = {}

    def process_files(self):
        if not self.fasta_data:
            messagebox.showwarning("Warning", "Please provide a FASTA .txt file")
            return

        fasta_data = {
            "Identifier": self.identifiers,
            "Sequence (5'-3')": self.upper_strands,
        }

        try:
            self.tb_model_name = self.kwargs["tb_model_name"]
            del self.kwargs["tb_model_name"]
        except:
            pass

        start_time = time.time()
        if self.lifetime_var.get():
            lifetime_dict = calc_lifetime_dict(
                self.upper_strands, self.tb_model_name, self.filename, **self.kwargs
            )
            fasta_data["Exciton Lifetime (fs)"] = list(lifetime_dict.values())
        end_time_lifetime = time.time()

        if self.dipole_var.get():
            dipole_dict = calc_dipole_dict(self.tb_model_name, self.filename)
            fasta_data["Charge Separation (A)"] = list(dipole_dict.values())

        if self.dipole_moment_var.get():
            dipole_moment_dict = calc_dipole_moment_dict(
                self.tb_model_name, self.filename
            )
            fasta_data["Dipole Moment (D)"] = list(dipole_moment_dict.values())

        if self.exciton_transfer_var.get():
            exciton_transfer_dict = calc_exciton_transfer_dict(
                self.tb_model_name, self.filename
            )
            fasta_data["Exciton Transfer (upper strand)"] = [
                val[0]["exciton"] for val in exciton_transfer_dict.values()
            ]
            fasta_data["Exciton Transfer (lower strand)"] = [
                val[1]["exciton"] for val in exciton_transfer_dict.values()
            ]
        end_time = time.time()

        print("-------------------")
        print(fasta_data)
        print("-------------------")
        df_fasta = pd.DataFrame(fasta_data)
        self.kwargs["Computation Time (s)"] = end_time - start_time
        self.kwargs["Computation Time Lifetime (s)"] = end_time_lifetime - start_time
        metadata_df = pd.DataFrame(
            list(self.kwargs.items()), columns=["Parameter", "Value"]
        )

        self.save_to_excel(df_fasta, metadata_df)

    def parse_fasta_data(self):
        for line in self.fasta_data:
            line = line.strip()
            if line.startswith(">"):
                self.identifiers.append(line)
            else:
                self.upper_strands.append(line)

    def save_to_excel(self, df_fasta, metadata_df):
        save_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")]
        )
        if not save_path:
            messagebox.showwarning("Warning", "No save path specified.")
            return

        try:
            with pd.ExcelWriter(save_path) as writer:
                df_fasta.to_excel(writer, sheet_name="Data", index=False)
                metadata_df.to_excel(writer, sheet_name="Metadata", index=False)

            messagebox.showinfo("Success", f"File saved to {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {e}")

    def estimate_comp_time(self):
        comp_time_dict = {
            "WM": [
                0.04849529266357422,
                0.04625391960144043,
                0.07857656478881836,
                0.13587403297424316,
                0.21225285530090332,
                0.3305344581604004,
            ],
            "LM": [
                0.2975034713745117,
                1.308180570602417,
                3.8845632076263428,
                9.669641017913818,
                19.188920259475708,
                35.40109610557556,
            ],
            "ELM": [
                0.30926942825317383,
                1.379915475845337,
                4.239468812942505,
                11.73042607307434,
                25.32820200920105,
                52.95796084403992,
            ],
            "FWM": [
                8.07165561,
                34.5582086,
                139.004404,
                364.221048,
                1015.64643,
                2094.59942,
            ],
            "FLM": [
                22.8213988,
                133.818269,
                523.606570,
                1572.64914,
                3480.68617,
                5173.29872,
            ],
            "FELM": [
                23.9733094,
                140.136105,
                553.018477,
                1760.12157,
                4562.22129,
                5655.94587,
            ],
        }
        comp_time = 0
        for upper_strand in self.upper_strands:
            comp_time += comp_time_dict[self.tb_model_name][len(upper_strand) - 2]
        comp_time /= min(self.num_cpu, len(self.upper_strands))
        print(
            "Estimated Computation Time "
            + format_time(comp_time)
            + f"\nUsing {self.num_cpu} kernels"
        )
