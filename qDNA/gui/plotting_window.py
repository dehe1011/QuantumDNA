# pylint: skip-file

import os

import customtkinter as ctk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ..io import save_figure

# -----------------------------------------------------------------------------------------------------


class PlottingFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the plotting_window. This means that the plotting_window is its master.
            The frame uses the commands save() and cancel() from the master

        Widgets with get() method:
            filename_entry
            directory_entry
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)

        # widgets
        self.filename_label = ctk.CTkLabel(self, text="Filename:")
        self.filename_label.grid(row=0, column=0, pady=10, padx=10)

        self.filename_entry = ctk.CTkEntry(self)
        self.filename_entry.grid(row=1, column=0, padx=10, pady=10)

        self.directory_label = ctk.CTkLabel(self, text="Directory:")
        self.directory_label.grid(row=0, column=1, pady=10, padx=10)

        self.directory_entry = ctk.CTkEntry(self)
        self.directory_entry.grid(row=1, column=1, padx=10, pady=10)

        self.subframe = ctk.CTkFrame(self)
        self.subframe.grid(row=3, column=0, columnspan=2, padx=10, pady=10)

        self.save_button = ctk.CTkButton(self, text="Save", command=master.save)
        self.save_button.grid(row=4, column=0, padx=10, pady=10)

        self.cancel_button = ctk.CTkButton(self, text="Cancel", command=master.cancel)
        self.cancel_button.grid(row=4, column=1, padx=10, pady=10)


class PlottingWindow(ctk.CTkToplevel):
    def __init__(self, master):
        """
        Notes:
            This window is a toplevel window of the main window. This means that the main window is its master.
            This window itself is the master for plotting_frame.
            The frame uses the command plot_options_kwargs, me_solver and tb_ham from the master
        """

        # initialization of the ctk.CTkToplevel class
        super().__init__(master)
        self.controller = master
        self.title("Plotting")

        self.plotting_frame = PlottingFrame(self)
        self.plotting_frame.grid(
            row=0, column=0, columnspan=2, padx=10, pady=10, sticky="nsew"
        )

        self.plot_option = self.controller.plot_kwargs["plot_option"]

        if self.plot_option == "Population":
            self.init_tb_site = self.controller.plot_kwargs["init_tb_site"]
            if self.init_tb_site == "Heatmap":
                self.plot_heatmap()
            elif self.init_tb_site == "All DNA Bases":
                self.plot_pops()
            else:
                self.plot_pop()

        if self.plot_option == "Coherence":
            self.plot_coh()

        if self.plot_option == "Spectrum":
            if self.controller.plot_kwargs["eigv_var"]:
                self.plot_eigv()
            if self.controller.plot_kwargs["eigs_var"]:
                self.eigenstate_idx = int(self.controller.plot_kwargs["eigenstate_idx"])
                self.plot_eigs()

        if self.plot_option == "Fourier":
            self.plot_fourier()
            if self.controller.plot_kwargs["average_pop_var"]:
                self.average_pop()
        self.plotting(self.fig)

    def average_pop(self):
        average_pop = self.controller.vis.get_average_pop(
            self.controller.plot_kwargs["init_state"],
            self.controller.plot_kwargs["end_state"],
        )
        print(f"Average population: {average_pop}")

    # ------------------------------------------------------------------

    def plot_eigv(self):
        self.fig, self.ax = self.controller.vis.plot_eigv()

    def plot_eigs(self):
        self.fig, self.ax = self.controller.vis.plot_eigs(self.eigenstate_idx)

    def plot_pop(self):
        self.fig, self.ax = self.controller.vis.plot_pop(self.init_tb_site)

    def plot_pops(self):
        self.fig, self.ax = self.controller.vis.plot_pops()

    def plot_heatmap(self):
        self.fig, self.ax = self.controller.vis.plot_heatmap(direction="horizontal")

    def plot_coh(self):
        self.fig, self.ax = self.controller.vis.plot_coh()

    def plot_fourier(self):
        init_state = self.controller.plot_kwargs["init_state"]
        end_state = self.controller.plot_kwargs["end_state"]
        x_axis = self.controller.plot_kwargs["x_axis"]
        self.fig, self.ax = self.controller.vis.plot_fourier(
            init_state, end_state, x_axis
        )

    # ------------------------------------------------------------------

    def save(self):
        filename = self.plotting_frame.filename_entry.get()
        directory = self.plotting_frame.directory_entry.get()
        filepath = os.path.join(directory, filename)
        save_figure(self.fig, filepath)
        self.destroy()

    def cancel(self):
        self.destroy()

    def plotting(self, fig):
        for widget in self.plotting_frame.subframe.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.plotting_frame.subframe)
        canvas.draw()
        canvas.get_tk_widget().pack()


# ----------------------------------------------------------------------


class PlottingWindow2(ctk.CTkToplevel):
    def __init__(self, master):
        """
        Notes:
            This window is a toplevel window of the main window. This means that the main window is its master.
            This window itself is the master for plotting_frame.
            The frame uses the command plot_options_kwargs, me_solver and tb_ham from the master
        """

        # initialization of the ctk.CTkToplevel class
        super().__init__(master)
        self.title("Plotting")
        self.controller = master

        self.plotting_frame = PlottingFrame(self)
        self.plotting_frame.grid(
            row=0, column=0, columnspan=2, padx=10, pady=10, sticky="nsew"
        )

        self.plot_couplings()
        self.plotting(self.fig)

    # ------------------------------------------------------------------

    def plot_couplings(self):
        self.fig, self.ax = self.controller.oligomer.plot_couplings(
            self.controller.particle,
            add_colorbar=False,
            dpi=120,
        )

    # ------------------------------------------------------------------

    def save(self):
        filename = self.plotting_frame.filename_entry.get()
        directory = self.plotting_frame.directory_entry.get()
        filepath = os.path.join(directory, filename)
        save_figure(self.fig, filepath)
        self.destroy()

    def cancel(self):
        self.destroy()

    def plotting(self, fig):
        for widget in self.plotting_frame.subframe.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.plotting_frame.subframe)
        canvas.draw()
        canvas.get_tk_widget().pack()


# ----------------------------------------------------------------------
