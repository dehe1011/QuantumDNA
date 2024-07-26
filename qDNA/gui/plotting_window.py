import customtkinter as ctk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from qDNA import plot_pop, plot_pops, plot_coh, plot_fourier
from qDNA.tools import save_figure


class PlottingFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        """
        Notes:
            This frame is located in the plotting_window. This means that the plotting_window is its master.
            The frame uses the commands save() and cancel() from the master

        Widgets with get() method:
            filename_entry
            directory_entry
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master, **kwargs)

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
    def __init__(self, master, **kwargs):
        """
        Notes:
            This window is a toplevel window of the main window. This means that the main window is its master.
            This window itself is the master for plotting_frame.
            The frame uses the command plot_options_kwargs, me_solver and tb_ham from the master
        """

        # initialization of the ctk.CTkToplevel class
        super().__init__(master, **kwargs)
        self.title("Plotting")

        self.plotting_frame = PlottingFrame(master=self)
        self.plotting_frame.grid(
            row=0, column=0, columnspan=2, padx=20, pady=20, sticky="nsew"
        )

        self.plot_options_kwargs = master.plot_options_kwargs
        self.plot_option = self.plot_options_kwargs["plot_option"]
        self.me_solver = master.me_solver
        self.tb_ham = master.tb_ham

        if self.plot_option == "Population":
            self.init_tb_site = self.plot_options_kwargs["init_tb_site"]
            if self.init_tb_site == "All":
                self.plot_pops()
            else:
                self.plot_pop()

        if self.plot_option == "Coherence":
            self.plot_coh()

        if self.plot_option == "Fourier":
            self.plot_fourier()
            if (
                master.plot_options_frame.plot_options_tab.fourier_frame.average_pop_var.get()
            ):
                self.average_pop()
        self.plotting(self.fig)

    def average_pop(self):
        average_pop = self.tb_ham.get_average_pop(
            self.plot_options_kwargs["init_state"],
            self.plot_options_kwargs["end_state"],
        )
        print(f"Average population: {average_pop}")

    def plot_pop(self):
        self.fig, self.ax = plt.subplots()
        plot_pop(self.ax, self.init_tb_site, self.me_solver)

    def plot_pops(self):
        self.fig, self.axes = plot_pops(self.me_solver)

    def plot_coh(self):
        self.fig, self.ax = plt.subplots()
        plot_coh(self.ax, self.me_solver)

    def plot_fourier(self):
        self.fig, self.ax = plt.subplots()
        init_state = self.plot_options_kwargs["init_state"]
        end_state = self.plot_options_kwargs["end_state"]
        x_axis = self.plot_options_kwargs["x_axis"]
        plot_fourier(self.ax, self.tb_ham, init_state, end_state, x_axis)

    def save(self):
        self.filename = self.plotting_frame.filename_entry.get()
        self.directory = self.plotting_frame.directory_entry.get()
        if self.directory:
            save_figure(self.fig, self.filename, directory=self.directory, format="pdf")
        else:
            save_figure(self.fig, self.filename, format="pdf")
        self.destroy()

    def cancel(self):
        self.destroy()

    def plotting(self, fig):
        for widget in self.plotting_frame.subframe.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.plotting_frame.subframe)
        canvas.draw()
        canvas.get_tk_widget().pack()
