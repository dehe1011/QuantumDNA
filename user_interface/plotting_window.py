import customtkinter as ctk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

from DNA import plot_pop
from utils import save_fig

class PlottingFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.filename_label = ctk.CTkLabel(self, text="Filename:")
        self.filename_label.grid(row=0, column=0, columnspan=2, pady=10, padx=10)
    
        self.filename_entry = ctk.CTkEntry(self)
        self.filename_entry.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

        self.subframe = ctk.CTkFrame(self)
        self.subframe.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

        self.save_button = ctk.CTkButton(self, text="Save", command=master.save)
        self.save_button.grid(row=11, column=0, padx=10, pady=10)
        
        self.cancel_button = ctk.CTkButton(self, text="Cancel", command=master.cancel)
        self.cancel_button.grid(row=11, column=1, padx=10, pady=10)

class PlottingWindow(ctk.CTkToplevel):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.title("Plotting")
        
        self.plotting_frame = PlottingFrame(master=self)
        self.plotting_frame.grid(row=0, column=0, columnspan=2, padx=20, pady=20, sticky="nsew")
        self.plot_pop(master)

    def plot_pop(self, master):
        self.fig, self.ax = plt.subplots()
        plot_pop(self.ax, master.plotting_kwargs["init_tb_site"], master.me_solver)
        self.plotting(self.fig)

    def save(self):
        self.filename = self.plotting_frame.filename_entry.get()
        save_fig(self.fig, self.filename, format='pdf')
        self.destroy()
    
    def cancel(self):
        self.destroy()

    def plotting(self, fig):
        for widget in self.plotting_frame.subframe.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.plotting_frame.subframe)
        canvas.draw()
        canvas.get_tk_widget().pack()
        