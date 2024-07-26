import sys
import pathlib

ROOT_DIR = str(pathlib.Path(__file__).absolute().parent.parent.parent)
if ROOT_DIR not in sys.path:
    del sys.path[0]
    sys.path.insert(0, ROOT_DIR)

import customtkinter as ctk
import webbrowser

from qDNA import TB_Ham, DNA_Seq, TB_Model, Lindblad_Diss, ME_Solver
from qDNA.tools import get_config
from qDNA.gui import (
    CustomWindow,
    InitialFrame,
    ConfigFrame,
    OptionsFrame,
    ScrollableConsoleFrame,
    PlotOptionsFrame,
    PlottingWindow,
)


class qDNA_app(ctk.CTk):
    def __init__(self):
        """
        Notes:
            This window is the main window.
            The window itself is master for initial_frame, config_frame, options_frame, plot_options_frame and scrollable_console_frame.
        """

        # initialization of the ctk.CTk class
        super().__init__()

        self.title("QuantumDNA")
        self.configs = get_config()
        self.kwargs = dict()
        self.configs = get_config()
        self.tb_basis = ["(0, 0)"]

        # frames
        self.initial_frame = InitialFrame(self)
        self.initial_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        self.config_frame = ConfigFrame(self, **self.kwargs)
        self.config_frame.grid(row=1, column=0, padx=20, pady=20, sticky="nsew")

        self.options_frame = OptionsFrame(self, self.configs, **self.kwargs)
        self.options_frame.grid(
            row=0, column=1, rowspan=2, padx=20, pady=20, sticky="nsew"
        )

        self.scrollable_console_frame = ScrollableConsoleFrame(self)
        self.scrollable_console_frame.grid(
            row=1, column=2, padx=20, pady=20, sticky="nsew"
        )

        self.plot_options_frame = PlotOptionsFrame(self, self.tb_basis, **self.kwargs)
        self.plot_options_frame.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")

        # the options_frame and plot_options_frame can not be manipulated when the window opens
        # (because at first the inputs in the initial_frame must be confirmed and is used to update these frames)
        self.options_frame.change_state("disabled")
        self.plot_options_frame.change_state("disabled")

    def open_github(self):
        webbrowser.open("https://github.com/dehe1011/QuantumDNA")

    def open_custom_window(self):
        self.custom_window = CustomWindow(self, self.configs)

    # ------------------------------------------

    def get_init_kwargs(self):
        # get the values from the initial_frame
        self.upper_strand = self.initial_frame.upper_strand_entry.get()
        self.tb_model_name = self.initial_frame.tb_model_combo.get()
        self.init_kwargs = {
            "upper_strand": self.upper_strand,
            "tb_model_name": self.tb_model_name,
        }
        self.kwargs.update(self.init_kwargs)

        # initialize dna_seq, tb_model and tb_basis
        self.dna_seq = DNA_Seq(self.upper_strand, self.tb_model_name)
        self.tb_model = TB_Model(self.dna_seq.tb_model_name, self.dna_seq.tb_dims)
        self.tb_basis = self.tb_model.tb_basis

    def get_options_kwargs(self):
        # get the values from the options_frame
        self.ham_kwargs = self.options_frame.options_tab.ham_frame.get_ham_kwargs()
        self.diss_kwargs = self.options_frame.options_tab.diss_frame.get_diss_kwargs()
        self.me_kwargs = self.options_frame.options_tab.dynamics_frame.get_me_kwargs()
        self.options_kwargs = dict(
            **self.ham_kwargs, **self.diss_kwargs, **self.me_kwargs
        )
        self.kwargs.update(self.options_kwargs)

        # initialize tb_ham, lindblad_diss, me_solver
        self.tb_ham = TB_Ham(self.dna_seq, **self.ham_kwargs)
        self.lindblad_diss = Lindblad_Diss(self.tb_ham, **self.diss_kwargs)
        self.me_solver = ME_Solver(self.tb_ham, self.lindblad_diss, **self.me_kwargs)

    def get_plot_options_kwargs(self):
        # get the values from the plot_options_frame
        self.plot_option = {
            "plot_option": self.plot_options_frame.plot_options_tab.get()
        }
        self.pop_kwargs = (
            self.plot_options_frame.plot_options_tab.pop_frame.get_pop_kwargs()
        )
        self.coh_kwargs = (
            self.plot_options_frame.plot_options_tab.coh_frame.get_coh_kwargs()
        )
        self.fourier_kwargs = (
            self.plot_options_frame.plot_options_tab.fourier_frame.get_fourier_kwargs()
        )
        self.plot_options_kwargs = dict(
            **self.plot_option,
            **self.pop_kwargs,
            **self.coh_kwargs,
            **self.fourier_kwargs
        )
        self.kwargs.update(self.plot_options_kwargs)

    # -----------------------------------------

    def press_first_confirm(self):
        """
        Event of the initial_frame.
        """

        self.get_init_kwargs()
        # update the options_frame
        self.options_frame = OptionsFrame(self, self.configs, **self.kwargs)
        self.options_frame.grid(
            row=0, column=1, rowspan=2, padx=20, pady=20, sticky="nsew"
        )
        # enable the options_frame
        self.enable_options_frame()

        # configure some widgets (since the TB basis is known)
        self.options_frame.options_tab.dynamics_frame.me_init_e_state_combo.configure(
            values=self.tb_model.tb_basis
        )
        self.options_frame.options_tab.dynamics_frame.me_init_h_state_combo.configure(
            values=self.tb_model.tb_basis
        )
        self.plot_options_frame.plot_options_tab.pop_frame.tb_site_combo.configure(
            values=self.tb_model.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.tb_site_combo.configure(
            values=self.tb_model.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.init_e_site_combo.configure(
            values=self.tb_model.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.init_h_site_combo.configure(
            values=self.tb_model.tb_basis
        )

    def press_second_confirm(self):
        """
        Event of the options_frame.
        """

        self.get_options_kwargs()
        # update the plot_options_frame
        self.plot_options_frame = PlotOptionsFrame(self, self.tb_basis, **self.kwargs)
        self.plot_options_frame.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")
        # enable the plot_options_frame
        self.enable_plotting_frame()

        # TODO: configure some widgets using the kwargs

    def submit(self):
        """
        Event of the plot_options_frame.
        """

        self.get_plot_options_kwargs()
        self.plotting_window = PlottingWindow(self)

    # ------------------------------------------------------------

    def enable_initial_frame(self):
        self.initial_frame.change_state("normal")
        self.options_frame.change_state("disabled")
        self.plot_options_frame.change_state("disabled")

    def enable_options_frame(self):
        self.initial_frame.change_state("disabled")
        self.options_frame.change_state("normal")
        self.plot_options_frame.change_state("disabled")

    def enable_plotting_frame(self):
        self.initial_frame.change_state("disabled")
        self.options_frame.change_state("disabled")
        self.plot_options_frame.change_state("normal")


if __name__ == "__main__":
    app = qDNA_app()
    app.mainloop()
