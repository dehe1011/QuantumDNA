# pylint: skip-file

import webbrowser

import customtkinter as ctk

from .. import DNA_Seq
from ..environment import Lindblad_Diss
from ..dynamics import ME_Solver
from ..hamiltonian import TB_Ham
from ..model import TB_Model
from ..tools import CONFIG

from .initial_frame import InitialFrame
from .config_frame import ConfigFrame, HelpFrame
from .pdb_window import PDBWindow
from .options_frame import OptionsFrame
from .plot_options_frame import PlotOptionsFrame
from .plotting_window import PlottingWindow
from .scrollable_console_frame import ScrollableConsoleFrame
from .fasta_window import FastaWindow

# --------------------------------------------------


class qDNA_app(ctk.CTk):
    def __init__(self):
        """
        Notes:
            This window is the main window.
            The window itself is master for initial_frame, config_frame, options_frame, plot_options_frame and scrollable_console_frame.
        """

        # initialization of the ctk.CTk class
        super().__init__()

        # self.tk.call('tk', 'scaling', 1.0)

        self.title("QuantumDNA")
        self.configs = CONFIG
        self.kwargs = dict()
        self.tb_basis = ["(0, 0)"]

        # Configure the grid layout for the root window
        self.grid_columnconfigure(0, weight=1)  # Column 0 takes 1 part
        self.grid_columnconfigure(1, weight=3)  # Column 1 takes 2 parts
        self.grid_columnconfigure(2, weight=1)  # Column 2 takes 1 part

        self.grid_rowconfigure(0, weight=5)  # Row 0 takes 2 parts
        self.grid_rowconfigure(1, weight=1)  # Row 1 takes 1 part
        self.grid_rowconfigure(2, weight=2)  # Row 2 takes 3 parts

        # frames
        self.initial_frame = InitialFrame(self)
        self.initial_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.config_frame = ConfigFrame(self, **self.kwargs)
        self.config_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

        self.help_frame = HelpFrame(self, **self.kwargs)
        self.help_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")

        self.options_frame = OptionsFrame(self, self.configs, **self.kwargs)
        self.options_frame.grid(
            row=0, column=1, rowspan=3, padx=10, pady=10, sticky="nsew"
        )

        self.scrollable_console_frame = ScrollableConsoleFrame(self)
        self.scrollable_console_frame.grid(
            row=1, column=2, padx=10, pady=10, rowspan=2, sticky="nsew"
        )

        self.plot_options_frame = PlotOptionsFrame(self, self.tb_basis, **self.kwargs)
        self.plot_options_frame.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")

        # the options_frame and plot_options_frame can not be manipulated when the window opens
        # (because at first the inputs in the initial_frame must be confirmed and is used to update these frames)
        self.options_frame.change_state("normal")
        self.plot_options_frame.change_state("normal")

    # ------------------------------------------

    def open_github(self):
        webbrowser.open("https://github.com/dehe1011/QuantumDNA")

    def open_documentation(self):
        webbrowser.open("https://quantumdna.readthedocs.io/en/latest/")

    def open_tutorials(self):
        webbrowser.open("https://github.com/dehe1011/QuantumDNA-notebooks")

    def open_pdb_window(self):
        self.pdb_window = PDBWindow(self, self.configs)

    def open_fasta_window(self):
        self.fasta_window = FastaWindow(self)

    # ------------------------------------------

    def get_init_kwargs(self):
        # get the values from the initial_frame
        self.upper_strand = self.initial_frame.upper_strand_entry.get()
        self.upper_strand = self.upper_strand.split("_")
        self.lower_strand = self.initial_frame.lower_strand_entry.get()
        if self.lower_strand == "auto complete":
            self.lower_strand = None
        else:
            self.lower_strand = self.lower_strand.split("_")

        self.tb_model_name = self.initial_frame.tb_model_combo.get()
        self.init_kwargs = {
            "upper_strand": self.upper_strand,
            "lower_strand": self.lower_strand,
            "tb_model_name": self.tb_model_name,
        }
        self.kwargs.update(self.init_kwargs)

        # initialize dna_seq, tb_model and tb_basis
        self.dna_seq = DNA_Seq(
            self.upper_strand, self.tb_model_name, lower_strand=self.lower_strand
        )
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
            **self.fourier_kwargs,
        )
        self.kwargs.update(self.plot_options_kwargs)

    # -----------------------------------------

    def press_first_confirm(self):
        """Event of the initial_frame."""

        self.get_init_kwargs()
        # update the options_frame
        self.options_frame = OptionsFrame(self, self.configs, **self.kwargs)
        self.options_frame.grid(
            row=0, column=1, rowspan=3, padx=10, pady=10, sticky="nsew"
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
        """Event of the options_frame."""

        self.get_options_kwargs()
        # update the plot_options_frame
        self.plot_options_frame = PlotOptionsFrame(self, self.tb_basis, **self.kwargs)
        self.plot_options_frame.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")
        # enable the plot_options_frame
        self.enable_plotting_frame()

        # TODO: configure some widgets using the kwargs

    def submit(self):
        """Event of the plot_options_frame."""

        self.get_plot_options_kwargs()
        self.plotting_window = PlottingWindow(self)

    # ------------------------------------------------------------

    def enable_initial_frame(self):
        self.initial_frame.change_state("normal")
        self.options_frame.change_state("normal")
        self.plot_options_frame.change_state("normal")

    def enable_options_frame(self):
        self.initial_frame.change_state("normal")
        self.options_frame.change_state("normal")
        self.plot_options_frame.change_state("normal")

    def enable_plotting_frame(self):
        self.initial_frame.change_state("normal")
        self.options_frame.change_state("normal")
        self.plot_options_frame.change_state("normal")
