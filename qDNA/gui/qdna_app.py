# pylint: skip-file

import webbrowser

import customtkinter as ctk

from ..model import TBModel
from ..hamiltonian import get_tb_sites
from ..visualization import Visualization
from ..evaluation import Evaluation
from ..io import OPTIONS, DEFAULTS

from .initial_frame import InitialFrame
from .advanced_frame import AdvancedFrame, HelpFrame
from .options_frame import OptionsFrame
from .plot_options_frame import PlotOptionsFrame
from .scrollable_console_frame import ScrollableConsoleFrame
from .plotting_window import PlottingWindow
from .pdb_window import PDBWindow
from .fasta_window import FastaWindow

# --------------------------------------------------


class QDNApp(ctk.CTk):
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
        self.kwargs = DEFAULTS["tb_model_kwargs_default"]
        self.tb_basis = ["(0, 0)"]
        self.options = {**OPTIONS, **DEFAULTS}

        # Configure the grid layout for the root window
        self.grid_columnconfigure(0, weight=1)  # Column 0 takes 1 part
        self.grid_columnconfigure(1, weight=3)  # Column 1 takes 2 parts
        self.grid_columnconfigure(2, weight=1)  # Column 2 takes 1 part

        self.grid_rowconfigure(0, weight=5)  # Row 0 takes 2 parts
        self.grid_rowconfigure(1, weight=1)  # Row 1 takes 1 part
        self.grid_rowconfigure(2, weight=2)  # Row 2 takes 3 parts

        # --------------------------------------------------------------

        # left frames
        self.initial_frame = InitialFrame(self)
        self.initial_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.advanced_frame = AdvancedFrame(self)
        self.advanced_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

        self.help_frame = HelpFrame(self)
        self.help_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")

        # --------------------------------------------------------------

        # middle frames
        self.options_frame = OptionsFrame(self)
        self.options_frame.grid(
            row=0, column=1, rowspan=3, padx=10, pady=10, sticky="nsew"
        )

        # --------------------------------------------------------------

        # right frames
        self.scrollable_console_frame = ScrollableConsoleFrame(self)
        self.scrollable_console_frame.grid(
            row=1, column=2, padx=10, pady=10, rowspan=2, sticky="nsew"
        )

        self.plot_options_frame = PlotOptionsFrame(self)
        self.plot_options_frame.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")

        # the options_frame and plot_options_frame can not be manipulated when the window opens
        # (because at first the inputs in the initial_frame must be confirmed and is used to update these frames)
        self.options_frame.change_state("disabled")
        self.plot_options_frame.change_state("disabled")

    # ------------------------------------------------------------------

    def open_github(self):
        webbrowser.open("https://github.com/dehe1011/QuantumDNA")

    def open_documentation(self):
        webbrowser.open("https://quantumdna.readthedocs.io/en/latest/")

    def open_tutorials(self):
        webbrowser.open("https://github.com/dehe1011/QuantumDNA-notebooks")

    def open_pdb_window(self):
        self.pdb_window = PDBWindow(self)

    def open_fasta_window(self):
        self.fasta_window = FastaWindow(self)

    # ------------------------------------------------------------------

    def get_init_kwargs(self):
        # get the values from the initial_frame

        self.upper_strand_input = self.initial_frame.upper_strand_entry.get()
        self.upper_strand = self.upper_strand_input.split("/")
        self.lower_strand_input = self.initial_frame.lower_strand_entry.get()
        if self.lower_strand_input != "auto complete":
            self.lower_strand = self.lower_strand_input.split("/")
        else:
            self.lower_strand = self.lower_strand_input
        self.tb_model_name = self.initial_frame.tb_model_combo.get()

        self.kwargs["tb_model_name"] = self.tb_model_name
        self.tb_model = TBModel(len(self.upper_strand), **self.kwargs)
        self.tb_sites = get_tb_sites(
            self.upper_strand,
            lower_strand=self.lower_strand,
            tb_model_name=self.tb_model_name,
        )
        self.tb_basis = self.tb_model.tb_basis

        if len(self.upper_strand_input) >= 8:
            print(
                "Info: This is a long sequence. The calculation may take some time."
                + "\n-------------------------------"
            )
        if self.tb_model_name in ["FWM", "FLM", "FELM", "TC", "FC"]:
            print(
                "Info: Many predefined TB parametrizations (sources) are not available for this model."
                + "\n-------------------------------"
            )

    # ------------------------------------------------------------------

    def get_options_kwargs(self):
        # get the values from the options_frame
        self.tb_ham_kwargs = (
            self.options_frame.options_tab.tb_ham_frame.get_tb_ham_kwargs()
        )
        self.lind_diss_kwargs = (
            self.options_frame.options_tab.lind_diss_frame.get_lind_diss_kwargs()
        )
        self.me_solver_kwargs = (
            self.options_frame.options_tab.me_solver_frame.get_me_solver_kwargs()
        )
        self.options_kwargs = dict(
            **self.tb_ham_kwargs, **self.lind_diss_kwargs, **self.me_solver_kwargs
        )
        self.kwargs.update(self.options_kwargs)

    # ------------------------------------------------------------------

    def get_plot_kwargs(self):
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
        self.spectrum_kwargs = (
            self.plot_options_frame.plot_options_tab.spectrum_frame.get_spectrum_kwargs()
        )
        self.fourier_kwargs = (
            self.plot_options_frame.plot_options_tab.fourier_frame.get_fourier_kwargs()
        )
        self.plot_kwargs = dict(
            **self.plot_option,
            **self.pop_kwargs,
            **self.coh_kwargs,
            **self.spectrum_kwargs,
            **self.fourier_kwargs,
        )

    # ------------------------------------------------------------------

    def press_first_confirm(self):
        """Event of the initial_frame."""

        self.get_init_kwargs()

        # update the options_frame
        # self.options_frame = OptionsFrame(self, self.tb_basis)
        # self.options_frame.grid(
        #     row=0, column=1, rowspan=3, padx=10, pady=10, sticky="nsew"
        # )
        # enable the options_frame
        self.enable_options_frame()

        # configure some widgets (since the TB basis is known)
        self.options_frame.options_tab.me_solver_frame.init_e_state_combo.configure(
            values=self.tb_basis
        )
        self.options_frame.options_tab.me_solver_frame.init_h_state_combo.configure(
            values=self.tb_basis
        )
        self.plot_options_frame.plot_options_tab.pop_frame.tb_site_combo.configure(
            values=self.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.tb_site_combo.configure(
            values=self.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.init_e_site_combo.configure(
            values=self.tb_basis
        )
        self.plot_options_frame.plot_options_tab.fourier_frame.init_h_site_combo.configure(
            values=self.tb_basis
        )

    def press_second_confirm(self):
        """Event of the options_frame."""

        self.get_options_kwargs()

        self.eva = Evaluation(self.tb_sites, **self.kwargs)

        if self.eva.description == "2P":
            dim = self.eva.num_sites**2
        else:
            dim = self.eva.num_sites
        self.plot_options_frame.plot_options_tab.spectrum_frame.idx_combo.configure(
            values=[str(i) for i in range(dim)]
        )

        # update the plot_options_frame
        # self.plot_options_frame = PlotOptionsFrame(self, self.tb_basis)
        # self.plot_options_frame.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")
        # enable the plot_options_frame
        self.enable_plotting_frame()

    def submit(self):
        """Event of the plot_options_frame."""

        self.get_plot_kwargs()
        self.vis = Visualization(self.tb_sites, **self.kwargs)
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

    # ----------------------------------------------------------------------

    def _calc_lifetime(self):
        assert (
            self.kwargs["description"] == "2P"
        ), "2P description is required for the calculation."
        assert (
            self.kwargs["relaxation"] == True
        ), "Groundstate is required for the calculation."

        lifetime = self.eva.calc_lifetime()
        if isinstance(lifetime, str):
            print(f"Exciton Lifetime: {lifetime}" "\n-------------------------------")
        else:
            print(
                f"Exciton Lifetime: {lifetime} fs" "\n-------------------------------"
            )

    def _calc_charge_separation(self):
        assert (
            self.kwargs["description"] == "2P"
        ), "2P description is required for the calculation."
        assert (
            self.kwargs["relaxation"] == True
        ), "Groundstate is required for the calculation."

        charge_separation = self.eva.calc_charge_separation()
        print(
            f"Charge Separation: {charge_separation} A"
            "\n-------------------------------"
        )

    def _calc_dipole_moment(self):
        assert (
            self.kwargs["description"] == "2P"
        ), "2P description is required for the calculation."
        assert (
            self.kwargs["relaxation"] == True
        ), "Groundstate is required for the calculation."

        dipole_moment = self.eva.calc_dipole_moment()
        print(f"Dipole Moment: {dipole_moment} D" "\n-------------------------------")

    def _calc_exciton_transfer(self):
        assert (
            self.kwargs["description"] == "2P"
        ), "2P description is required for the calculation."
        assert (
            self.kwargs["relaxation"] == True
        ), "Groundstate is required for the calculation."
        assert (
            "exciton" in self.kwargs["particles"]
        ), "Exciton must be selected in particles."

        avg_pop_upper, avg_pop_lower = self.eva.calc_exciton_transfer().values()
        avg_pop_upper, avg_pop_lower = (
            avg_pop_upper["exciton"],
            avg_pop_lower["exciton"],
        )
        print(
            f"Average Exciton Population (upper strand): {avg_pop_upper}"
            + f"\nAverage Exciton Population (lower strand): {avg_pop_lower}"
            + "\n-------------------------------"
        )

    # ------------------------------------------------------------------


# alias for legacy code compatibility
class qDNA_app(QDNApp):
    pass


# ----------------------------------------------------------------------

if __name__ == "__main__":
    app = QDNApp()
    app.geometry("1200x800")
    app.resizable(True, True)
    app.mainloop()

# ----------------------------------------------------------------------
