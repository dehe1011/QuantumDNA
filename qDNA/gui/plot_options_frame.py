# pylint: skip-file

import customtkinter as ctk

from .gui_utils import change_state_all_widgets

# ----------------------------------------------------------------------


class PopFrame(ctk.CTkFrame):
    def __init__(self, master, controller):
        """
        Notes:
            This frame is located in the plot_options tab. This means that the plot_options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.controller = controller

        self.explain_label = ctk.CTkLabel(
            self, text="Calculate the population of each \nbase of the DNA sequence."
        )
        self.explain_label.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

        # tb site label and combo box
        self.tb_site_label = ctk.CTkLabel(self, text="TB site:")
        self.tb_site_label.grid(row=1, column=0, padx=10, pady=10)

        self.tb_site_combo = ctk.CTkComboBox(
            self, values=["Heatmap", "All DNA Bases"] + self.controller.tb_basis
        )
        self.tb_site_combo.grid(row=1, column=1, padx=10, pady=10)
        self.tb_site_combo.set("Heatmap")

    def get_pop_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""
        pop_kwargs = {
            "init_tb_site": self.tb_site_combo.get(),
        }
        return pop_kwargs


# ----------------------------------------------------------------------


class CohFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the plot_options tab. This means that the plot_options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)

        self.explain_label = ctk.CTkLabel(
            self, text="Calculate the coherence of the DNA sequence."
        )
        self.explain_label.grid(row=0, column=0, padx=10, pady=10)

    def get_coh_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""
        coh_kwargs = {}
        return coh_kwargs


# ----------------------------------------------------------------------


class SpectrumFrame(ctk.CTkFrame):
    def __init__(self, master, controller):
        """
        Notes:
            This frame is located in the plot_options tab. This means that the plot_options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.controller = controller

        self.explain_label = ctk.CTkLabel(self, text="Calculate the eigenspectrum.")
        self.explain_label.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

        # eigenenergies checkbox
        self.eigv_var = ctk.BooleanVar(value=False)
        self.eigv_check = ctk.CTkCheckBox(
            self, text="Eigenenergies", variable=self.eigv_var
        )
        self.eigv_check.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

        # eigenstates checkbox
        self.eigs_var = ctk.BooleanVar(value=False)
        self.eigs_check = ctk.CTkCheckBox(
            self, text="Eigenstates", variable=self.eigs_var
        )
        self.eigs_check.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

        # tb site label and combo box
        self.idx_label = ctk.CTkLabel(self, text="Eigenstate index:")
        self.idx_label.grid(row=3, column=0, padx=10, pady=10)

        self.idx_combo = ctk.CTkComboBox(self, values=["0"])
        self.idx_combo.grid(row=3, column=1, padx=10, pady=10)
        self.idx_combo.set(0)

    def get_spectrum_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""
        spectrum_kwargs = {
            "eigenstate_idx": self.idx_combo.get(),
            "eigv_var": self.eigv_var.get(),
            "eigs_var": self.eigs_var.get(),
        }
        return spectrum_kwargs


# ----------------------------------------------------------------------


class FourierFrame(ctk.CTkFrame):
    def __init__(self, master, controller):
        """
        Notes:
            This frame is located in the plot_options tab. This means that the plot_options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.controller = controller

        # tb site label and combo box
        self.tb_site_label = ctk.CTkLabel(self, text="TB site:")
        self.tb_site_label.grid(row=0, column=0, padx=10, pady=10)

        self.tb_site_combo = ctk.CTkComboBox(self, values=self.controller.tb_basis)
        self.tb_site_combo.grid(row=0, column=1, padx=10, pady=10)

        # initial e site label and combo box
        self.init_e_site_label = ctk.CTkLabel(self, text="Initial electron \n TB site:")
        self.init_e_site_label.grid(row=1, column=0, padx=10, pady=10)

        self.init_e_site_combo = ctk.CTkComboBox(self, values=self.controller.tb_basis)
        self.init_e_site_combo.grid(row=1, column=1, padx=10, pady=10)

        # initial h site label and combo box
        self.init_h_site_label = ctk.CTkLabel(self, text="Initial hole \n TB site:")
        self.init_h_site_label.grid(row=2, column=0, padx=10, pady=10)

        self.init_h_site_combo = ctk.CTkComboBox(self, values=self.controller.tb_basis)
        self.init_h_site_combo.grid(row=2, column=1, padx=10, pady=10)

        # x axis label and combo box
        self.x_axis_combo = ctk.CTkComboBox(self, values=["period", "frequency"])
        self.x_axis_combo.grid(row=3, column=0, padx=10, pady=10)
        self.x_axis_combo.set("period")

        # average population checkbox
        self.average_pop_var = ctk.BooleanVar(value=False)
        self.average_pop_check = ctk.CTkCheckBox(
            self, text="Calculate average population", variable=self.average_pop_var
        )
        self.average_pop_check.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

        # --------------------------------------------------------------

    def get_fourier_kwargs(self):
        """Returns the values of widgets with get() method in dictionary format."""

        if self.controller.kwargs["description"] == "1P":
            if self.controller.kwargs["particles"] == ["electron"]:
                init_state = self.init_e_site_combo.get()
            if self.controller.kwargs["particles"] == ["hole"]:
                init_state = self.init_h_site_combo.get()
        if self.controller.kwargs["description"] == "2P":
            init_state = (self.init_e_site_combo.get(), self.init_h_site_combo.get())

        fourier_kwargs = {
            "end_state": self.tb_site_combo.get(),
            "init_state": init_state,
            "x_axis": self.x_axis_combo.get(),
            "average_pop_var": self.average_pop_var.get(),
        }
        return fourier_kwargs


# ----------------------------------------------------------------------


class ExcitonFrame(ctk.CTkFrame):
    def __init__(self, master, controller):
        """
        Notes:
            This frame is located in the plot_options tab. This means that the plot_options tab is its master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        self.pack(fill="both", expand=True)
        self.controller = controller

        # widgets
        self.lifetime_button = ctk.CTkButton(
            self,
            text="Calculate Exciton Lifetime",
            command=self.controller._calc_lifetime,
        )
        self.lifetime_button.grid(
            row=0, column=0, columnspan=2, padx=10, pady=10, sticky="ew"
        )
        self.dipole_button = ctk.CTkButton(
            self,
            text="Calculate Charge Separation",
            command=self.controller._calc_charge_separation,
        )
        self.dipole_button.grid(
            row=1, column=0, columnspan=2, padx=10, pady=10, sticky="ew"
        )
        self.dipole_moment_button = ctk.CTkButton(
            self,
            text="Calculate Dipole Moment",
            command=self.controller._calc_dipole_moment,
        )
        self.dipole_moment_button.grid(
            row=2, column=0, columnspan=2, padx=10, pady=10, sticky="ew"
        )
        self.exciton_transfer_button = ctk.CTkButton(
            self,
            text="Calculate Exciton Population",
            command=self.controller._calc_exciton_transfer,
        )
        self.exciton_transfer_button.grid(
            row=3, column=0, padx=10, pady=10, columnspan=2, sticky="ew"
        )


# ----------------------------------------------------------------------


class PlotOptionsTab(ctk.CTkTabview):
    def __init__(self, master, controller):
        """
        Notes:
            This tab is located in the plot_options frame. This means that the plot_options frame is its master.
            The tab itself is the master for four further frames.
        """

        # initialization of the ctk.CTkTabview class
        super().__init__(master)

        self.pop_tab = self.add("Population")
        self.coh_tab = self.add("Coherence")
        self.spectrum_tab = self.add("Spectrum")
        self.fourier_tab = self.add("Fourier")
        self.exciton_tab = self.add("Exciton")
        self.set("Population")

        self.pop_frame = PopFrame(self.pop_tab, controller)
        self.coh_frame = CohFrame(self.coh_tab)
        self.spectrum_frame = SpectrumFrame(self.spectrum_tab, controller)
        self.fourier_frame = FourierFrame(self.fourier_tab, controller)
        self.exciton_frame = ExcitonFrame(self.exciton_tab, controller)


class PlotOptionsFrame(ctk.CTkFrame):
    def __init__(self, master):
        """
        Notes:
            This frame is located in the main window. This means that the main window is its master.
            The frame itself is the master for the tab plot_options_tab.
            The frame uses the press_second_confirm(), enable_initial_frame() from the master.
        """

        # initialization of the ctk.CTkFrame class
        super().__init__(master)
        controller = master

        self.label = ctk.CTkLabel(
            self, text="Plotting", font=ctk.CTkFont(size=20, weight="bold")
        )
        self.label.grid(row=0, column=0, columnspan=2, pady=10, padx=10)

        self.plot_options_tab = PlotOptionsTab(self, controller)
        self.plot_options_tab.grid(row=1, column=0, columnspan=2, pady=10, padx=10)

        self.back_button = ctk.CTkButton(
            self, text="Back", command=master.enable_options_frame
        )
        self.back_button.grid(row=2, column=0, pady=10, padx=10)

        self.submit_button = ctk.CTkButton(self, text="Submit", command=master.submit)
        self.submit_button.grid(row=2, column=1, pady=10, padx=10)

    def change_state(self, state):
        """Changes the state of certain widgets (between 'normal' and 'disabled')."""
        change_state_all_widgets(self.plot_options_tab.pop_frame, state=state)
        change_state_all_widgets(self.plot_options_tab.coh_frame, state=state)
        change_state_all_widgets(self.plot_options_tab.spectrum_frame, state=state)
        change_state_all_widgets(self.plot_options_tab.fourier_frame, state=state)
        change_state_all_widgets(self.plot_options_tab.exciton_frame, state=state)
        self.submit_button.configure(state=state)
        self.back_button.configure(state=state)


# ----------------------------------------------------------------------
