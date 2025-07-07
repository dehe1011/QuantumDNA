import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ..evaluation import Evaluation
from ..hamiltonian import get_pop_fourier
from ..dynamics import get_reduced_dm_eigs
from ..utils import get_conversion
from ..io import OPTIONS, load_color_palette

# ----------------------------------------------------------------------


def _get_colors():
    dna_bases = OPTIONS["dna_bases"]
    particles = OPTIONS["particles"]

    # A, T, G, C, F
    color_palette = load_color_palette("seaborn")["icefire7"]
    colors_dna_bases = dict(zip(dna_bases, [color_palette[i] for i in [5, 6, 1, 0, 3]]))

    # electron, hole, exciton
    # color_palette = load_color_palette("seaborn")["icefire5"]
    # colors_particles = dict(zip(particles, [color_palette[i] for i in [0, 4, 2]]))
    colors_particles = dict(zip(particles, ["#459DD9", "#B74244", "#4EB572"]))

    return colors_dna_bases, colors_particles


COLORS_DNA_BASES, COLORS_PARTICLES = _get_colors()

# ----------------------------------------------------------------------


class Visualization(Evaluation):
    def __init__(self, tb_sites, **kwargs):

        self.kwargs = kwargs
        self.tb_sites = tb_sites
        super().__init__(self.tb_sites, **kwargs)

    # ----------------------------------------------------------------------

    def plot_heatmap(
        self, heatmap_type="seaborn", fig=None, ax=None, dpi=None, **plot_kwargs
    ):
        """
        Plots a heatmap of populations for each particle in the system.

        Parameters
        ----------
        me_solver : object
            An instance of a solver that contains the time evolution data and Hamiltonian information.
        vmax_list : list of float, optional
            A list of maximum values for the color scale of each particle's heatmap. Default is [1, 1, 1].
        heatmap_type : str, optional
            The type of heatmap to plot. Options are "seaborn" or "matplotlib". Default is "seaborn".

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the heatmaps.
        ax : numpy.ndarray of matplotlib.axes._subplots.AxesSubplot
            An array of axes objects for each particle's heatmap.

        Notes
        -----
        .. note::

            - The function uses seaborn or matplotlib to generate heatmaps.
            - The color maps used are "Blues", "Reds", and "Greens" for different particles.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        direction = plot_kwargs.get("direction", "horizontal")
        if "direction" in plot_kwargs:
            del plot_kwargs["direction"]

        num_particles = len(self.particles)
        if direction == "vertical":
            x_num, y_num = num_particles, 1
        else:
            x_num, y_num = 1, num_particles

        if fig is None:
            if x_num == 1 and y_num == 1:
                fig, ax = plt.subplots(dpi=dpi)
            else:
                fig, ax = plt.subplots(
                    x_num,
                    y_num,
                    figsize=(3.4 * y_num, 2.1 * x_num),
                    sharex=True,
                    sharey=True,
                    dpi=dpi,
                )

        ax = ax.reshape((x_num, y_num))

        # ---

        pop_dict = self.get_pop()
        cmaps = {"electron": "Blues", "hole": "Reds", "exciton": "Greys"}
        cmaps = {"electron": "Blues", "hole": "Reds", "exciton": "Greens"}

        for i in range(x_num):
            for j in range(y_num):
                particle = self.particles[i + j]
                particle_pop = np.array(
                    [
                        value
                        for key, value in pop_dict.items()
                        if key.startswith(particle)
                    ]
                )
                vmax = 1
                if particle == "exciton":
                    vmax = np.max(particle_pop)

                # seaborn heatmap (looks prettier in my opinion)
                if heatmap_type == "seaborn":
                    heatmap = sns.heatmap(
                        particle_pop,
                        xticklabels=[],
                        yticklabels=[],
                        cmap=cmaps[particle],
                        ax=ax[i, j],
                        cbar=False,
                        vmax=vmax,
                        **plot_kwargs,
                    )
                    heatmap.figure.colorbar(heatmap.collections[0], ax=ax[i, j])

                # matplotlib heatmap
                if heatmap_type == "matplotlib":
                    im = ax[i, j].imshow(
                        particle_pop,
                        cmap=cmaps[particle],
                        aspect="auto",
                        vmax=vmax,
                        **plot_kwargs,
                    )
                    im.figure.colorbar(im, ax=ax[i])

                ax[i, j].set_ylabel(particle.capitalize())

                # ticks
                ax[i, j].set_xticks([])
                ax[i, j].set_yticks([])
                y_len, x_len = particle_pop.shape
                xticks = np.linspace(0, x_len, 7)
                ax[i, j].set_xticks(
                    xticks, labels=np.round(np.linspace(0, self.t_end, 7), 0)
                )
                yticks = np.arange(y_len) + 0.5
                ax[i, j].set_yticks(yticks, labels=self.tb_sites_flattened)

        for j in range(y_num):
            ax[-1, j].set_xlabel("Time [" + self.t_unit + "]")
        return fig, ax

    def plot_pop(
        self, tb_site, fig=None, ax=None, dpi=None, add_legend=True, **plot_kwargs
    ):
        """Plots population of one base.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_site : str
            The tight-binding site.
        me_solver : ME_Solver
            The master equation solver instance.
        add_legend : bool, optional
            Whether to add a legend, by default True.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        pop = self.get_pop()

        # plotting
        for particle in self.particles:

            plot_kwargs["color"] = COLORS_PARTICLES[particle]
            plot_kwargs["label"] = particle

            ax.plot(
                self.times,
                pop[particle + "_" + tb_site],
                **plot_kwargs,
            )

        dna_base = self.tb_basis_sites_dict[tb_site]
        x_center = self.t_end / 2
        y_center = 0.8
        ax.text(
            x_center,
            y_center,
            dna_base,
            ha="center",
            va="center",
            color="grey",
            fontsize=15,
            # fontweight='bold',
            bbox={
                "facecolor": "white",
                "edgecolor": "white",
                "boxstyle": "round,pad=0.2",
            },
        )

        # plot settings
        ax.set_ylim(0, 1.02)
        if add_legend:
            ax.set_ylabel("Population")
            ax.set_xlabel("Time [" + self.t_unit + "]")
            ax.legend()
        return fig, ax

    def plot_pops(self, fig=None, ax=None, dpi=None, **plot_kwargs):
        """Plots populations of all bases.

        Parameters
        ----------
        me_solver : ME_Solver
            The master equation solver instance.

        Returns
        -------
        tuple
            A tuple containing the figure and axes.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        direction = plot_kwargs.get("direction", "horizontal")
        if "direction" in plot_kwargs:
            del plot_kwargs["direction"]

        if direction == "vertical":
            x_num, y_num = self.num_sites_per_strand, self.num_channels
        else:
            x_num, y_num = self.num_channels, self.num_sites_per_strand

        if fig is None:
            if x_num == 1 and y_num == 1:
                fig, ax = plt.subplots(dpi=dpi)
            else:
                fig, ax = plt.subplots(
                    x_num,
                    y_num,
                    figsize=(3.4 * y_num, 2.1 * x_num),
                    sharex=True,
                    sharey=True,
                    dpi=dpi,
                )

        ax = ax.reshape((x_num, y_num))
        for i in range(x_num):
            ax[i, 0].set_ylabel("Population")
        for j in range(y_num):
            ax[-1, j].set_xlabel("Time [" + self.t_unit + "]")

        # ---

        for i in range(x_num):
            for j in range(y_num):
                if direction == "vertical":
                    tb_site = f"({j}, {i})"
                else:
                    tb_site = f"({i}, {j})"
                _, ax[i, j] = self.plot_pop(
                    tb_site, fig, ax[i, j], dpi, add_legend=False
                )

        ax[0, 0].legend(self.particles, loc="upper right")

        return fig, ax

    def plot_pop_fourier(
        self,
        init_state,
        end_state,
        times,
        t_unit,
        fig=None,
        ax=None,
        dpi=None,
        add_legend=True,
        **plot_kwargs,
    ):
        """Plots population using the Fourier decomposition.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        init_state : str
            The initial state.
        end_state : str
            The end state.
        times : list of float
            The time points for plotting.
        t_unit : str
            The unit of time.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        # calculation
        self.unit = "rad/" + t_unit
        amplitudes_dict, frequencies_dict, average_pop_dict = self.get_fourier(
            init_state, end_state, ["amplitude", "frequency", "average_pop"]
        )
        pop_dict = {}
        for particle in self.particles:
            amplitudes = amplitudes_dict[particle]
            frequencies = frequencies_dict[particle]
            average_pop = average_pop_dict[particle]
            pop_dict[particle] = [
                get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times
            ]

        # plotting
        for particle in self.particles:
            ax.plot(
                times,
                pop_dict[particle],
                label=particle,
                color=COLORS_PARTICLES[particle],
                **plot_kwargs,
            )

        # plot settings
        ax.set_ylim(0, 1.02)
        if add_legend:
            ax.set_ylabel("Population")
            ax.set_xlabel("Time [" + self.t_unit + "]")
            ax.legend(self.particles)
        return fig, ax

    def plot_coh(self, fig=None, ax=None, dpi=None, **plot_kwargs):
        """Plots a measure of coherence as the sum of absolute values of the off- diagonal
        elements.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        me_solver : ME_Solver
            The master equation solver instance.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        # calculation
        coh = self.get_coh()
        if plot_kwargs is None:
            plot_kwargs = {}

        # plotting
        for particle in self.particles:

            plot_kwargs["color"] = COLORS_PARTICLES[particle]
            plot_kwargs["label"] = particle

            ax.plot(
                self.times,
                coh[particle],
                **plot_kwargs,
            )

        # plot settings
        ax.set_ylabel("Coherence")
        ax.set_xlabel("Time [" + self.t_unit + "]")
        ax.legend()
        return fig, ax

    def plot_test_fourier(self, tb_site, fig=None, ax=None, dpi=None, **plot_kwargs):
        """Plots a comparison of the result obtained by Fourier decomposition and the result
        obtained solving the Schrödinger equation with the qutip mesolver.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_site : str
            The tight-binding site.
        me_solver : ME_Solver
            The master equation solver instance.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        self.plot_pop_fourier(
            self.init_states[0],
            tb_site,
            self.times,
            self.t_unit,
            fig,
            ax,
            dpi,
            **plot_kwargs,
        )
        self.plot_pop(tb_site, fig, ax, dpi, **plot_kwargs)
        return fig, ax

    # ----------------------------------------------------------------------

    def plot_eigv(self, energy_unit="eV", fig=None, ax=None, dpi=None, color=None):
        """Plots the eigenenergies.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        energy_unit : str, optional
            The unit of energy to plot, by default '100meV'.
        """

        if fig is None:
            fig, ax = plt.subplots(figsize=(3.4, 3.4), dpi=dpi)

        # calculation
        eigv, _ = self.get_eigensystem()
        eigv *= get_conversion(self.unit, energy_unit)

        # plotting
        x_start, x_end = 0, 1

        if color is None:
            color = "black"
        for e in eigv:
            ax.hlines(y=e, xmin=x_start, xmax=x_end, color=color)

        # Optional: Layout anpassen
        ax.set_xlim(x_start, x_end)
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.grid(True, axis="y", linestyle=":", alpha=0.4)
        ax.set_ylabel("Energy in " + energy_unit)
        return fig, ax

    def plot_eigs(self, eigenstate_idx, fig=None, ax=None, dpi=None):
        """Plots the distribution of eigenstates over tight-binding sites.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        eigenstate_idx : int
            The index of the eigenstate to plot.
        """

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        # calculation
        _, eigs = self.get_eigensystem()

        dm = None
        for particle in self.particles:
            if self.description == "2P":
                dm = get_reduced_dm_eigs(eigs, eigenstate_idx, self.tb_basis, particle)

            elif self.description == "1P":
                dm = np.outer(eigs[:, eigenstate_idx], eigs[:, eigenstate_idx].conj())

            eigs_distribution = np.diag(dm).real
            if particle != "exciton":
                assert np.allclose(
                    sum(eigs_distribution), 1, atol=1e-2
                ), "The distribution does not sum to 1."

            ax.plot(
                range(self.num_sites),
                eigs_distribution,
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
            ax.set_xticks(range(self.num_sites))
            ax.set_xticklabels(self.tb_sites_flattened)
        ax.set_ylim(0, 1.02)
        ax.set_title(f"Distribution of Eigenstate {eigenstate_idx}")
        ax.legend()
        return fig, ax

    def plot_fourier(self, init_state, end_state, x_axis, fig=None, ax=None, dpi=None):
        """Plots the Fourier transform of the transition amplitudes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        init_state : Any
            The initial state.
        end_state : Any
            The end state.
        x_axis : str
            The x-axis type ('frequency' or 'period').
        """

        if fig is None:
            fig, ax = plt.subplots(dpi=dpi)

        # calculation
        amplitudes_dict = self.get_amplitudes(init_state, end_state)
        frequencies_dict = self.get_frequencies(init_state, end_state)
        # transform frequencies to rad/ps

        markers = {"electron": "^", "hole": "v", "exciton": "*"}
        for particle in self.particles:
            conversion = get_conversion(self.unit, "rad/ps") / (2 * np.pi)
            frequencies_dict[particle] = (
                np.array(frequencies_dict[particle]) * conversion
            )

            amplitudes = amplitudes_dict[particle]
            frequencies = frequencies_dict[particle]
            # transform frequency to period (in fs)
            periods = 1e3 / frequencies

            # frequency as x-axis
            if x_axis.lower() == "frequency":
                ax.plot(
                    frequencies,
                    np.abs(amplitudes),
                    ls="",
                    marker=markers[particle],
                    label=particle,
                    color=COLORS_PARTICLES[particle],
                    markersize=10,
                    alpha=0.8,
                )

            # period as x-axis
            elif x_axis.lower() == "period":
                ax.plot(
                    periods,
                    amplitudes,
                    ls="",
                    marker=markers[particle],
                    label=particle,
                    color=COLORS_PARTICLES[particle],
                    markersize=10,
                    alpha=0.8,
                )

        # plot settings
        if x_axis.lower() == "frequency":
            ax.set_xlabel("Frequency in rad/ps")
        elif x_axis.lower() == "period":
            ax.set_xlabel("Period in fs")
        ax.set_ylabel("Amplitude")
        # ax.set_ylim(0.02)
        ax.legend()
        return fig, ax

    # ----------------------------------------------------------------------

    def get_cumulative_average_pop(self, J_list, J_unit):
        """Computes the cumulative average population.

        Parameters
        ----------
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        J_list : list
            List of interaction parameters.

        Returns
        -------
        numpy.ndarray
            Cumulative average population.
        """

        # pop_list contains the average population for each particle, J, and tb_site
        pop_list = np.zeros((len(self.particles), len(J_list), self.num_sites))

        # calculate the average population for each particle, J, and tb_site using tb_ham.get_average_pop
        self.unit = J_unit
        for J_idx, J in enumerate(J_list):
            self.coulomb_param = J
            for tb_site_idx, tb_site in enumerate(self.tb_basis):
                average_pop = self.get_average_pop(self.init_states[0], tb_site)
                for particle_idx, particle in enumerate(self.particles):
                    pop_list[particle_idx][J_idx][tb_site_idx] = average_pop[particle]

        # calculate the cumulative average population
        cumulative_pop_list = [0] * (self.num_sites + 1)

        # add zero population
        running_pop_list = np.zeros((len(self.particles), len(J_list)))
        cumulative_pop_list[-1] = np.array(running_pop_list)

        # add cumulative population
        for tb_basis_idx in range(self.num_sites):
            running_pop_list += pop_list[:, :, tb_basis_idx]
            cumulative_pop_list[tb_basis_idx] = np.array(running_pop_list)

        return np.array(cumulative_pop_list)

    # cumulative_average_pop = get_cumulative_average_pop(tb_ham, J_list)
    def plot_average_pop(
        self, J_list, J_unit="100meV", fig=None, ax=None, dpi=None, **plot_kwargs
    ):
        """Plots the average population for each interaction parameter.

        Parameters
        ----------
        ax : list of matplotlib.axes.Axes
            The axes to plot on.
        tb_ham : TB_Ham
            The tight-binding Hamiltonian.
        J_list : list
            List of interaction parameters.
        J_unit : str
            The unit of the interaction parameters.
        """

        if plot_kwargs is None:
            plot_kwargs = {}

        direction = plot_kwargs.get("direction", "horizontal")
        if "direction" in plot_kwargs:
            del plot_kwargs["direction"]

        num_particles = len(self.particles)
        if direction == "vertical":
            x_num, y_num = num_particles, 1
        else:
            x_num, y_num = 1, num_particles

        if fig is None:
            if x_num == 1 and y_num == 1:
                fig, ax = plt.subplots(dpi=dpi)
            else:
                fig, ax = plt.subplots(
                    x_num,
                    y_num,
                    figsize=(3.4 * y_num, 3.4 * x_num),
                    sharex=True,
                    sharey=True,
                    dpi=dpi,
                )

        ax = ax.reshape((x_num, y_num))
        # ---

        dna_seq = self.tb_sites_flattened
        cumulative_average_pop = self.get_cumulative_average_pop(J_list, J_unit)

        for i in range(x_num):
            for j in range(y_num):
                for dna_base_idx, dna_base in enumerate(dna_seq):
                    # black lines
                    ax[i, j].plot(
                        J_list,
                        cumulative_average_pop[dna_base_idx][:][i + j],
                        color="k",
                        lw=1.5,
                    )

                    # fill between the black lines
                    ax[i, j].fill_between(
                        J_list,
                        cumulative_average_pop[dna_base_idx - 1][:][i + j],
                        cumulative_average_pop[dna_base_idx][:][i + j],
                        color=COLORS_DNA_BASES[dna_base],
                        alpha=0.3,
                    )
                    particle = self.particles[i + j]
                    ax[i, j].set_title(particle.capitalize(), fontsize=15)

        # plot settings
        for i in range(x_num):
            ax[i, 0].set_ylabel("Acc. Population")
        for j in range(y_num):
            # plot the bottom line
            ax[-1, j].plot(J_list, [0] * len(J_list), color="k", lw=1.5)
            ax[-1, j].set_xlabel("J [" + self.unit + "]")
            # ax[particle_idx].set_ylim(0, 1.05)

        # plot the DNA bases as letters
        for dna_base_idx, dna_base in enumerate(dna_seq):
            ax[0, 0].text(
                0,
                dna_base_idx / len(dna_seq),
                dna_base,
                fontsize=15,
                color="k",
                alpha=0.8,
            )


# ----------------------------------------------------------------------
