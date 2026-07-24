import os
import shutil
import copy

import numpy as np
import matplotlib.pyplot as plt

from ..utils import check_lcao_kwargs
from ..io import load_lcao_param, pdb_to_xyz, find_xyz, save_tb_params, DEFAULTS
from ..model import (
    get_tb_couplings,
    get_tb_energies,
    TBModel,
    str_to_tuple,
    tuple_to_int,
)
from .slater_koster import calc_H_inter
from .dipolar import calc_dipolar_coupling
from .monomer import Monomer

# ----------------------------------------------------------------------


class Oligomer(TBModel):
    """
    Oligomer class for calculating TB parameters of DNA oligomers.
    This class extends the TBModel class and provides functionality to calculate
    tight-binding parameters from PDB files. It processes the PDB file to generate
    XYZ files (components) and groups the components into monomers aranged in a geometry
    determined by the required TB model.

    Attributes
    ----------
    filepath_pdb : str
        Path to the input PDB file.
    directory : str
        Directory where XYZ files are stored.
    kwargs : dict
        Configuration parameters for the model.
    filename_pdb : str
        Base name of the PDB file without extension.
    sites_bases : list
        List of base site filenames.
    sites_backbone : list
        List of backbone site filenames.
    sites : list
        List of site groups based on the number of channels.
    sites_id : numpy.ndarray
        Array of site identifiers reshaped to match tight-binding dimensions.
    filepaths : list
        File paths to XYZ files for each site.
    num_sites : int
        Total number of sites in the oligomer.
    monomers : numpy.ndarray
        Array of Monomer objects representing individual sites.
    couplings : list
        Tight-binding coupling information.
    energies : list
        Tight-binding energy information.
    tb_params : dict or None
        Calculated tight-binding parameters, including hole, electron, and exciton
        properties.

    Methods
    -------
    calc_tb_couplings(monomer1, monomer2)
        Calculate tight-binding couplings between two monomers.
    calc_tb_energies(monomer)
        Calculate tight-binding energies for a monomer.
    calc_tb_params()
        Compute and return tight-binding parameters for the oligomer.
    save_tb_params(directory=None)
        Save tight-binding parameters to a file.
    plot_couplings(particle, add_colorbar=False, add_label=True, fig=None, ax=None, dpi=None, max_coupling=None)
        Plot tight-binding couplings for a specified particle type.
    clean()
        Clean up temporary files and directories created during processing.

    Notes
    -----
    .. note::
        - The class assumes a specific structure for the input PDB file and generates XYZ files accordingly.
        - Tight-binding parameters are calculated based on the provided LCAO model and configuration.

    """

    def __init__(self, filepath_pdb, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(DEFAULTS["lcao_kwargs_default"])
        self.kwargs.update(kwargs)
        # check_lcao_kwargs(**self.kwargs) # TODO

        # load LCAO parameters
        self.lcao_param = load_lcao_param(self.kwargs.get("param_id"))  # self.param_id # TODO

        # prepare xyz files
        self.filepath_pdb = filepath_pdb
        self.filename_pdb = os.path.splitext(os.path.basename(self.filepath_pdb))[0]
        self.directory = self.filepath_pdb.split(".")[0]
        if os.path.exists(self.directory):
            self.clean()
        pdb_to_xyz(self.filepath_pdb, **self.kwargs)
        self.filenames_xyz = find_xyz(self.directory)

        # init TB model
        num_sites_per_strand = len(self.filenames_xyz) // 4
        super().__init__(num_sites_per_strand, **self.kwargs)
        self.couplings = get_tb_couplings(self.tb_model_name, self.num_sites_per_strand)
        self.energies = get_tb_energies(self.num_channels, self.num_sites_per_strand)

        # init monomers
        self.sites_all = self._get_sites_all()
        self.sites = self._get_sites()
        self.sites_id = self._get_sites_id()
        self.monomers = self._get_monomers()

        self.tb_params = None

        if self.kwargs.get("auto_clean"):
            self.clean()

    # ------------------------------------------------------------------

    def __repr__(self):
        return f"Oligomer({self.filename_pdb})"

    def _get_sites_all(self):
        """Orders XYZ filenames (components) in four channels."""

        sites_bases, sites_backbone = [], []
        for filename in self.filenames_xyz:
            if "B" in filename:
                sites_backbone.append(filename)
            else:
                sites_bases.append(filename)

        n = self.num_sites_per_strand
        sites_all = np.array(
            [
                sites_backbone[:n],
                sites_bases[:n],
                sites_bases[n:][::-1],
                sites_backbone[n:][::-1],
            ],
        )

        return sites_all

    def _get_sites(self):
        """Adapts to the TB model, i.e., groups the monomers accordingly."""

        sites = []
        n = self.num_sites_per_strand
        backbone = self.tb_model_props["backbone"]

        if self.num_channels == 1:
            if backbone:
                sites = [self.sites_all[0:4, i].tolist() for i in range(n)]
            else:
                sites = [self.sites_all[1:3, i].tolist() for i in range(n)]

        if self.num_channels == 2:
            if backbone:
                sites_upper = [self.sites_all[0:2, i].tolist() for i in range(n)]
                sites_lower = [self.sites_all[2:4, i].tolist() for i in range(n)]
            else:
                sites_upper = [self.sites_all[1:2, i].tolist() for i in range(n)]
                sites_lower = [self.sites_all[2:3, i].tolist() for i in range(n)]
            sites = sites_upper + sites_lower

        if self.num_channels == 3:
            sites_upper = [self.sites_all[0:1, i].tolist() for i in range(n)]
            sites_middle = [self.sites_all[1:3, i].tolist() for i in range(n)]
            sites_lower = [self.sites_all[3:4, i].tolist() for i in range(n)]
            sites = sites_upper + sites_middle + sites_lower

        if self.num_channels == 4:
            sites = [self.sites_all[0:4, i].tolist() for i in range(n)]

        return sites

    def _get_sites_id(self):
        """Generates identifiers for the monomers."""

        sites_id = [str(site[0]) for site in self.sites]
        sites_id = np.array(sites_id).reshape(self.tb_dims)
        return sites_id

    def _get_monomers(self):
        """Generates a list of monomers in shape of the TB model."""

        filepaths = [
            [os.path.join(self.directory, site + ".xyz") for site in row] for row in self.sites
        ]
        monomers = [Monomer(filepaths, self.lcao_param) for filepaths in filepaths]
        return np.array(monomers).reshape(self.tb_dims)

    # ------------------------------------------------------------------

    def calc_H(self):
        """
        Calculates the LCAO Hamiltonian matrix for the oligomer in the basis of the
        atmomic orbitals.
        """

        n = self.num_sites_per_strand
        dims = []

        # A contains the block matrices
        A = np.zeros((self.num_sites, self.num_sites), dtype=object)

        # add inter monomer blocks
        for coupling in self.couplings:
            _, monomer1_idx, monomer2_idx = coupling
            monomer1 = self.monomers[monomer1_idx]
            monomer2 = self.monomers[monomer2_idx]
            i = monomer1_idx[0] * n + monomer1_idx[1]
            j = monomer2_idx[0] * n + monomer2_idx[1]
            H_inter = calc_H_inter(self.lcao_param, monomer1, monomer2)
            A[i, j] = H_inter
            A[j, i] = H_inter.conj().T

        # add intra monomer blocks
        for energy in self.energies:
            _, monomer_idx, _ = energy
            monomer = self.monomers[monomer_idx]
            i = monomer_idx[0] * n + monomer_idx[1]
            A[i, i] = monomer.H
            dims.append(monomer.num_orbitals)

        # add zero matrices for non-coupled monomers
        for i in range(self.num_sites):
            for j in range(self.num_sites):
                if isinstance(A[i, j], int):
                    A[i, j] = np.zeros((dims[i], dims[j]))
        H = np.block(A.tolist())
        return H, dims

    def calc_tb_params2(self):
        """
        Alternative method to calculate tight-binding parameters by first calculating
        the full Hamiltonian matrix in the basis of the atomic orbitals.
        This method can be used as consistency check.
        """

        H, dims = self.calc_H()
        cum_dims = np.cumsum([0] + dims)

        n = self.num_sites_per_strand
        HOMO_dict, LUMO_dict = {}, {}

        # TB couplings
        for coupling in self.couplings:
            key, monomer1_idx, monomer2_idx = coupling
            coupling_id = (
                key + "_" + self.sites_id[monomer1_idx] + "_" + self.sites_id[monomer2_idx]
            )

            HOMO_1 = self.monomers[monomer1_idx].HOMO
            HOMO_2 = self.monomers[monomer2_idx].HOMO
            LUMO_1 = self.monomers[monomer1_idx].LUMO
            LUMO_2 = self.monomers[monomer2_idx].LUMO
            i = monomer1_idx[0] * n + monomer1_idx[1]
            j = monomer2_idx[0] * n + monomer2_idx[1]

            t_HOMO = (
                HOMO_1 @ H[cum_dims[i] : cum_dims[i + 1], cum_dims[j] : cum_dims[j + 1]] @ HOMO_2
            )
            t_LUMO = (
                LUMO_1 @ H[cum_dims[i] : cum_dims[i + 1], cum_dims[j] : cum_dims[j + 1]] @ LUMO_2
            )
            HOMO_dict[coupling_id] = t_HOMO
            LUMO_dict[coupling_id] = t_LUMO

        # TB energies
        for energy in self.energies:
            key, monomer_idx, _ = energy
            coupling_id = key + "_" + self.sites_id[monomer_idx]

            HOMO = self.monomers[monomer_idx].HOMO
            LUMO = self.monomers[monomer_idx].LUMO
            i = monomer_idx[0] * n + monomer_idx[1]

            E_HOMO = HOMO @ H[cum_dims[i] : cum_dims[i + 1], cum_dims[i] : cum_dims[i + 1]] @ HOMO
            E_LUMO = LUMO @ H[cum_dims[i] : cum_dims[i + 1], cum_dims[i] : cum_dims[i + 1]] @ LUMO
            HOMO_dict[coupling_id] = E_HOMO
            LUMO_dict[coupling_id] = E_LUMO

        tb_params = {"hole": HOMO_dict, "electron": LUMO_dict}
        return tb_params

    # ------------------------------------------------------------------

    def calc_tb_couplings(self, monomer1, monomer2):
        """Calculate TB coupling between two monomers."""

        H_inter = calc_H_inter(self.lcao_param, monomer1, monomer2)

        t_HOMO = monomer1.HOMO @ H_inter @ monomer2.HOMO
        t_LUMO = monomer1.LUMO @ H_inter @ monomer2.LUMO
        t_EXC = calc_dipolar_coupling(monomer1, monomer2)
        return round(t_HOMO, 5), round(t_LUMO, 5), round(t_EXC, 5)

    def calc_tb_energies(self, monomer):
        """Calculate TB energy of a monomer."""
        return round(monomer.E_HOMO, 3), round(monomer.E_LUMO, 3), 0

    def calc_tb_params(self):
        """Calculates the tight-binding parameters for the system."""

        if self.tb_params is not None:
            return self.tb_params

        HOMO_dict, LUMO_dict, EXC_dict = {}, {}, {}

        # TB couplings
        for coupling in self.couplings:
            key, monomer1_idx, monomer2_idx = coupling
            coupling_id = (
                key + "_" + self.sites_id[monomer1_idx] + "_" + self.sites_id[monomer2_idx]
            )

            monomer1 = self.monomers[monomer1_idx]
            monomer2 = self.monomers[monomer2_idx]
            t_HOMO, t_LUMO, t_EXC = self.calc_tb_couplings(monomer1, monomer2)
            HOMO_dict[coupling_id] = t_HOMO
            LUMO_dict[coupling_id] = t_LUMO
            EXC_dict[coupling_id] = t_EXC

        # TB energies
        for energy in self.energies:
            key, monomer_idx, _ = energy
            coupling_id = key + "_" + self.sites_id[monomer_idx]

            monomer = self.monomers[monomer_idx]
            E_HOMO, E_LUMO, E_EXC = self.calc_tb_energies(monomer)
            HOMO_dict[coupling_id] = E_HOMO
            LUMO_dict[coupling_id] = E_LUMO
            EXC_dict[coupling_id] = E_EXC

        self.tb_params = {"hole": HOMO_dict, "electron": LUMO_dict, "exciton": EXC_dict}
        return self.tb_params

    def save_tb_params(self, directory=None, filename=None):
        """Save tight-binding parameters to a specified directory or default location."""

        if filename is None:
            filename = self.filename_pdb

        tb_params = self.calc_tb_params()
        save_tb_params(
            tb_params,
            filename,
            self.tb_model_name,
            directory=directory,
            unit="eV",
        )

    # ------------------------------------------------------------------

    def plot_couplings(
        self,
        particle,
        add_colorbar=False,
        add_label=True,
        fig=None,
        ax=None,
        dpi=None,
        max_coupling=None,
    ):
        """Plots coupling interactions between sites for a given particle type."""

        lw, fs, ms = 5, 11, 22
        if fig is None:
            overhead_x = (ms - lw) / 72 - (self.num_sites_per_strand - 2) * lw / 72
            overhead_y = (ms - lw) / 72
            overhead_label = 0
            if add_label:
                overhead_label = 0.3
            pad = ms / 72
            x_size = 3.4 / 4 * (self.num_sites_per_strand - 1) + overhead_label
            y_size = (
                (x_size - overhead_x - overhead_label - 2 * pad)
                * (self.num_channels - 1)
                / (self.num_sites_per_strand - 1)
                + overhead_y
                + 2 * pad
            )
            if add_colorbar:
                x_size += 0.5

            fig, ax = plt.subplots(figsize=(x_size, y_size), dpi=dpi)

        tb_params = self.calc_tb_params()[particle]
        cmap_dict = {
            "hole": plt.cm.get_cmap("Reds"),
            "electron": plt.cm.get_cmap("Blues"),
            "exciton": plt.cm.get_cmap("Greens"),
        }
        cmap = cmap_dict[particle]
        labels = [site_id[2] for site_id in self.sites_id.flatten()]

        tb_basis = [str_to_tuple(element) for element in self.tb_basis]
        positions = [(element[1], self.num_strands - 1 - element[0]) for element in tb_basis]
        # x_coords, y_coords = zip(*positions)
        # x_margin = 0.3
        # y_margin = 0.3
        # ax.set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
        # ax.set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)

        couplings = []
        for tb_str, old_state, new_state in self.couplings:
            tb_str = f"{tb_str}_{self.sites_id[old_state]}_{self.sites_id[new_state]}"
            if tb_str[0] in ["h", "r"] and tb_str not in tb_params:
                tb_str = (
                    f"{tb_str.split('_')[0]}_{self.sites_id[new_state]}_{self.sites_id[old_state]}"
                )
            tb_val = 1e3 * abs(tb_params[tb_str])  # convert to meV

            old_idx = tuple_to_int(self.tb_dims, old_state)
            new_idx = tuple_to_int(self.tb_dims, new_state)
            couplings.append((old_idx, new_idx, tb_val))

        if max_coupling is None:
            strengths = [c[2] for c in couplings]
            max_coupling = max(strengths)
        norm = plt.Normalize(0, max_coupling)

        for i, j, strength in couplings:
            x1, y1 = positions[i]
            x2, y2 = positions[j]
            color = cmap(norm(strength))
            ax.plot([x1, x2], [y1, y2], color=color, lw=lw)

        for (x, y), label in zip(positions, labels):

            if label == "M":
                color = "#780B1C"
            else:
                color = "#4d4d4dff"
            ax.plot(x, y, "o", markersize=ms, color=color)
            ax.text(
                x,
                y,
                rf"$\mathbf{{{label}}}$",
                color="white",
                fontsize=fs,
                ha="center",
                va="center",
                weight="bold",
            )

        if add_colorbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = ax.figure.colorbar(sm, ax=ax, orientation="vertical", aspect=10)
            cbar.set_label("Interaction [meV]", fontsize=fs)
            cbar.ax.tick_params(labelsize=fs)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if add_label:
            ax.set_ylabel(particle, fontsize=fs)

        pad = ms / 72
        ax.set_xlim(-pad, self.num_sites_per_strand - 1 + pad)
        ax.set_ylim(-pad, self.num_channels - 1 + pad)
        return fig, ax

    def clean(self):
        """Clean up temporary files and directories created during processing."""

        shutil.rmtree(self.directory)


# ----------------------------------------------------------------------
