import os
import shutil
import copy

import numpy as np
import scipy.constants as c
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
from .monomer import Monomer

# ----------------------------------------------------------------------


def calc_dipolar_coupling(monomer1, monomer2):
    R1, R2 = monomer1.center_of_mass, monomer2.center_of_mass
    if np.allclose(R1, R2):
        return 0
    mu1, mu2 = monomer1.dipole_moment, monomer2.dipole_moment
    prefactor = 1 / (4 * np.pi * c.epsilon_0) * 1 / np.linalg.norm(R1 - R2) ** 3
    return (
        prefactor
        * (
            mu1 @ mu2
            - 3 / np.linalg.norm(R1 - R2) ** 2 * (mu1 @ (R2 - R1)) * (mu2 @ (R2 - R1))
        )
        * 1e10
        / c.e
    )


class Oligomer(TBModel):
    def __init__(self, filepath_pdb, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(DEFAULTS["lcao_kwargs_default"])
        self.kwargs.update(kwargs)
        check_lcao_kwargs(**self.kwargs)

        self.filepath_pdb = filepath_pdb
        self.directory = self.filepath_pdb.split(".")[0]

        # if self.directory exists delete it
        if os.path.exists(self.directory):
            shutil.rmtree(self.directory)

        pdb_to_xyz(self.filepath_pdb, **self.kwargs)
        self.xyz_filenames = find_xyz(self.directory)

        num_sites_per_strand = len(self.xyz_filenames) // 4
        super().__init__(num_sites_per_strand, **self.kwargs)

        # arguments
        self.backbone = self.tb_model_props["backbone"]
        self.auto_clean = self.kwargs.get("auto_clean")
        self.param_id = self.kwargs.get("param_id")
        self.lcao_param = load_lcao_param(self.param_id)

        self.filename_pdb = os.path.splitext(os.path.basename(self.filepath_pdb))[0]
        self.sites_bases, self.sites_backbone = self._get_sites_all()
        self.sites, self.sites_id = self._get_sites()
        self.sites_id = np.array(self.sites_id).reshape(self.tb_dims)
        self.filepaths = self._get_filepaths()

        self.num_sites = self.num_channels * self.num_sites_per_strand
        self.monomers = [
            Monomer(filepaths, **self.kwargs) for filepaths in self.filepaths
        ]
        self.monomers = np.array(self.monomers).reshape(self.tb_dims)

        self.couplings = get_tb_couplings(self.tb_model_name, self.num_sites_per_strand)
        self.energies = get_tb_energies(self.num_channels, self.num_sites_per_strand)

        self.tb_params = None

        if self.auto_clean:
            self.clean()

    def __repr__(self):
        return f"Oligomer({self.filename_pdb})"

    def _get_filepaths(self):
        return [
            [os.path.join(self.directory, site + ".xyz") for site in row]
            for row in self.sites
        ]

    def _get_sites_all(self):
        sites_bases, sites_backbone = [], []
        for filename in self.xyz_filenames:
            if "B" in filename:
                sites_backbone.append(filename)
            else:
                sites_bases.append(filename)

        return sites_bases, sites_backbone

    def _get_sites(self):
        n = len(self.sites_bases) // 2

        sites_all = np.array(
            [
                self.sites_backbone[:n],
                self.sites_bases[:n],
                self.sites_bases[n:][::-1],
                self.sites_backbone[n:][::-1],
            ],
            dtype=object,
        )
        sites, sites_id = [], []
        if self.num_channels == 1:
            if self.backbone:
                sites = [sites_all[:, i] for i in range(n)]
            else:
                sites = [sites_all[1:3, i] for i in range(n)]
            sites_id = [sites_all[k, i] for k in [1] for i in range(n)]

        if self.num_channels == 2:
            if self.backbone:
                sites_upper = [sites_all[:2, i] for i in range(n)]
                sites_lower = [sites_all[2:, i] for i in range(n)]
            else:
                sites_upper = [[sites_all[1, i]] for i in range(n)]
                sites_lower = [[sites_all[2, i]] for i in range(n)]
            sites = sites_upper + sites_lower
            sites_id = [sites_all[k, i] for k in [1, 2] for i in range(n)]

        if self.num_channels == 3:
            sites_upper = [sites_all[:1, i] for i in range(n)]
            sites_middle = [sites_all[1:3, i] for i in range(n)]
            sites_lower = [sites_all[3:, i] for i in range(n)]
            sites = sites_upper + sites_middle + sites_lower
            sites_id = [sites_all[k, i] for k in [0, 1, 3] for i in range(n)]

        if self.num_channels == 4:
            sites = [[site] for site in sites_all.flatten()]
            sites_id = [sites_all[k, i] for k in [0, 1, 2, 3] for i in range(n)]

        return sites, sites_id

    def calc_tb_couplings(self, monomer1, monomer2):
        H_inter = calc_H_inter(self.lcao_param, monomer1, monomer2)

        t_HOMO = monomer1.HOMO @ H_inter @ monomer2.HOMO
        t_LUMO = monomer1.LUMO @ H_inter @ monomer2.LUMO
        t_EXC = calc_dipolar_coupling(monomer1, monomer2)
        return round(t_HOMO, 3), round(t_LUMO, 3), round(t_EXC, 3)

    def calc_tb_energies(self, monomer):
        return round(monomer.E_HOMO, 3), round(monomer.E_LUMO, 3), 0

    def calc_tb_params(self):

        if self.tb_params is not None:
            return self.tb_params

        HOMO_dict, LUMO_dict, EXC_dict = {}, {}, {}
        for coupling in self.couplings:
            key, monomer1_idx, monomer2_idx = coupling
            coupling_id = (
                key
                + "_"
                + self.sites_id[monomer1_idx]
                + "_"
                + self.sites_id[monomer2_idx]
            )

            monomer1 = self.monomers[monomer1_idx]
            monomer2 = self.monomers[monomer2_idx]
            t_HOMO, t_LUMO, t_EXC = self.calc_tb_couplings(monomer1, monomer2)
            HOMO_dict[coupling_id] = t_HOMO
            LUMO_dict[coupling_id] = t_LUMO
            EXC_dict[coupling_id] = t_EXC

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

    def save_tb_params(self, directory=None):
        tb_params = self.calc_tb_params()
        save_tb_params(
            tb_params,
            self.filename_pdb,
            self.tb_model_name,
            directory=directory,
            unit="eV",
        )

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

        lw, fs, fs2, ms = 5, 11, 22, 22
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
        positions = [
            (element[1], self.num_strands - 1 - element[0]) for element in tb_basis
        ]
        # x_coords, y_coords = zip(*positions)
        # x_margin = 0.3
        # y_margin = 0.3
        # ax.set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
        # ax.set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)

        couplings = []
        for tb_str, old_state, new_state in self.couplings:
            tb_str = f"{tb_str}_{self.sites_id[old_state]}_{self.sites_id[new_state]}"
            if tb_str[0] in ["h", "r"] and tb_str not in tb_params:
                tb_str = f"{tb_str.split('_')[0]}_{self.sites_id[new_state]}_{self.sites_id[old_state]}"
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
        shutil.rmtree(self.directory)


# ----------------------------------------------------------------------
