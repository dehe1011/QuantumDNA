import os

from .io import pdb_to_xyz, save_tb_params, DATA_DIR
from .lcao import Oligomer
from .io import load_tb_params
from .hamiltonian import get_tb_sites, TBHam
from .environment import LindDiss
from .dynamics import MeSolver
from .evaluation import Evaluation
from .visualization import Visualization

# ----------------------------------------------------------------------


def convert_pdb_to_xyz(filepath_pdb):
    pdb_to_xyz(filepath_pdb)


def calc_tb_params(directories, tb_model, double_stranded=True):
    HOMO_dict, LUMO_dict = {}, {}
    for directory in directories:
        oligomer = Oligomer(
            directory + ".pdb", tb_model_name=tb_model, auto_clean=False
        )
        tb_params = oligomer.calc_tb_params()
        HOMO_dict_new, LUMO_dict_new = tb_params["hole"], tb_params["electron"]
        HOMO_dict.update(HOMO_dict_new)
        LUMO_dict.update(LUMO_dict_new)
    return HOMO_dict, LUMO_dict


def wrap_save_tb_params(
    tb_params, source, particle, tb_model_name, unit=None, notes=None
):
    unit = "eV"
    tb_params = {particle: tb_params}
    directory = os.path.join(DATA_DIR, "tb_params")
    filename = f"{source}_{tb_model_name}"
    if os.path.exists(os.path.join(directory, filename + ".json")):
        tb_params_old = load_tb_params(source, tb_model_name)
        tb_params_new = {**tb_params_old, **tb_params}
    else:
        tb_params_new = tb_params
    save_tb_params(
        tb_params_new,
        source,
        tb_model_name,
        directory=None,
        unit=unit,
        notes=notes,
        override=True,
    )


class DNA_Seq:
    def __init__(
        self, upper_strand, tb_model_name, methylated=True, lower_strand="auto_complete"
    ):
        if lower_strand == "auto_complete":
            lower_strand = "auto complete"
        self.tb_sites = get_tb_sites(
            upper_strand, tb_model_name=tb_model_name, lower_strand=lower_strand
        )


class TB_Ham(TBHam):
    def __init__(self, dna_seq, **kwargs):
        super().__init__(dna_seq.tb_sites, **kwargs)


class Lindblad_Diss(LindDiss):
    def __init__(self, tb_ham, **kwargs):
        super().__init__(tb_ham.tb_sites, **{**tb_ham.kwargs, **kwargs})


class ME_Solver(MeSolver):
    def __init__(self, tb_ham, lindblad_diss, **kwargs):
        super().__init__(lindblad_diss.tb_sites, **{**lindblad_diss.kwargs, **kwargs})


def plot_pops_heatmap(me_solver):
    vis = Visualization(me_solver.tb_sites, **me_solver.kwargs)
    return vis.plot_heatmap()


def calc_lifetime(upper_strand, tb_model, **kwargs):
    if "lower_strand" not in kwargs:
        kwargs["lower_strand"] = "auto complete"

    tb_sites = get_tb_sites(
        upper_strand, tb_model_name=tb_model, lower_strand=kwargs["lower_strand"]
    )
    eva = Evaluation(tb_sites, **kwargs)
    return eva.calc_lifetime()


def calc_dipole(upper_strand, tb_model, **kwargs):
    if "lower_strand" not in kwargs:
        kwargs["lower_strand"] = "auto complete"

    tb_sites = get_tb_sites(
        upper_strand, tb_model_name=tb_model, lower_strand=kwargs["lower_strand"]
    )
    eva = Evaluation(tb_sites, **kwargs)
    return eva.calc_charge_separation()


# ----------------------------------------------------------------------
