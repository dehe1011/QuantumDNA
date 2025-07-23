import copy

from ..utils import check_tb_model_kwargs
from ..io import DEFAULTS, OPTIONS, load_tb_model_props
from .tb_basis import get_tb_basis, get_eh_basis
from .tb_config import get_tb_config

# ----------------------------------------------------------------------


class TBModel:
    """
    Tight-Binding Model class for representing DNA-like structures.

    Parameters
    ----------
    num_sites_per_strand : int
        Number of sites per strand in the model.

    Attributes
    ----------
    num_sites_per_strand : int
        Number of sites per strand in the model.
    tb_model_name : str
        Name of the tight-binding model.
    tb_model_props : dict
        Properties of the selected tight-binding model.
    backbone : str
        Backbone type of the model.
    double_stranded : bool
        Indicates if the model is double-stranded.
    num_channels : int
        Number of channels (strands) in the model.
    num_strands : int
        Number of strands in the model (same as `num_channels`).
    tb_dims : tuple
        Dimensions of the tight-binding model (channels, sites per strand).
    tb_config : dict
        Configuration of the tight-binding model.
    num_sites : int
        Total number of sites in the model.
    tb_basis : ndarray
        Basis for the tight-binding model.
    eh_basis : ndarray
        Electron-hole basis for the model.

    Notes
    -----
    .. note::

        - If a custom tight-binding model is used, `kwargs` must include `num_channels` and `tb_config`.

    """

    def __init__(self, num_sites_per_strand, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(kwargs)
        self.kwargs.update(DEFAULTS["tb_model_kwargs_default"])
        self.kwargs.update(kwargs)

        self.num_sites_per_strand = num_sites_per_strand
        self.tb_model_name = self.kwargs.get("tb_model_name")

        if self.tb_model_name in OPTIONS["tb_models"]:
            check_tb_model_kwargs(**self.kwargs)
            self.tb_model_props = load_tb_model_props()[self.tb_model_name]
            self.backbone = self.tb_model_props["backbone"]
            self.double_stranded = self.tb_model_props["double_stranded"]
            self.num_channels = self.tb_model_props["num_strands"]
            self.num_strands = self.num_channels
            self.tb_dims = (self.num_channels, self.num_sites_per_strand)
            self.tb_config = get_tb_config(self.tb_model_name, self.tb_dims)

        else:
            print(f"Using custom TB model: {self.tb_model_name}")
            self.num_channels = self.kwargs["num_channels"]
            self.tb_dims = (self.num_channels, self.num_sites_per_strand)
            self.tb_config = self.kwargs["tb_config"]

        self.num_sites = self.num_channels * self.num_sites_per_strand
        self.tb_basis = get_tb_basis(self.tb_dims)
        self.eh_basis = get_eh_basis(self.tb_dims)


# ----------------------------------------------------------------------
