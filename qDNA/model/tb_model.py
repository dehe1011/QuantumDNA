import copy

from ..utils import check_tb_model_kwargs
from ..io import DEFAULTS, OPTIONS, load_tb_model_props
from .tb_basis import get_tb_basis, get_eh_basis
from .tb_config import get_tb_config

# ----------------------------------------------------------------------


class TBModel:
    def __init__(self, num_sites_per_strand, **kwargs):

        # check kwargs
        self.kwargs = copy.copy(kwargs)
        self.kwargs.update(DEFAULTS["tb_model_kwargs_default"])
        self.kwargs.update(kwargs)
        check_tb_model_kwargs(**self.kwargs)

        self.num_sites_per_strand = num_sites_per_strand
        self.tb_model_name = self.kwargs.get("tb_model_name")

        if self.tb_model_name in OPTIONS["tb_models"]:
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
