from typing import List, Optional, Type, Tuple

from utils import get_config
from .tb_basis import get_tb_basis, get_eh_basis
from .tb_config import get_tb_config
from DNA import TB_MODELS_PROPS

# Shortcuts
# tb: tight-binding
# dims: dimensions

__all__ = ['Custom_TB_Model', 'CustomTBModelType', 'TB_Model', 'TBModelType']

# ------------------------------------------------------------------------------------------------------------

class Custom_TB_Model:
    def __init__(self, tb_model_name: str, tb_dims: Tuple[int], tb_basis: List[str], tb_config: List[str]):

        assert isinstance(tb_model_name, str), "tb_model_name must be of type str"
        assert isinstance(tb_dims, tuple) and len(tb_dims)==2, "tb_dims must be of type tuple and of length two"
        assert all([isinstance(tb_dim, int) and tb_dim > 0 for tb_dim in tb_dims]), "elements of tb_dim must be of type int and > 0"
        assert isinstance(tb_basis, list), "tb_basis must be of type list"
        assert all([isinstance(tb_basis_element, str) for tb_basis_element in tb_basis]), "elements of tb_basis must be strings" 
        assert isinstance(tb_config, list), "tb_config must be of type list"
        assert all([isinstance(tb_config_element, tuple) and len(tb_config_element) == 3 for tb_config_element in tb_config]), "elements of tb_config must be of type tuple and of length three"

        self.verbose = get_config()['verbose']
        if self.verbose:
            print("Successfully checked all inputs of the Custom_TB_Model instance.")

        self.tb_model_name = tb_model_name 
        self.tb_dims = tb_dims
        self.num_strands, self.num_sites_per_strand = self.tb_dims
        self.num_sites = self.num_strands * self.num_sites_per_strand
        self.tb_config = tb_config
        self.tb_basis = tb_basis

        if self.verbose:
            print("Successfully initialized the Custom_TB_Model instance.")

    def __vars__(self) -> dict:
        return vars(self)

    def __repr__(self) -> str:
        return f"Custom_TB_Model({self.tb_model_name}, {self.tb_dims}, {self.tb_basis}, {self.tb_config})"

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

CustomTBModelType = Type[Custom_TB_Model]
        
# ---------------------------------------------------------------------------------------------------------

class TB_Model:
    def __init__(self, tb_model_name: str, tb_dims: Tuple[int]):

        assert isinstance(tb_model_name, str), "tb_model_name must be of type str"
        assert isinstance(tb_dims, tuple) and len(tb_dims)==2, "tb_dims must be of type tuple and of length two"
        assert all([isinstance(tb_dim, int) and tb_dim > 0 for tb_dim in tb_dims]), "elements of tb_dim must be of type int and > 0"
        self.verbose = get_config()['verbose']
        if self.verbose:
            print("Successfully checked all inputs of the TB_Model instance.")
            
        self.tb_model_name = tb_model_name 
        self.tb_dims = tb_dims
        self.num_strands, self.num_sites_per_strand = self.tb_dims
        self.num_sites = self.num_strands * self.num_sites_per_strand

        TB_MODELS = list(TB_MODELS_PROPS.keys())
        assert self.tb_model_name in TB_MODELS, f"tb_model_name be a predefined model {TB_MODELS}"
        self.tb_model_props = TB_MODELS_PROPS[self.tb_model_name]
        assert self.num_strands == self.tb_model_props['num_strands'], f"provided number of strands {self.num_strands} does not match the number of strands required for the predefined model {self.tb_model_props['num_strands']}"

        self.tb_config = get_tb_config(self.tb_model_name, self.tb_dims)
        self.tb_basis = get_tb_basis(self.tb_dims)
        self.eh_basis = get_eh_basis(self.tb_dims)

        if self.verbose:
            print("Successfully initialized the TB_Model instance.")

    def __vars__(self) -> dict:
        return vars(self)

    def __repr__(self) -> str:
        return f"TB_Model({self.tb_model_name}, {self.tb_dims})"

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

TBModelType = Type[TB_Model]

# ------------------------------------------------------------------------------------------------------------
