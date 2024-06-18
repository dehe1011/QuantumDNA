from .tb_basis import get_tb_basis
from .tb_config import TB_MODELS, get_tb_config, get_num_strands
from typing import List, Optional, Type, Tuple

# ------------------------------------------------------------------------------------------------------------

class TB_Model:
    """
    Represents a tight-binding model.

    Attributes:
        tb_model_name (str): The name of the tight-binding model.
            Example: 'WM', 'LM', 'ELM'
        tb_sites (list): A list containing the string representation of the tight-binding sites.
            Example: ['A', 'G', 'C'] or 'AGC'
        tb_site_basis (list, optional): A list containing the string representation of each basis element.
            Example: ['(0, 0)', '(0, 1)', ...]
        tb_config (list, optional): A list containing the string representation of the tight-binding parameter and the basis elements it connects.
            Example: [('t_GC', ('0', '2')), ...]
        num_strands (int, optional): Number of strands in the tight-binding model.

    Notes:
        Usually the DNA string is used for tb_sites.
    """

    def __init__(
        self,
        tb_model_name: str,
        tb_dims: Tuple[int, int],
        tb_basis: Optional[List[str]] = None,
        tb_config: Optional[List[str]] = None,
    ):
        """
        Initialize a TB_Model instance.

        Args:
            num_tb_sites (int): Number of tight-binding sites.
            tb_model_name (str): Name of the tight-binding model.
            tb_site_basis (list, optional): String representation of each basis element. Defaults to None.
            tb_config (list, optional): String representation of the tight-binding parameters and the basis elements they connect. Defaults to None.
            num_strands (int, optional): Number of strands in the tight-binding model. Required if the model is not predefined. Defaults to None.
        """
        self.tb_model_name = tb_model_name
        self.tb_dims = tb_dims
        self.num_strands, self.num_sites_per_strand = self.tb_dims
        self.num_sites = self.num_strands * self.num_sites_per_strand
        if self.num_strands < 1:
            raise ValueError("Number of tight-binding strands must be >= 1")
        if self.num_sites_per_strand < 1:
            raise ValueError("Number of tight-binding sites per strand must be >= 1")

        if self.tb_model_name in TB_MODELS:
            if self.num_strands != get_num_strands(self.tb_model_name):
                print(f"Note: Changed number of strands from input {self.num_strands} to {get_num_strands(self.tb_model_name)} as required for the selected predefined model.")
                self.num_strands = get_num_strands(self.tb_model_name)

            self.tb_config = get_tb_config(self.num_sites, self.num_strands, self.tb_model_name)
            self.tb_basis = get_tb_basis(self.tb_dims)
        else:
            if not tb_config:
                raise ValueError("Please provide a tight-binding configuration.")
            self.tb_config = tb_config
            self.tb_basis = ( tb_basis if tb_basis else get_tb_basis(self.tb_dims) )

    def __str__(self) -> str:
        return vars(self)

TBModelType = Type[TB_Model]

# ------------------------------------------------------------------------------------------------------------
