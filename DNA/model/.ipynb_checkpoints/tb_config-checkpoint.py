from typing import List, Tuple

__all__ = ["TB_MODELS", "get_tb_config", "get_num_strands"]

""" 
Predefined tight-binding models
If you want to add a tight-binding model you have to provide a tb_model_name, tb_config and num_strands
"""
TB_MODELS = ["WM", "LM", "ELM", "FWM", "FLM", "FELM", "FC"]

# ----------------------------------------------------------------------------------------------------------------------


def get_WM_config(
    num_sites_per_strand: int, strand: int = 0, reversed_direction: bool = False
) -> List[Tuple[str, str, str]]:
    """
    Generate a WM (wire model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.
        strand (int, optional): Strand index. Defaults to 0.
        reversed_direction (bool, optional): Whether the direction is reversed. Defaults to False.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    E_list = [
        ("E", str((strand, tb_site)), str((strand, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    t_list = [
        (
            "t",
            str((strand, tb_site + (1 if not reversed_direction else 0))),
            str((strand, tb_site + (0 if not reversed_direction else 1))),
        )
        for tb_site in range(num_sites_per_strand - 1)
    ]
    return E_list + t_list


def get_LM_config(
    num_sites_per_strand: int, strand1: int = 0, strand2: int = 1
) -> List[Tuple[str, str, str]]:
    """
    Generate an LM (ladder model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.
        strand1 (int, optional): First strand index. Defaults to 0.
        strand2 (int, optional): Second strand index. Defaults to 1.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    WM_config_strand1 = get_WM_config(num_sites_per_strand, strand=strand1)
    WM_config_strand2 = get_WM_config(
        num_sites_per_strand, strand=strand2, reversed_direction=True
    )
    h_list = [
        ("h", str((strand1, site)), str((strand2, site)))
        for site in range(num_sites_per_strand)
    ]
    return WM_config_strand1 + WM_config_strand2 + h_list


def get_ELM_config(
    num_sites_per_strand: int, strand1: int = 0, strand2: int = 1
) -> List[Tuple[str, str, str]]:
    """
    Generate an ELM (extended ladder model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.
        strand1 (int, optional): First strand index. Defaults to 0.
        strand2 (int, optional): Second strand index. Defaults to 1.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    LM_config = get_LM_config(num_sites_per_strand, strand1=strand1, strand2=strand2)
    r_plus_list = [
        ("r+", str((strand1, site)), str((strand2, site + 1)))
        for site in range(num_sites_per_strand - 1)
    ]
    r_minus_list = [
        ("r-", str((strand1, site + 1)), str((strand2, site)))
        for site in range(num_sites_per_strand - 1)
    ]
    return LM_config + r_plus_list + r_minus_list


def get_F_config(
    num_sites_per_strand: int, strand1: int = 0, strand2: int = 2
) -> List[Tuple[str, str, str]]:
    """
    Generate an F (fishbone) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.
        strand1 (int, optional): First strand index. Defaults to 0.
        strand2 (int, optional): Second strand index. Defaults to 2.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    E_list = [
        ("E", str((strand1, tb_site)), str((strand1, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    E_list += [
        ("E", str((strand2, tb_site)), str((strand2, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    h_list1 = [
        ("h", str((strand1, tb_site)), str((strand1 + 1, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    h_list2 = [
        ("h", str((strand2 - 1, tb_site)), str((strand2, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    return h_list1 + h_list2


def get_FWM_config(num_sites_per_strand: int) -> List[Tuple[str, str, str]]:
    """
    Generate an FWM (fishbone wire model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=2)
    WM_config = get_WM_config(num_sites_per_strand, strand=1)
    return F_config + WM_config


def get_FLM_config(num_sites_per_strand: int) -> List[Tuple[str, str, str]]:
    """
    Generate an FLM (fishbone ladder model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=3)
    LM_config = get_LM_config(num_sites_per_strand, strand1=1, strand2=2)
    return F_config + LM_config


def get_FELM_config(num_sites_per_strand: int) -> List[Tuple[str, str, str]]:
    """
    Generate an FELM (fishbone extended ladder model) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=3)
    ELM_config = get_ELM_config(num_sites_per_strand, strand1=1, strand2=2)
    return F_config + ELM_config


def get_FC_config(num_sites_per_strand: int) -> List[Tuple[str, str, str]]:
    """
    Generate an FC (combined model including backbone charge transfer) configuration for tight-binding.

    Args:
        num_sites_per_strand (int): Number of sites per strand.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.
    """
    FELM_config = get_FELM_config(num_sites_per_strand)
    t_list1 = [
        ("t", str((0, tb_site)), str((0, tb_site + 1)))
        for tb_site in range(num_sites_per_strand - 1)
    ]
    t_list2 = [
        ("t", str((3, tb_site)), str((3, tb_site + 1)))
        for tb_site in range(num_sites_per_strand - 1)
    ]
    return FELM_config + t_list1 + t_list2


# ----------------------------------------------------------------------------------------------------------------------

TB_CONFIGS = dict(
    zip(
        TB_MODELS,
        [
            get_WM_config,
            get_LM_config,
            get_ELM_config,
            get_FWM_config,
            get_FLM_config,
            get_FELM_config,
            get_FC_config,
        ],
    )
)


def get_tb_config(
    num_tb_sites: int, num_strands: int, tb_model_name: str
) -> List[Tuple[str, str, str]]:
    """
    Get the configuration for a specified tight-binding model.

    Args:
        num_tb_sites (int): Number of tight-binding sites.
        num_strands (int): Number of strands.
        tb_model_name (str): Name of the tight-binding model.

    Returns:
        List[Tuple[str, str, str]]: List of configurations.

    Raises:
        ValueError: If an unknown tight-binding model name is provided.
    """
    if num_tb_sites < 1:
        raise ValueError("Number of tight-binding sites must be >= 1")
    if num_strands < 1:
        raise ValueError("Number of strands must be >= 1")

    tb_model_name = tb_model_name.upper()
    num_sites_per_strand = num_tb_sites // num_strands

    if tb_model_name in TB_MODELS:
        return TB_CONFIGS[tb_model_name](num_sites_per_strand)
    else:
        raise ValueError(
            f"Unknown tight-binding model name: {tb_model_name}. Predefined models: {TB_MODELS}"
        )


# ----------------------------------------------------------------------------------------------------------------------

NUM_STRANDS = dict(zip(TB_MODELS, [1, 2, 2, 3, 4, 4, 4]))


def get_num_strands(tb_model_name: str) -> int:
    """
    Get the number of strands for a specified tight-binding model.

    Args:
        tb_model_name (str): Name of the tight-binding model.

    Returns:
        int: Number of strands.

    Raises:
        ValueError: If an unknown tight-binding model name is provided.
    """

    if tb_model_name in TB_MODELS:
        return NUM_STRANDS[tb_model_name]
    else:
        raise ValueError(
            f"Unknown tight-binding model name: {tb_model_name}. Predefined models: {TB_MODELS}"
        )
