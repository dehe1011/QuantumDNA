from typing import List, Tuple
from utils import get_config
from DNA import TB_MODELS_PROPS

# Shortcuts:
# WM: wire model
# LM: ladder model
# ELM: extended ladder model
# F: fishbone
# FC: fully connected 

# __all__ = ['get_tb_config', 'TB_MODELS_PROPS']

# ----------------------------------------------------------------------------------------------------------------------

def get_WM_config(num_sites_per_strand: int, strand: int = 0, reversed_direction: bool = False) -> List[Tuple[str]]:
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


def get_LM_config(num_sites_per_strand: int, strand1: int = 0, strand2: int = 1) -> List[Tuple[str]]:
    WM_config_strand1 = get_WM_config(num_sites_per_strand, strand=strand1)
    WM_config_strand2 = get_WM_config(
        num_sites_per_strand, strand=strand2, reversed_direction=True
    )
    h_list = [
        ("h", str((strand1, site)), str((strand2, site)))
        for site in range(num_sites_per_strand)
    ]
    return WM_config_strand1 + WM_config_strand2 + h_list


def get_ELM_config(num_sites_per_strand: int, strand1: int = 0, strand2: int = 1) -> List[Tuple[str]]:

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


def get_F_config(num_sites_per_strand: int, strand1: int = 0, strand2: int = 2) -> List[Tuple[str]]:

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


def get_FWM_config(num_sites_per_strand: int) -> List[Tuple[str]]:

    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=2)
    WM_config = get_WM_config(num_sites_per_strand, strand=1)
    return F_config + WM_config


def get_FLM_config(num_sites_per_strand: int) -> List[Tuple[str]]:
    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=3)
    LM_config = get_LM_config(num_sites_per_strand, strand1=1, strand2=2)
    return F_config + LM_config


def get_FELM_config(num_sites_per_strand: int) -> List[Tuple[str]]:

    F_config = get_F_config(num_sites_per_strand, strand1=0, strand2=3)
    ELM_config = get_ELM_config(num_sites_per_strand, strand1=1, strand2=2)
    return F_config + ELM_config


def get_FC_config(num_sites_per_strand: int) -> List[Tuple[str]]:

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

# --------------------------------------------------------------

TB_MODELS = get_config()['TB_MODELS']
TB_CONFIGS = dict(zip(TB_MODELS,[get_WM_config, get_LM_config, get_ELM_config, get_FWM_config, get_FLM_config, get_FELM_config, get_FC_config]))

def get_tb_config(tb_model_name: str, tb_dims: Tuple[int]) -> List[Tuple[str]]:
    tb_model_name = tb_model_name.upper()
    _, num_sites_per_strand = tb_dims
    return TB_CONFIGS[tb_model_name](num_sites_per_strand)
    