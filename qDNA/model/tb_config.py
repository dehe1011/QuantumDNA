"""
This module provides functions to generate configurations for different tight-binding models.
"""

from qDNA.tools import get_config
from qDNA import TB_MODELS_PROPS

__all__ = ["get_tb_config", "TB_MODELS_PROPS"]

# Shortcuts:
# WM: wire model
# LM: ladder model
# ELM: extended ladder model
# F: fishbone
# FC: fully connected

# --------------------------------------------------------------------------


def get_wm_config(num_sites_per_strand, strand=0, reversed_direction=False):
    """
    Generate the configuration for the Wire Model (WM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.
    strand : int, optional
        The strand index, by default 0.
    reversed_direction : bool, optional
        If True, reverses the direction of the connections, by default False.

    Returns
    -------
    list of tuple
        The configuration of the Wire Model.
    """
    e_list = [
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
    return e_list + t_list


def get_lm_config(num_sites_per_strand, strand1=0, strand2=1):
    """
    Generate the configuration for the Ladder Model (LM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.
    strand1 : int, optional
        The first strand index, by default 0.
    strand2 : int, optional
        The second strand index, by default 1.

    Returns
    -------
    list of tuple
        The configuration of the Ladder Model.
    """
    wm_config_strand1 = get_wm_config(num_sites_per_strand, strand=strand1)
    wm_config_strand2 = get_wm_config(
        num_sites_per_strand, strand=strand2, reversed_direction=True
    )
    h_list = [
        ("h", str((strand1, site)), str((strand2, site)))
        for site in range(num_sites_per_strand)
    ]
    return wm_config_strand1 + wm_config_strand2 + h_list


def get_elm_config(num_sites_per_strand, strand1=0, strand2=1):
    """
    Generate the configuration for the Extended Ladder Model (ELM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.
    strand1 : int, optional
        The first strand index, by default 0.
    strand2 : int, optional
        The second strand index, by default 1.

    Returns
    -------
    list of tuple
        The configuration of the Extended Ladder Model.
    """
    lm_config = get_lm_config(num_sites_per_strand, strand1=strand1, strand2=strand2)
    r_plus_list = [
        ("r+", str((strand1, site)), str((strand2, site + 1)))
        for site in range(num_sites_per_strand - 1)
    ]
    r_minus_list = [
        ("r-", str((strand1, site + 1)), str((strand2, site)))
        for site in range(num_sites_per_strand - 1)
    ]
    return lm_config + r_plus_list + r_minus_list


def get_f_config(num_sites_per_strand, strand1=0, strand2=2):
    """
    Generate the configuration for the Fishbone Model (F).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.
    strand1 : int, optional
        The first strand index, by default 0.
    strand2 : int, optional
        The second strand index, by default 2.

    Returns
    -------
    list of tuple
        The configuration of the Fishbone Model.
    """
    e_list = [
        ("E", str((strand1, tb_site)), str((strand1, tb_site)))
        for tb_site in range(num_sites_per_strand)
    ]
    e_list += [
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


def get_fwm_config(num_sites_per_strand):
    """
    Generate the configuration for the Fishbone Wire Model (FWM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.

    Returns
    -------
    list of tuple
        The configuration of the Fishbone Wire Model.
    """
    f_config = get_f_config(num_sites_per_strand, strand1=0, strand2=2)
    wm_config = get_wm_config(num_sites_per_strand, strand=1)
    return f_config + wm_config


def get_flm_config(num_sites_per_strand):
    """
    Generate the configuration for the Fishbone Ladder Model (FLM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.

    Returns
    -------
    list of tuple
        The configuration of the Fishbone Ladder Model.
    """
    f_config = get_f_config(num_sites_per_strand, strand1=0, strand2=3)
    lm_config = get_lm_config(num_sites_per_strand, strand1=1, strand2=2)
    return f_config + lm_config


def get_felm_config(num_sites_per_strand):
    """
    Generate the configuration for the Fishbone Extended Ladder Model (FELM).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.

    Returns
    -------
    list of tuple
        The configuration of the Fishbone Extended Ladder Model.
    """
    f_config = get_f_config(num_sites_per_strand, strand1=0, strand2=3)
    elm_config = get_elm_config(num_sites_per_strand, strand1=1, strand2=2)
    return f_config + elm_config


def get_fc_config(num_sites_per_strand):
    """
    Generate the configuration for the Fully Connected Model (FC).

    Parameters
    ----------
    num_sites_per_strand : int
        The number of sites per strand.

    Returns
    -------
    list of tuple
        The configuration of the Fully Connected Model.
    """
    felm_config = get_felm_config(num_sites_per_strand)
    t_list1 = [
        ("t", str((0, tb_site)), str((0, tb_site + 1)))
        for tb_site in range(num_sites_per_strand - 1)
    ]
    t_list2 = [
        ("t", str((3, tb_site)), str((3, tb_site + 1)))
        for tb_site in range(num_sites_per_strand - 1)
    ]
    return felm_config + t_list1 + t_list2


# --------------------------------------------------------------

TB_MODELS = get_config()["TB_MODELS"]
TB_CONFIGS = dict(
    zip(
        TB_MODELS,
        [
            get_wm_config,
            get_lm_config,
            get_elm_config,
            get_fwm_config,
            get_flm_config,
            get_felm_config,
            get_fc_config,
        ],
    )
)


def get_tb_config(tb_model_name, tb_dims):
    """
    Get the tight-binding configuration for a specified model.

    Parameters
    ----------
    tb_model_name : str
        The name of the tight-binding model.
    tb_dims : tuple
        The dimensions of the model (number of strands, number of sites per strand).

    Returns
    -------
    list of tuple
        The configuration of the specified tight-binding model.
    """
    tb_model_name = tb_model_name.upper()
    _, num_sites_per_strand = tb_dims
    return TB_CONFIGS[tb_model_name](num_sites_per_strand)
