from utils import my_save, my_load
from typing import Dict, Optional
from copy import deepcopy

__all__ = [
    "save_tb_params",
    "load_tb_params",
    "wrap_save_tb_params",
    "wrap_load_tb_params",
]


def save_tb_params(
    tb_param_dict: Dict,
    info_dict: Dict[str, str],
    directory: str = "stored_data/tb_params",
    notes: Optional[str] = None,
) -> None:
    """
    Save tight-binding parameters to a file.

    Args:
        tb_param_dict (dict): Dictionary containing the tight-binding parameters.
        info_dict (dict): Dictionary with metadaata, e.g., 'source', 'particle', and 'tb_model_name'.
        directory (str, optional): Directory to save the file. Defaults to 'stored_data/tb_params'.
        notes (str, optional): Additional notes to include in the info_dict. Defaults to None.

    Example:
        >>> save_tb_params(
                {'couple': 5, 'bff': 3, 'enemies': -2, 'very_cool': 3.5},
                {'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'WM'},
                directory='stored_data/my_params',
                notes='here you can add some notes.'
            )
    """
    filename = "_".join(info_dict.values())
    info_dict_notes = deepcopy(info_dict)
    if notes:
        info_dict_notes["notes"] = notes
    my_save(tb_param_dict, info_dict_notes, filename, directory=directory, save_excel=False)


def load_tb_params(
    info_dict: Dict[str, str],
    directory: str = "stored_data/tb_params",
    load_metadata: bool = False,
) -> Dict:
    """
    Load tight-binding parameters from a file.

    Args:
        info_dict (dict): Dictionary with metadaata, e.g., 'source', 'particle', and 'tb_model_name'.
        directory (str, optional): Directory to load the file from. Defaults to 'stored_data/tb_params'.
        load_metadata (bool, optional): Whether to load metadata along with the parameters. Defaults to False.

    Returns:
        dict: Loaded tight-binding parameters.

    Example:
        >>> load_tb_params(
                {'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'WM'},
                directory='stored_data/my_params',
                load_metadata=True
            )
    """
    filename = "_".join(info_dict.values())
    return my_load(filename, load_metadata=load_metadata, directory=directory)


def wrap_save_tb_params(
    tb_param_dict: Dict,
    source: str,
    particle: str,
    tb_model_name: str,
    directory: str = "stored_data/tb_params",
    notes: Optional[str] = None,
) -> None:
    """
    Wrapper function for save_tb_params().

    Args:
        tb_param_dict (dict): Dictionary containing the tight-binding parameters.
        source (str): Source of the parameters, e.g., 'Hawke2010'.
        particle (str): Type of particle, e.g., 'electron', 'hole', 'exciton'.
        tb_model_name (str): Name of the tight-binding model.
        directory (str, optional): Directory to save the file. Defaults to 'stored_data/tb_params'.
        notes (str, optional): Additional notes to include in the info_dict. Defaults to None.
    """
    info_dict = {"source": source, "particle": particle, "tb_model_name": tb_model_name}
    save_tb_params(tb_param_dict, info_dict, directory=directory, notes=notes)


def wrap_load_tb_params(
    source: str,
    particle: str,
    tb_model_name: str,
    directory: str = "stored_data\\tb_params",
    load_metadata: bool = False,
) -> Dict:
    """
    Wrapper function for load_tb_params().

    Args:
        source (str): Source of the parameters, e.g., 'Hawke2010'.
        particle (str): Type of particle, e.g., 'electron', 'hole', 'exciton'.
        tb_model_name (str): Name of the tight-binding model.
        directory (str, optional): Directory to load the file from. Defaults to 'stored_data/tb_params'.
        load_metadata (bool, optional): Whether to load metadata along with the parameters. Defaults to False.

    Returns:
        dict: Loaded tight-binding parameters.
    """
    info_dict = {"source": source, "particle": particle, "tb_model_name": tb_model_name}
    return load_tb_params(info_dict, directory=directory, load_metadata=load_metadata)
