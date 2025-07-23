import os

from .. import DATA_DIR
from .io_json import load_json, save_json, modify_json

# -----------------------------------------------------------------------


def save_tb_params(
    tb_params,
    source,
    tb_model_name,
    directory=None,
    unit="meV",
    notes=None,
    override=False,
):
    """
    Save tight-binding parameters to a JSON file and update metadata.

    Parameters
    ----------
    tb_params : dict
        Dictionary containing tight-binding parameters.
    source : str
        Identifier for the source of the parameters.
    tb_model_name : str
        Name of the tight-binding model.
    directory : str, optional
        Directory to save the JSON file. Defaults to a subdirectory in `DATA_DIR`.
    unit : str, optional
        Unit of the parameters. Defaults to "meV".
    notes : str, optional
        Additional notes about the parameters. Defaults to "No notes provided."
    override : bool, optional
        Whether to override the file if it already exists. Defaults to False.

    Returns
    -------
    None

    """

    if directory is None:
        directory = os.path.join(DATA_DIR, "tb_params")
    if notes is None:
        notes = "No notes provided."
    metadata = {
        "source": source,
        "tb_model_name": tb_model_name,
        "unit": unit,
        "notes": notes,
    }
    data = {"data": tb_params, "metadata": metadata}

    filename = "_".join([source, tb_model_name])
    filepath = os.path.join(directory, f"{filename}.json")
    save_json(data, filepath, override=override)
    modify_json(os.path.join(DATA_DIR, "options.json"), "sources", metadata["source"])


def load_tb_params(
    source,
    tb_model_name,
    directory=None,
    load_metadata=False,
):
    """
    Load tight-binding parameters from a JSON file.

    Parameters
    ----------
    source : str
        The source identifier for the parameters.
    tb_model_name : str
        The name of the tight-binding model.
    directory : str, optional
        The directory containing the parameter files. Defaults to a subdirectory
        within `DATA_DIR` named "tb_params".
    load_metadata : bool, optional
        If True, returns both data and metadata. Defaults to False.

    Returns
    -------
    dict or tuple
        If `load_metadata` is False, returns the parameter data as a dictionary.
        If `load_metadata` is True, returns a tuple containing the data dictionary
        and metadata dictionary.

    """

    if directory is None:
        directory = os.path.join(DATA_DIR, "tb_params")
    filename = "_".join([source, tb_model_name])
    filepath = os.path.join(directory, f"{filename}.json")
    data = load_json(filepath)
    if load_metadata:
        return data["data"], data["metadata"]
    return data["data"]


def delete_tb_params(
    source,
    tb_model_name,
    directory=None,
):
    """
    Deletes a tight-binding parameter file.

    Parameters
    ----------
    source : str
        The source identifier for the parameter file.
    tb_model_name : str
        The name of the tight-binding model.
    directory : str, optional
        The directory containing the parameter files. Defaults to
        a subdirectory 'tb_params' within DATA_DIR.

    Returns
    -------
    bool
        True if the file was successfully deleted, False otherwise.
    """

    if directory is None:
        directory = os.path.join(DATA_DIR, "tb_params")
    filename = "_".join([source, tb_model_name])
    filepath = os.path.join(directory, f"{filename}.json")

    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    return False


# -----------------------------------------------------------------------
