import json
import logging
import os
import random
import yaml

from qDNA import ROOT_DIR

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    filename=os.path.join(ROOT_DIR, "qDNA", "data", "logging.log"),
    encoding="utf-8",
    level=logging.INFO,
)


def load_yaml(filepath):
    """
    Load a YAML file and return its contents.

    Parameters
    ----------
    filepath : str
        The path to the YAML file to be loaded.

    Returns
    -------
    dict
        The contents of the YAML file as a dictionary.
    """
    with open(filepath, "r") as yaml_file:
        return yaml.safe_load(yaml_file)


def get_config():
    """
    Retrieves the configuration settings from the 'config.yaml' file located in the 'qDNA' directory.

    Returns:
        dict: The configuration settings loaded from the YAML file.
    """
    filepath = os.path.join(ROOT_DIR, "qDNA", "config.yaml")
    return load_yaml(filepath)


CONFIG = get_config()


def save_json(data, metadata, filename, directory):
    """
    Save data and metadata to a JSON file.
    Parameters
    ----------
    data : dict
        The main data to be saved.
    metadata : dict
        Additional metadata to be saved alongside the data.
    filename : str
        The name of the file (without extension) to save the data to.
    directory : str
        The directory where the file should be saved.
    Returns
    -------
    None
    Notes
    -----
    If a file with the same name already exists in the specified directory,
    a warning will be logged and the function will return without saving.
    """
    os.makedirs(directory, exist_ok=True)

    filepath = os.path.join(directory, filename + ".json")
    if os.path.exists(filepath):
        logger.warning(f"File {filepath} already exists.")
        if CONFIG["verbose"]:
            print(f"File {filepath} already exists")
        return

    combined_data = dict(data=data, metadata=metadata)
    with open(filepath, "w") as f:
        json.dump(combined_data, f)

    logger.info(f"Data saved as {filepath}")
    if CONFIG["verbose"]:
        print(f"Data saved as {filepath}")


def load_json(filename, directory, load_metadata=False):
    """
    Load data from a JSON file.
    Parameters
    ----------
    filename : str
        The name of the JSON file (without the .json extension).
    directory : str
        The directory where the JSON file is located.
    load_metadata : bool, optional
        If True, also load metadata from the JSON file. Default is False.
    Returns
    -------
    dict or tuple
        If `load_metadata` is False, returns a dictionary containing the data.
        If `load_metadata` is True, returns a tuple containing the data and metadata.
    Raises
    ------
    FileNotFoundError
        If the specified JSON file does not exist.
    Notes
    -----
    Logs a warning if the file does not exist and an info message when the data is successfully loaded.
    """
    filepath = os.path.join(directory, filename + ".json")
    if not os.path.exists(filepath):
        logger.warning(f"File {filepath} does not exist.")
        if CONFIG["verbose"]:
            print(f"File {filepath} does not exist.")
        return

    logger.info(f"Data loaded from {filepath}")
    if CONFIG["verbose"]:
        print(f"Data loaded from {filepath}")

    with open(filepath, "r") as f:
        combined_data = json.load(f)
    if not load_metadata:
        return combined_data["data"]
    else:
        return combined_data["data"], combined_data["metadata"]


def save_figure(fig, filename, directory, extension="svg"):
    """
    Save a figure to a specified directory with a given filename and extension.
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure object to be saved.
    filename : str
        The name of the file (without extension) to save the figure as.
    directory : str
        The directory where the figure will be saved.
    extension : str, optional
        The file extension for the saved figure (default is "svg").
    Notes
    -----
    If a file with the same name already exists in the specified directory,
    a random index will be appended to the filename to avoid overwriting.
    The function logs the saving process and prints messages if the verbose
    configuration is enabled.
    Raises
    ------
    OSError
        If the directory does not exist or is not writable.
    """
    os.makedirs(directory, exist_ok=True)

    filepath = os.path.join(directory, filename + "." + extension)
    if os.path.exists(filepath):
        rand_idx = random.randint(0, 1000)
        logger.warning(f"Figure {filepath} already exists. Filepath changed.")
        if CONFIG["verbose"]:
            print(f"Figure {filepath} already exists. Filepath changed.")
        filepath = os.path.join(directory, filename, f"_{rand_idx}_", ".", extension)

    fig.savefig(filepath)
    logger.info(f"Figure saved as {filepath}")
    if CONFIG["verbose"]:
        print(f"Figure saved as {filepath}")
