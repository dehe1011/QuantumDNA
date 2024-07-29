import pandas as pd
import json
import yaml
import os
import pathlib
import logging
import pickle

ROOT_DIR = str(pathlib.Path(__file__).absolute().parent.parent)

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    filename=os.path.join(ROOT_DIR, "data", "logging.log"),
    encoding="utf-8",
    level=logging.INFO,
)

filepath = os.path.join(ROOT_DIR, "config.yaml")
with open(filepath, "r") as file:
    config = yaml.safe_load(file)
verbose = config["verbose"]

# ------------------------------------------------- save and load data --------------------------------------------


def my_save(
    data,
    metadata,
    filename,
    directory=os.path.join(ROOT_DIR, "data"),
    save_excel=False,
    version_index: bool = True,
):
    """
    Saves data and metadata to a json file. The data (without metadata) can additionally be saved as excel file.

    Note:
    -----
        Usually the filepath starts with 'data/...'
    """

    os.makedirs(directory, exist_ok=True)

    if version_index:
        filename += "_version_"
        index = 0
        filepath = os.path.join(directory, filename)
        while os.path.exists(filepath + str(index) + ".json"):
            index += 1
        filename += str(index) + ".json"

    else:
        filename += ".json"

    filepath = os.path.join(directory, filename)

    logger.info(f"Data saved as {filepath}")
    if verbose:
        print(f"Data saved as {filepath}")

    data = (data, metadata)
    with open(filepath, "w") as f:
        json.dump(data, f)
    if save_excel:
        df = pd.DataFrame(data)
        df.to_excel("data.xlsx", index=False)


def my_load(filename, load_metadata=False, directory=os.path.join(ROOT_DIR, "data")):
    """
    Loads the data and the metadata (if wanted) from a json file.

    Note:
    -----
        Usually the filepath starts with 'data/...'
    """
    os.makedirs(directory, exist_ok=True)
    filename += ".json"
    filepath = os.path.join(directory, filename)

    logger.info(f"Data loaded from {filepath}")
    if verbose:
        print(f"Data loaded from {filepath}")

    with open(filepath, "r") as f:
        data = json.load(f)
    if not load_metadata:
        return data[0]
    else:
        return data[0], data[1]


def convert_pickle_to_json(filepath):
    """
    Converts a pickle file to a json file with better readablity.
    """
    with open(filepath, "rb") as f:
        data = pickle.load(f)
    with open(filepath + ".json", "w") as f:
        json.dump(data, f)


# ----------------------------------------- load configuration file with defualt values --------------------------


def get_config(filename="config", directory=ROOT_DIR):
    """
    Load data stored in a .yaml configuration file (e.g., global variables).
    """
    filename += ".yaml"
    filepath = os.path.join(directory, filename)
    with open(filepath, "r") as file:
        config = yaml.safe_load(file)
    return config


# ---------------------------------------------------- save figures ------------------------------------


def save_figure(
    fig,
    filename: str,
    directory: str = os.path.join(ROOT_DIR, "data", "figures"),
    format: str = "pdf",
):
    """
    Save Figures.
    """
    os.makedirs(directory, exist_ok=True)

    filename += "_version_"
    index = 0
    while os.path.exists(os.path.join(directory, filename + str(index) + "." + format)):
        index += 1
    filename += str(index) + "." + format
    filepath = os.path.join(directory, filename)
    fig.savefig(filepath)
    logger.info(f"Figure saved as {filepath}")
    if verbose:
        print(f"Figure saved as {filepath}")
