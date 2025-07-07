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
    """Deletes the tight-binding parameters file."""

    if directory is None:
        directory = os.path.join(DATA_DIR, "tb_params")
    filename = "_".join([source, tb_model_name])
    filepath = os.path.join(directory, f"{filename}.json")

    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    return False


# -----------------------------------------------------------------------
