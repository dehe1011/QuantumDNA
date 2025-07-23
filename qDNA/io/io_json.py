import json
from .helpers import get_non_overwriting_path

# ----------------------------------------------------------------------


def save_json(data, filepath, override=False):
    """
    Save data to a JSON file, optionally ensuring the file does not overwrite existing files.
    """

    if not override:
        filepath = get_non_overwriting_path(filepath)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)


def load_json(filepath):
    """
    Load and return JSON data from the specified file.
    """

    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data


def modify_json(filepath, key, value, override=False):
    """
    Modify a JSON file by updating or appending a value to a specified key.
    """

    data = load_json(filepath)
    if data is not None:
        if override:
            data[key] = value
        else:
            if not value in data[key]:
                data[key].append(value)
        save_json(data, filepath, override=True)


# ----------------------------------------------------------------------
