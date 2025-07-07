import json
from .helpers import get_non_overwriting_path

# ----------------------------------------------------------------------


def save_json(data, filepath, override=False):

    if not override:
        filepath = get_non_overwriting_path(filepath)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)


def load_json(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data


def modify_json(filepath, key, value, override=False):
    data = load_json(filepath)
    if data is not None:
        if override:
            data[key] = value
        else:
            if not value in data[key]:
                data[key].append(value)
        save_json(data, filepath, override=True)


# ----------------------------------------------------------------------
