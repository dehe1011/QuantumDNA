import yaml
from .helpers import get_non_overwriting_path

# ----------------------------------------------------------------------


def save_yaml(filepath, data):
    """
    Save data to a YAML file at the specified filepath, ensuring no overwriting.
    """

    filepath = get_non_overwriting_path(filepath)
    with open(filepath, "w", encoding="utf-8") as file:
        yaml.safe_dump(data, file, indent=4, default_flow_style=False)


def load_yaml(filepath):
    """
    Load and parse a YAML file, returning its contents as a Python object.
    """

    with open(filepath, "r", encoding="utf-8") as file:
        data = yaml.safe_load(file)
        return data


def modify_yaml(filepath, key, value, override=False):
    """
    Modify a YAML file by updating or appending a value to a specified key.
    """

    data = load_yaml(filepath)
    if data is not None:
        if override:
            data[key] = value
        else:
            if not value in data[key]:
                data[key].append(value)
        save_yaml(filepath, data)


# ----------------------------------------------------------------------
