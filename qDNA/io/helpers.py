import os

# ----------------------------------------------------------------------


def get_non_overwriting_path(filepath):
    """
    Generate a unique file path by appending a counter if the given path already exists.
    """

    base, ext = os.path.splitext(filepath)  # ext: .yaml, .json, .xyz, etc.
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    return new_filepath


# ----------------------------------------------------------------------
