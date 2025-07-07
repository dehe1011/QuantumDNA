import os

# ----------------------------------------------------------------------


def get_non_overwriting_path(filepath):
    base, ext = os.path.splitext(filepath)  # ext: .yaml, .json, .xyz, etc.
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    return new_filepath


# ----------------------------------------------------------------------
