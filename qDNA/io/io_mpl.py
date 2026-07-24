import os
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt

from .. import ROOT_DIR, DATA_DIR
from .io_json import load_json
from .helpers import get_non_overwriting_path

# ----------------------------------------------------------------------


def load_mpl_style(filepath=None):
    """
    Load and apply a custom Matplotlib style file to the style library.
    """

    config_dir = mpl.get_configdir()
    stylelib_dir = os.path.join(config_dir, "stylelib")
    os.makedirs(stylelib_dir, exist_ok=True)
    if filepath is None:
        filepath = os.path.join(ROOT_DIR, "qDNA", "qDNA-default.mplstyle")
    shutil.copy(filepath, stylelib_dir)
    plt.style.reload_library()


def save_figure(fig, filepath):
    """
    Save a matplotlib figure to the specified filepath without overwriting existing files.
    """

    filepath = get_non_overwriting_path(filepath)
    fig.savefig(filepath, transparent=True, pad_inches=0)


def load_color_palette(color_palette_id):
    """
    Load a color palette from a JSON file based on the given palette ID.
    """

    filepath = os.path.join(DATA_DIR, "color_palettes", f"{color_palette_id}.json")
    return load_json(filepath)


def plot_color_palette(color_palette):
    """
    Visualizes a color palette as a horizontal bar of colors.
    """

    colors = color_palette
    n_colors = len(colors)

    fig, ax = plt.subplots(figsize=(n_colors, 1))

    for i, color in enumerate(colors):
        ax.add_patch(plt.Rectangle((i / n_colors, 0), 1 / n_colors, 1, color=color))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    return fig, ax


# ----------------------------------------------------------------------
