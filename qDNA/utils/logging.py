import os
import logging
from .. import DATA_DIR

# ----------------------------------------------------------------------


def setup_logging(filepath=None):
    """Sets up logging for the application."""

    if filepath is None:
        filepath = os.path.join(DATA_DIR, "qDNA.log")

    logging.basicConfig(
        filename=filepath,
        level=logging.INFO,
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
    )


setup_logging()
# use logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
