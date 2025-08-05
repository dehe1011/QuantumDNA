# pylint: disable=wrong-import-position
__version__ = "0.1.10"

import pathlib
import os

ROOT_DIR = str(pathlib.Path(__file__).absolute().parent.parent)
DATA_DIR = os.path.join(ROOT_DIR, "qDNA", "data")

from .io import *
from .utils import *

from .lcao import *
from .model import *
from .hamiltonian import *
from .environment import *
from .dynamics import *
from .evaluation import *
from .visualization import *
from .legacy import *

# some other comment...
# one more comment...
# some nicer comment...
