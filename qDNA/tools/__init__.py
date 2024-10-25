from .save_load import *
from .helpers import *

CONFIGS = load_json("config", os.path.join(DATA_DIR, "raw"), load_metadata=False)
CONFIG = {**DEFAULTS, **CONFIGS}
TB_MODELS_PROPS = load_json("tb_models", os.path.join(DATA_DIR, "raw"))

DNA_BASES = CONFIGS["DNA_BASES"]
TB_MODELS = CONFIGS["TB_MODELS"]
SOURCES = CONFIGS["SOURCES"]
DESCRIPTIONS = CONFIGS["DESCRIPTIONS"]
PARTICLES = CONFIGS["PARTICLES"]
UNITS = CONFIGS["UNITS"]
T_UNITS = CONFIGS["T_UNITS"]
SPECTRAL_DENSITIES = CONFIGS["SPECTRAL_DENSITIES"]

from .check_input import *
