from .save_load import *
from .helpers import *

CONFIGS: dict = load_json("config", os.path.join(DATA_DIR, "raw"), load_metadata=False)
CONFIG: dict = {**DEFAULTS, **CONFIGS}
TB_MODELS_PROPS: dict = load_json("tb_models", os.path.join(DATA_DIR, "raw"))

DNA_BASES: list = CONFIGS["DNA_BASES"]
TB_MODELS: list = CONFIGS["TB_MODELS"]
SOURCES: list = CONFIGS["SOURCES"]
DESCRIPTIONS: list = CONFIGS["DESCRIPTIONS"]
PARTICLES: list = CONFIGS["PARTICLES"]
UNITS: list = CONFIGS["UNITS"]
T_UNITS: list = CONFIGS["T_UNITS"]
SPECTRAL_DENSITIES: list = CONFIGS["SPECTRAL_DENSITIES"]

from .check_input import *
