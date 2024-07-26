from .save_load import get_config


def check_ham_kwargs(**ham_kwargs):
    config = get_config()
    assert not None in ham_kwargs.values(), "None is not allowed as value"
    kwargs = ham_kwargs

    string_keys = ["source", "description", "unit"]
    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    float_keys = ["interaction_param"]
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"

    bool_keys = ["relaxation", "nn_cutoff"]
    for key in bool_keys:
        assert isinstance(kwargs.get(key), bool), f"{key} must be of type bool"

    assert isinstance(kwargs.get("particles"), list), "particles must be of type list"

    assert all(
        isinstance(particle, str) for particle in kwargs["particles"]
    ), "elements of particles must be of type str"
    assert all(
        [particle in config["PARTICLES"] for particle in kwargs["particles"]]
    ), f"all particles must be in {config['PARTICLES']}"
    if kwargs["description"] == "1P":
        assert kwargs["particles"] in [
            ["electron"],
            ["hole"],
        ], "in the one particle description the particles must be either ['electron'] or ['hole']"
    assert (
        kwargs["source"] in config["SOURCES"]
    ), f"source must be in {config['SOURCES']}"
    assert (
        kwargs["description"] in config["DESCRIPTIONS"]
    ), f"description must be in {config['DESCRIPTIONS']}"
    assert kwargs["unit"] in config["UNITS"], f"unit must be in {config['UNITS']}"


# --------------------------------------------------


def check_diss_kwargs(**diss_kwargs):
    config = get_config()
    assert not None in diss_kwargs.values(), "None is not allowed as value"
    kwargs = diss_kwargs

    string_keys = ["spectral_density"]
    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    float_keys = [
        "loc_deph_rate",
        "glob_deph_rate",
        "deph_rate",
        "relax_rate",
        "cutoff_freq",
        "reorg_energy",
        "temperature",
        "exponent",
    ]
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"

    bool_keys = ["loc_therm", "glob_therm", "uniform_relaxation"]
    for key in bool_keys:
        assert isinstance(kwargs.get(key), bool), f"{key} must be of type bool"

    assert isinstance(kwargs["relax_rates"], dict), f"relax_rates must be of form dict"
    DNA_SITES = config["DNA_BASES"] + ["B"]
    assert all(
        [
            key in DNA_SITES and isinstance(value, (int, float))
            for key, value in kwargs["relax_rates"].items()
        ]
    )
    assert (
        kwargs["spectral_density"] in config["SPECTRAL_DENSITIES"]
    ), f"spectral_density must be in {config['SPECTRAL_DENSITIES']}"
    assert not (
        kwargs["loc_deph_rate"] != 0 and kwargs["glob_deph_rate"] != 0
    ), "dephasing must either be local or global"
    assert not (
        kwargs["loc_therm"] and kwargs["glob_therm"]
    ), "thermalizing must either be local or global"


# --------------------------------------------------


def check_me_kwargs(**me_kwargs):
    config = get_config()
    assert not None in me_kwargs.values(), "None is not allowed as value"
    kwargs = me_kwargs

    string_keys = ["init_e_state", "init_h_state", "t_unit"]
    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    float_keys = ["t_steps", "t_end"]
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"

    assert (
        kwargs["t_unit"] in config["T_UNITS"]
    ), f"t_unit must be in {config['T_UNITS']}"
