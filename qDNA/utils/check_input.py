from ..io import OPTIONS

# ----------------------------------------------------------------------


def check_lcao_kwargs(**lcao_kwargs):
    """
    Validate and check the keyword arguments for LCAO parameters.
    """

    # check for None values
    assert not None in lcao_kwargs.values(), "None is not allowed as value"

    # check datatypes
    kwargs = lcao_kwargs
    string_keys = ["param_id"]

    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    # Validate
    param_id = kwargs.get("param_id")
    param_id_opts = OPTIONS["param_ids"]

    # Validate tb_model_name
    assert (
        param_id in param_id_opts
    ), f"Invalid TB model name: '{param_id}'. Must be one of {param_id_opts}"


# ----------------------------------------------------------------------


def check_tb_model_kwargs(**tb_model_kwargs):
    """
    Validates the keyword arguments for the TBModel class.
    """

    # check for None values
    assert not None in tb_model_kwargs.values(), "None is not allowed as value"

    # check datatypes
    kwargs = tb_model_kwargs
    string_keys = ["tb_model_name"]

    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    # Validate
    tb_model_name = kwargs.get("tb_model_name")
    tb_model_name_opts = OPTIONS["tb_models"]

    # Validate tb_model_name
    assert (
        tb_model_name in tb_model_name_opts
    ), f"Invalid TB model name: '{tb_model_name}'. Must be one of {tb_model_name_opts}"


# ----------------------------------------------------------------------


def check_tb_ham_kwargs(**ham_kwargs):
    """
    Validates the keyword arguments for the TBHam class.
    """

    # check for None values
    assert not None in ham_kwargs.values(), "None is not allowed as value"

    # check datatypes
    kwargs = ham_kwargs
    string_keys = ["source", "description", "unit"]
    float_keys = ["coulomb_param", "exchange_param"]
    bool_keys = ["relaxation", "nn_cutoff"]

    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"
    for key in bool_keys:
        assert isinstance(kwargs.get(key), bool), f"{key} must be of type bool"

    # Validate
    particles = kwargs.get("particles")
    description = kwargs.get("description")
    unit = kwargs.get("unit")
    # source = kwargs.get("source")  # Uncomment if needed
    particles_opts = OPTIONS["particles"]
    description_opts = OPTIONS["descriptions"]
    unit_opts = OPTIONS["units"]
    # sources_opts = CONFIG["SOURCES"]

    # Validate particles
    assert isinstance(
        particles, list
    ), f"`particles` must be a list, got {type(particles).__name__}"
    for particle in particles:
        assert isinstance(
            particle, str
        ), f"Each particle must be a string, got {type(particle).__name__}"
        assert (
            particle in particles_opts
        ), f"Invalid particle: '{particle}'. Must be one of {particles_opts}"

    # Validate description
    assert (
        description in description_opts
    ), f"Invalid description: '{description}'. Must be one of {description_opts}"

    # Description-specific particle constraints
    if description == "1P":
        assert particles in [
            ["electron"],
            ["hole"],
            ["exciton"],
        ], f"For 1P description, `particles` must be ['electron'] or ['hole'] or ['exciton'], got {particles}"

    # Validate unit
    assert unit in unit_opts, f"Invalid unit: '{unit}'. Must be one of {unit_opts}"

    # Validate source
    # assert source in sources_opts, f"Invalid source: '{source}'. Must be one of {sources_opts}"


# ----------------------------------------------------------------------


def check_lind_diss_kwargs(**diss_kwargs):
    """
    Validates the keyword arguments for the LindDiss class.
    """

    # Check for None values
    assert None not in diss_kwargs.values(), "None is not allowed as a value"

    # Check datatypes
    kwargs = diss_kwargs
    string_keys = ["spectral_density"]
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
    bool_keys = ["loc_therm", "glob_therm", "uniform_relaxation"]

    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"

    for key in float_keys:
        assert isinstance(
            kwargs.get(key), (float, int)
        ), f"{key} must be of type float or int"

    for key in bool_keys:
        assert isinstance(kwargs.get(key), bool), f"{key} must be of type bool"

    relax_rates = kwargs.get("relax_rates")
    assert isinstance(relax_rates, dict), "`relax_rates` must be a dictionary"

    # Validate
    spectral_density = kwargs.get("spectral_density")
    spectral_densities_opts = OPTIONS["spectral_densities"]
    loc_deph_rate = kwargs.get("loc_deph_rate")
    glob_deph_rate = kwargs.get("glob_deph_rate")
    loc_therm = kwargs.get("loc_therm")
    glob_therm = kwargs.get("glob_therm")

    # Validate spectral density
    assert (
        spectral_density in spectral_densities_opts
    ), f"Invalid spectral density: '{spectral_density}'. Must be one of {spectral_densities_opts}"

    # Mutual exclusions
    assert not (
        loc_deph_rate != 0 and glob_deph_rate != 0
    ), "Dephasing must either be local or global, not both"
    assert not (
        loc_therm and glob_therm
    ), "Thermalization must either be local or global, not both"


# ----------------------------------------------------------------------


def check_me_solver_kwargs(**me_kwargs):
    """
    Validates the keyword arguments for the MeSolver class.
    """

    # check for None values
    assert not None in me_kwargs.values(), "None is not allowed as value"

    # Check datatypes
    kwargs = me_kwargs
    string_keys = ["t_unit"]
    float_keys = ["t_steps", "t_end"]
    list_keys = ["init_e_states", "init_h_states", "init_ex_states"]
    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"
    for key in list_keys:
        assert isinstance(kwargs.get(key), list), f"{key} must be of type list"

    # Validate
    t_unit_opts = OPTIONS["t_units"]
    t_unit = kwargs.get("t_unit")
    # check values
    assert t_unit in t_unit_opts, f"t_unit must be in {t_unit_opts}"


# ----------------------------------------------------------------------


def check_kwargs(**kwargs):
    """
    Validates the keyword arguments.
    """

    check_tb_model_kwargs(**kwargs)
    check_tb_ham_kwargs(**kwargs)
    check_lind_diss_kwargs(**kwargs)
    check_me_solver_kwargs(**kwargs)


# ----------------------------------------------------------------------
