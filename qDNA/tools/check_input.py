"""This module provides functions to validate keyword arguments for various
configurations related to Hamiltonian, dissipation parameters, and other miscellaneous
parameters."""

from . import CONFIG

# --------------------------------------------------


def check_ham_kwargs(**ham_kwargs):
    """Validates the keyword arguments for Hamiltonian configuration.

    Parameters
    ----------
    **ham_kwargs : dict
        Arbitrary keyword arguments representing Hamiltonian configuration. Expected keys are:
        - 'source' (str): The source of the Hamiltonian.
        - 'description' (str): Description of the Hamiltonian.
        - 'unit' (str): Unit of measurement.
        - 'interaction_param' (float or int): Interaction parameter.
        - 'relaxation' (bool): Relaxation flag.
        - 'nn_cutoff' (bool): Nearest neighbor cutoff flag.
        - 'particles' (list of str): List of particles involved.
    Raises
    ------
    AssertionError
        If any of the following conditions are not met:
        - None is not allowed as a value.
        - 'source', 'description', and 'unit' must be of type str.
        - 'interaction_param' must be of type float or int.
        - 'relaxation' and 'nn_cutoff' must be of type bool.
        - 'particles' must be a list of strings.
        - All elements of 'particles' must be in CONFIG["PARTICLES"].
        - If 'description' is "1P", 'particles' must be either ["electron"] or ["hole"].
        - 'source' must be in CONFIG["SOURCES"].
        - 'description' must be in CONFIG["DESCRIPTIONS"].
        - 'unit' must be in CONFIG["UNITS"].
    """

    # check for None values
    assert not None in ham_kwargs.values(), "None is not allowed as value"

    # check datatypes
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

    # check values
    assert all(
        particle in CONFIG["PARTICLES"] for particle in kwargs["particles"]
    ), f"all particles must be in {CONFIG['PARTICLES']}"
    if kwargs["description"] == "1P":
        assert kwargs["particles"] in [
            ["electron"],
            ["hole"],
        ], "in the one particle description the particles must be either ['electron'] or ['hole']"
    assert (
        kwargs["source"] in CONFIG["SOURCES"]
    ), f"source must be in {CONFIG['SOURCES']}"
    assert (
        kwargs["description"] in CONFIG["DESCRIPTIONS"]
    ), f"description must be in {CONFIG['DESCRIPTIONS']}"
    assert kwargs["unit"] in CONFIG["UNITS"], f"unit must be in {CONFIG['UNITS']}"


# --------------------------------------------------


def check_diss_kwargs(**diss_kwargs):
    """Validates the keyword arguments for dissipation parameters.

    Parameters
    ----------
    **diss_kwargs : dict
        Keyword arguments representing dissipation parameters. The expected keys and their types are:
        - spectral_density : str
        - loc_deph_rate : float or int
        - glob_deph_rate : float or int
        - deph_rate : float or int
        - relax_rate : float or int
        - cutoff_freq : float or int
        - reorg_energy : float or int
        - temperature : float or int
        - exponent : float or int
        - loc_therm : bool
        - glob_therm : bool
        - uniform_relaxation : bool
        - relax_rates : dict
            A dictionary where keys are DNA sites and values are either int or float.
    Raises
    ------
    AssertionError
        If any of the following conditions are not met:
        - No value in `diss_kwargs` is None.
        - `spectral_density` is of type str.
        - All float keys are of type float or int.
        - All bool keys are of type bool.
        - `relax_rates` is a dictionary with keys in DNA_SITES and values of type int or float.
        - `spectral_density` is in CONFIG["SPECTRAL_DENSITIES"].
        - Either `loc_deph_rate` or `glob_deph_rate` is zero.
        - Either `loc_therm` or `glob_therm` is False.
    """

    # check for None values
    assert not None in diss_kwargs.values(), "None is not allowed as value"

    # check datatypes
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
    assert isinstance(kwargs["relax_rates"], dict), "relax_rates must be of form dict"

    # check values
    assert (
        kwargs["spectral_density"] in CONFIG["SPECTRAL_DENSITIES"]
    ), f"spectral_density must be in {CONFIG['SPECTRAL_DENSITIES']}"
    assert not (
        kwargs["loc_deph_rate"] != 0 and kwargs["glob_deph_rate"] != 0
    ), "dephasing must either be local or global"
    assert not (
        kwargs["loc_therm"] and kwargs["glob_therm"]
    ), "thermalizing must either be local or global"


# --------------------------------------------------


def check_me_kwargs(**me_kwargs):
    """Validates the keyword arguments provided to ensure they meet the required
    criteria.

    Parameters
    ----------
    **me_kwargs : dict
        Arbitrary keyword arguments to be validated. Expected keys and their types are:
        - "init_e_state" (str): Initial electronic state.
        - "init_h_state" (str): Initial hole state.
        - "t_unit" (str): Time unit, must be one of the values in CONFIG["T_UNITS"].
        - "t_steps" (float or int): Number of time steps.
        - "t_end" (float or int): End time.
    Raises
    ------
    AssertionError
        If any of the following conditions are not met:
        - None is not allowed as a value.
        - "init_e_state", "init_h_state", and "t_unit" must be of type str.
        - "t_steps" and "t_end" must be of type float or int.
        - "t_unit" must be in CONFIG["T_UNITS"].
    """

    # check for None values
    assert not None in me_kwargs.values(), "None is not allowed as value"

    # check datatypes
    kwargs = me_kwargs
    string_keys = ["init_e_state", "init_h_state", "t_unit"]
    for key in string_keys:
        assert isinstance(kwargs.get(key), str), f"{key} must be of type str"
    float_keys = ["t_steps", "t_end"]
    for key in float_keys:
        assert isinstance(kwargs.get(key), (float, int)), f"{key} must be of type float"

    # check values
    assert (
        kwargs["t_unit"] in CONFIG["T_UNITS"]
    ), f"t_unit must be in {CONFIG['T_UNITS']}"
