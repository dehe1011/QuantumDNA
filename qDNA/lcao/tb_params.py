from .save_load import load_xyz, find_xyz_files
from .base import Base
from .base_pair import BasePair
from .dimer import Dimer

# ----------------------------------------------------


def calc_tb_energies_monomers(directory):
    """
    Calculate the tight-binding energies for monomers in a given directory.
    This function searches for XYZ files in the specified directory, loads the data,
    and calculates the HOMO and LUMO energies for each monomer. The results are
    returned as dictionaries.

    Parameters
    ----------
    directory : str
        The directory containing the XYZ files.

    Returns
    -------
    HOMO_dict : dict
        A dictionary where the keys are identifiers of the form 'E_<base_identifier>'
        and the values are the HOMO energies in eV rounded to two decimal places.
    LUMO_dict : dict
        A dictionary where the keys are identifiers of the form 'E_<base_identifier>'
        and the values are the LUMO energies in eV rounded to two decimal places.
    """

    filenames = find_xyz_files(directory)  # e.g., ['A1', 'A2', 'T3', 'T4']

    HOMO_dict, LUMO_dict = {}, {}

    for filename in filenames:
        xyz_identifier, xyz_data = load_xyz(filename, directory)
        base = Base(xyz_identifier, xyz_data)

        HOMO_dict["E_" + base.identifier] = round(float(base.E_HOMO), 2)
        LUMO_dict["E_" + base.identifier] = round(float(base.E_LUMO), 2)

    return HOMO_dict, LUMO_dict


def calc_tb_params_dimer(bases, tb_model_name, double_stranded=True):
    """
    Calculate tight-binding (TB) parameters for a dimer in meV.

    Parameters
    ----------
    bases : list
        List of base objects representing the bases in the dimer.
    tb_model_name : str
        Name of the tight-binding model. Can be "WM", "LM", or "ELM".
    double_stranded : bool, optional
        Flag indicating whether the dimer is double-stranded (default is True).

    Returns
    -------
    tuple
        A tuple containing two dictionaries:
        - HOMO_dict: Dictionary with HOMO (Highest Occupied Molecular Orbital) parameters.
        - LUMO_dict: Dictionary with LUMO (Lowest Unoccupied Molecular Orbital) parameters.

    Raises
    ------
    AssertionError
        If `tb_model_name` is not "WM" for single-stranded mode.
        If the length of `bases` is not 2 for single-stranded mode.
    """

    HOMO_dict, LUMO_dict = {}, {}

    if not double_stranded:
        assert (
            tb_model_name == "WM"
        ), "For single-stranded mode, the TB model must be 'WM'."
        assert len(bases) == 2, "Single-stranded dimers should only contain two bases. "

        dimer = Dimer(bases[0], bases[1], identifier=None)
        HOMO_dict["t_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["t_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        for base in bases:
            HOMO_dict["E_" + base.identifier] = round(float(base.E_HOMO) * 1e3, 2)
            LUMO_dict["E_" + base.identifier] = round(float(base.E_LUMO) * 1e3, 2)

        return HOMO_dict, LUMO_dict

    if tb_model_name in ["LM", "ELM"]:
        dimer = Dimer(bases[0], bases[1], identifier=None)
        HOMO_dict["t_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["t_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        dimer = Dimer(bases[0], bases[3], identifier=None)
        HOMO_dict["h_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["h_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        dimer = Dimer(bases[1], bases[2], identifier=None)
        HOMO_dict["h_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["h_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        dimer = Dimer(bases[2], bases[3], identifier=None)
        HOMO_dict["t_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["t_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        for base in bases:
            HOMO_dict["E_" + base.identifier] = round(float(base.E_HOMO) * 1e3, 2)
            LUMO_dict["E_" + base.identifier] = round(float(base.E_LUMO) * 1e3, 2)

    if tb_model_name == "ELM":
        dimer = Dimer(bases[0], bases[2], identifier=None)
        HOMO_dict["r+_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["r+_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        dimer = Dimer(bases[1], bases[3], identifier=None)
        HOMO_dict["r-_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["r-_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

    if tb_model_name == "WM":
        basepairs = [BasePair(bases[0], bases[3]), BasePair(bases[1], bases[2])]

        dimer = Dimer(basepairs[0], basepairs[1], identifier=None)
        HOMO_dict["t_" + dimer.identifier] = round(float(dimer.t_HOMO) * 1e3, 2)
        LUMO_dict["t_" + dimer.identifier] = round(float(dimer.t_LUMO) * 1e3, 2)

        for base in basepairs:
            HOMO_dict["E_" + base.identifier] = round(float(base.E_HOMO) * 1e3, 2)
            LUMO_dict["E_" + base.identifier] = round(float(base.E_LUMO) * 1e3, 2)

    return HOMO_dict, LUMO_dict


def calc_tb_params(directories, tb_model_name, double_stranded=True):
    """
    Calculate tight-binding parameters in meV for a set of DNA base pairs.

    Parameters
    ----------
    directories : list of str
        List of directories containing the XYZ files for the DNA bases.
    tb_model_name : str
        Name of the tight-binding model to be used.
    double_stranded : bool, optional
        If True, consider the DNA as double-stranded. If False, consider it as single-stranded. Default is True.

    Returns
    -------
    HOMO_dict : dict
        Dictionary containing the Highest Occupied Molecular Orbital (HOMO) parameters for each dimer.
    LUMO_dict : dict
        Dictionary containing the Lowest Unoccupied Molecular Orbital (LUMO) parameters for each dimer.
    """

    HOMO_dict, LUMO_dict = {}, {}
    for directory in directories:
        filenames = find_xyz_files(directory)  # e.g., ['A1', 'A2', 'T3', 'T4']
        num_bases = len(filenames)

        bases = []
        for filename in filenames:
            xyz_identifier, xyz_data = load_xyz(filename, directory)
            bases.append(Base(xyz_identifier, xyz_data))

        dimers = None
        if not double_stranded:
            dimers = [[bases[i], bases[i + 1]] for i in range(num_bases - 1)]

        if double_stranded:
            dimers = [
                [
                    bases[i],
                    bases[i + 1],
                    bases[num_bases - 2 - i],
                    bases[num_bases - i - 1],
                ]
                for i in range(num_bases // 2 - 1)
            ]
        assert (
            dimers is not None
        ), "No dimers found. Check the input directories and files."

        for dimer in dimers:
            HOMO_dict_new, LUMO_dict_new = calc_tb_params_dimer(
                dimer, tb_model_name, double_stranded=double_stranded
            )
            HOMO_dict.update(HOMO_dict_new)
            LUMO_dict.update(LUMO_dict_new)

    return HOMO_dict, LUMO_dict
