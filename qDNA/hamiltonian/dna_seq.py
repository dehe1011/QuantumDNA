from ..model import TBModel

# ----------------------------------------------------------------------


def get_tb_sites(upper_strand, **kwargs):
    """
    Generate tight-binding (TB) sites for a DNA sequence based on the given parameters.

    Parameters
    ----------
    upper_strand : list
        The upper strand of the DNA sequence.
    **kwargs : dict
        Additional keyword arguments:
        - lower_strand : list or str, optional
            The lower strand of the DNA sequence. If 'auto complete', the complementary
            strand is generated automatically based on the upper strand. Default is 'auto complete'.
        - Other parameters are passed to the TBModel constructor.

    Returns
    -------
    tb_sites : list of lists
        A list of lists representing the TB sites. The structure depends on the input parameters:
        - If `double_stranded` is True and `is_pdb` is False:
            Includes the upper strand, lower strand (auto-completed if not provided), and optionally
            backbone sites.
        - If `double_stranded` is False and `is_pdb` is False:
            Includes only the upper strand and optionally backbone sites.
        - If `double_stranded` is True and `is_pdb` is True:
            Includes the upper strand, lower strand (must be provided), and optionally backbone sites
            with PDB formatting.
        - If `double_stranded` is False and `is_pdb` is True:
            Includes only the upper strand and optionally backbone sites with PDB formatting.

    Notes
    -----
    - The `TBModel` class is used to determine the configuration of the TB sites.
    - The `complementary_dict` is used to auto-complete the lower strand for non-PDB inputs.
    - Backbone sites are added if required by the model.

    Examples
    --------
    >>> upper_strand = ['A', 'T', 'G', 'C']
    >>> tb_sites = get_tb_sites(upper_strand, tb_model_name = 'FELM')
    >>> print(tb_sites)
    [['B', 'B', 'B', 'B'], ['A', 'T', 'G', 'C'], ['T', 'A', 'C', 'G'], ['B', 'B', 'B', 'B']]
    """

    lower_strand = kwargs.get("lower_strand", "auto complete")
    tb_model = TBModel(len(upper_strand), **kwargs)
    is_pdb = len(upper_strand[0]) == 3
    complementary_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    if tb_model.double_stranded and not is_pdb:
        if lower_strand == "auto complete":
            lower_strand = [complementary_dict[key] for key in upper_strand]
        tb_sites = [list(upper_strand), list(lower_strand)]

        if tb_model.backbone:
            B = ["B"] * len(upper_strand)
            tb_sites = [B, *tb_sites, B]

    if not tb_model.double_stranded and not is_pdb:
        tb_sites = [list(upper_strand)]

        if tb_model.backbone:
            B = ["B"] * len(upper_strand)
            tb_sites = [B, *tb_sites, B]

    if tb_model.double_stranded and is_pdb:
        assert (
            lower_strand != "auto complete"
        ), "Please provide a lower strand for PDB files."
        tb_sites = [list(upper_strand), list(lower_strand)]

        if tb_model.backbone:
            B_upper = [f"{e[:2]}B" for e in upper_strand]
            B_lower = [f"{e[:2]}B" for e in lower_strand]
            tb_sites = [B_upper, *tb_sites, B_lower]

    if not tb_model.double_stranded and is_pdb:
        tb_sites = [list(upper_strand)]

        if tb_model.backbone:
            n = len(upper_strand)
            B_upper = [f"{e[:2]}B" for e in upper_strand]
            B_lower = [f"{2*n+1-int(e[:2])}".zfill(2) + "B" for e in upper_strand]
            tb_sites = [B_upper, *tb_sites, B_lower]

    return tb_sites


# ----------------------------------------------------------------------
