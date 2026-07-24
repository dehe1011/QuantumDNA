from ..model import TBModel

# ----------------------------------------------------------------------


def get_tb_sites(upper_strand, **kwargs):
    """
    Generate tight-binding (TB) sites for a DNA sequence based on the given parameters.

    Parameters
    ----------
    upper_strand : list
        The upper strand of the DNA sequence.

    Returns
    -------
    list of lists
        TIght-binding sites

    Notes
    -----
    .. note::

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

    tb_model = TBModel(len(upper_strand), **kwargs)
    is_double_stranded = tb_model.double_stranded
    is_pdb = len(upper_strand[0]) == 3
    lower_strand = kwargs.get("lower_strand", "auto complete")
    backbone = tb_model.backbone

    complementary_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # double stranded and not PDB input
    if is_double_stranded and not is_pdb:
        if lower_strand == "auto complete":
            lower_strand = [complementary_dict[key] for key in upper_strand]
        tb_sites = [list(upper_strand), list(lower_strand)]

        if backbone:
            B = ["B"] * len(upper_strand)
            tb_sites = [B, *tb_sites, B]

    # single stranded and not PDB input
    if not is_double_stranded and not is_pdb:
        tb_sites = [list(upper_strand)]

        if backbone:
            B = ["B"] * len(upper_strand)
            tb_sites = [B, *tb_sites, B]

    # double stranded and PDB input
    if is_double_stranded and is_pdb:
        assert lower_strand != "auto complete", "Please provide a lower strand for PDB files."
        tb_sites = [list(upper_strand), list(lower_strand)]

        if backbone:
            B_upper = [f"{e[:2]}B" for e in upper_strand]
            B_lower = [f"{e[:2]}B" for e in lower_strand]
            tb_sites = [B_upper, *tb_sites, B_lower]

    # single stranded and PDB input
    if not is_double_stranded and is_pdb:
        tb_sites = [list(upper_strand)]

        if backbone:
            n = len(upper_strand)
            B_upper = [f"{e[:2]}B" for e in upper_strand]
            B_lower = [f"{2*n+1-int(e[:2])}".zfill(2) + "B" for e in upper_strand]
            tb_sites = [B_upper, *tb_sites, B_lower]

    return tb_sites


# ----------------------------------------------------------------------
