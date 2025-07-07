# ----------------------------------------------------------------------


# pylint: disable=too-many-statements
def get_tb_couplings(tb_model_name, num_sites_per_strand):
    n = num_sites_per_strand

    if tb_model_name == "WM":
        t_upper = [["t", (0, i), (0, i + 1)] for i in range(n - 1)]
        couplings = t_upper

    elif tb_model_name == "LM":
        t_upper = [["t", (0, i), (0, i + 1)] for i in range(n - 1)]
        t_lower = [["t", (1, i + 1), (1, i)] for i in range(n - 1)]
        h = [["h", (0, i), (1, i)] for i in range(n)]
        couplings = t_upper + t_lower + h

    elif tb_model_name == "ELM":
        t_upper = [["t", (0, i), (0, i + 1)] for i in range(n - 1)]
        r_plus = [["r+", (0, i), (1, i + 1)] for i in range(n - 1)]
        r_minus = [["r-", (1, i), (0, i + 1)] for i in range(n - 1)]
        t_lower = [["t", (1, i + 1), (1, i)] for i in range(n - 1)]
        h = [["h", (0, i), (1, i)] for i in range(n)]
        couplings = t_upper + r_plus + r_minus + t_lower + h

    elif tb_model_name == "FWM":
        t = [["t", (1, i), (1, i + 1)] for i in range(n - 1)]
        h1 = [["h", (0, i), (1, i)] for i in range(n)]
        h2 = [["h", (1, i), (2, i)] for i in range(n)]
        couplings = t + h1 + h2

    elif tb_model_name == "FLM":
        t_upper = [["t", (1, i), (1, i + 1)] for i in range(n - 1)]
        t_lower = [["t", (2, i + 1), (2, i)] for i in range(n - 1)]
        h = [["h", (1, i), (2, i)] for i in range(n)]
        h1 = [["h", (0, i), (1, i)] for i in range(n)]
        h2 = [["h", (2, i), (3, i)] for i in range(n)]
        couplings = t_upper + t_lower + h + h1 + h2

    elif tb_model_name == "FELM":
        t_upper = [["t", (1, i), (1, i + 1)] for i in range(n - 1)]
        r_plus = [["r+", (1, i), (2, i + 1)] for i in range(n - 1)]
        r_minus = [["r-", (2, i), (1, i + 1)] for i in range(n - 1)]
        t_lower = [["t", (2, i + 1), (2, i)] for i in range(n - 1)]
        h = [["h", (1, i), (2, i)] for i in range(n)]
        h1 = [["h", (0, i), (1, i)] for i in range(n)]
        h2 = [["h", (2, i), (3, i)] for i in range(n)]
        couplings = t_upper + r_plus + r_minus + t_lower + h + h1 + h2

    elif tb_model_name == "TC":
        t = [["t", (1, i), (1, i + 1)] for i in range(n - 1)]
        h1 = [["h", (0, i), (1, i)] for i in range(n)]
        h2 = [["h", (1, i), (2, i)] for i in range(n)]
        t1 = [["t", (0, i), (0, i + 1)] for i in range(n - 1)]
        t2 = [["t", (2, i), (2, i + 1)] for i in range(n - 1)]
        couplings = t + h1 + h2 + t1 + t2

    elif tb_model_name == "FC":
        t_upper = [["t", (1, i), (1, i + 1)] for i in range(n - 1)]
        r_plus = [["r+", (1, i), (2, i + 1)] for i in range(n - 1)]
        r_minus = [["r-", (2, i), (1, i + 1)] for i in range(n - 1)]
        t_lower = [["t", (2, i + 1), (2, i)] for i in range(n - 1)]
        h = [["h", (1, i), (2, i)] for i in range(n)]
        h1 = [["h", (0, i), (1, i)] for i in range(n)]
        h2 = [["h", (2, i), (3, i)] for i in range(n)]
        t1 = [["t", (0, i), (0, i + 1)] for i in range(n - 1)]
        t2 = [["t", (3, i), (3, i + 1)] for i in range(n - 1)]
        couplings = t_upper + r_plus + r_minus + t_lower + h + h1 + h2 + t1 + t2
    else:
        couplings = []
        raise ValueError(f"Unknown tight-binding model: {tb_model_name}")
    return couplings


def get_tb_energies(num_channels, num_sites_per_strand):
    energies = []
    for i in range(num_channels):
        for j in range(num_sites_per_strand):
            energies.append(["E", (i, j), (i, j)])
    return energies


def get_tb_config(tb_model_name, tb_dims):
    num_channels, num_sites_per_strand = tb_dims
    energies = get_tb_energies(num_channels, num_sites_per_strand)
    couplings = get_tb_couplings(tb_model_name, num_sites_per_strand)
    return energies + couplings


# ----------------------------------------------------------------------
