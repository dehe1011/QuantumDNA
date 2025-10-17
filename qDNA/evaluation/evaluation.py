import time
import multiprocessing
import os

from tqdm import tqdm
import numpy as np

from ..lcao import Oligomer
from ..dynamics import MeSolver
from ..io import DEFAULTS, save_json, find_pdb, delete_tb_params
from ..model import get_eh_distance
from ..utils import convert_to_debye

# ----------------------------------------------------------------------


class Evaluation(MeSolver):
    """
    Evaluation class for analyzing exciton dynamics in quantum DNA models.

    Attributes
    ----------
    lifetime : float or None
        Estimated exciton lifetime in femtoseconds.
    charge_separation : float or None
        Average charge separation based on electron-hole distances.
    dipole_moment : float or None
        Dipole moment of the system in Debye units.

    Methods
    -------
    calc_lifetime()
        Calculate the exciton lifetime based on ground state population.
    calc_charge_separation(average=True)
        Compute the charge separation based on electron-hole distances.
    calc_dipole_moment()
        Calculate the dipole moment of the system.
    calc_backbone_transfer()
        Compute average transfer for backbone sites in Fishbone models.
    calc_exciton_transfer()
        Calculate exciton transfer for double-stranded, non-backbone models.

    Parameters
    ----------
    tb_sites : list
        List of tight-binding sites defining the system.
    **kwargs : dict
        Additional parameters for the MeSolver superclass.
    """

    def __init__(self, tb_sites, **kwargs):
        super().__init__(tb_sites, **kwargs)
        self.lifetime = None
        self.charge_separation = None
        self.dipole_moment = None

    def calc_lifetime(self):
        """
        Calculate the estimated exciton lifetime (in fs) based on ground state population.

        Returns
        -------
        float or str
            The calculated lifetime in the specified time unit, or a message
            indicating no relaxation occurred within the given time.

        Examples
        --------
        >>> eva = Evaluation([list('GCG'), list('CGC')], relax_rate=3., unit="rad/ps")
        >>> eva.calc_lifetime()
        775.5511022044088
        """

        # Check if lifetime has already been calculated
        if self.lifetime is not None:
            return self.lifetime

        start_time = time.time()
        gs_pop = self.get_groundstate_pop()["groundstate"]
        try:
            _, index = next((val, i) for i, val in enumerate(gs_pop) if val >= 1 - 1 / np.e)
            self.lifetime = self.times[index]
            if self.t_unit == "ps":
                self.lifetime *= 1000
            end_time = time.time()
            if DEFAULTS["verbose"]:
                print(f"Calculation time: {end_time-start_time}")
            return self.lifetime
        except StopIteration:
            return "no relaxation in the given time"

    def _back_to_unitary(self):
        """Sets the simulation time to the calculated lifetime and sets relaxation rates to zero."""

        # calculate the lifetime in ps
        lifetime = self.calc_lifetime() / 1000

        # Set the simulation time to the calculated lifetime
        if isinstance(lifetime, float) and self.t_end != lifetime:
            self.t_end = lifetime
            self.t_steps = int(1000 * lifetime // 2 + 2)

        # Set relaxation rates to zero
        self.relax_rates = dict(zip(self.tb_sites_flattened, [0] * self.num_sites))

    def calc_charge_separation(self, average=True):
        """
        Calculate the charge separation based on electron-hole distances.

        Parameters
        ----------
        average : bool, optional
            If True, returns the average charge separation. Defaults to True.

        Returns
        -------
        float or list
            The calculated charge separation. Returns a single float if `average`
            is True, otherwise returns a list of charge separations.

        Examples
        --------
        >>> eva = Evaluation([list('GCG'), list('CGC')], relax_rate=0.)
        >>> eva.calc_charge_separation()
        2.951734389657976
        """

        # Check if charge separation has already been calculated
        if self.charge_separation is not None:
            return self.charge_separation

        # Perform calculation within lifetime timescale
        self._back_to_unitary()

        # Calculate the electron-hole distance
        distance_list = 3.4 * get_eh_distance(self.eh_basis)
        self.charge_separation = [distance_list @ dm.diag()[1:] for dm in self.get_result()]
        if average:
            self.charge_separation = np.mean(self.charge_separation).real
        return self.charge_separation

    def calc_dipole_moment(self):
        """
        Calculates the dipole moment of the system.

        Returns
        -------
        float
            The dipole moment in Debye units.

        Examples
        --------
        >>> eva = Evaluation([list('GCG'), list('CGC')], relax_rate=0.)
        >>> eva.calc_dipole_moment()
        14.177784530660903
        """

        # Check if dipole moment has already been calculated
        if self.dipole_moment is not None:
            return self.dipole_moment

        # Perform calculation within lifetime timescale
        self._back_to_unitary()

        # Calculate the dipole moment
        self.dipole_moment = convert_to_debye(self.calc_charge_separation())
        return self.dipole_moment

    def _calc_average_transfer(self, sites, average=True):
        """Calculates the average populations on the given TB sites in a given time period."""

        # Perform calculation within lifetime timescale
        self._back_to_unitary()

        # Calculate the average population for each particle on the given sites
        average_pop = dict(zip(self.particles, [0] * len(self.particles)))
        for particle in self.particles:
            site_populations = [self.get_pop()[f"{particle}_{site}"] for site in sites]
            total_pop = np.sum(site_populations, axis=0)
            if average:
                avg_pop = np.mean(total_pop)
            else:
                avg_pop = total_pop
            average_pop[particle] = avg_pop
        return average_pop

    def calc_backbone_transfer(self, average=True):
        """
        Calculates the average transfer for backbone sites in a Fishbone model.

        Returns
        -------
        float
            The average transfer value across the backbone sites.

        Raises
        ------
        AssertionError
            If the backbone attribute is not set, indicating the model is not a Fishbone model.
        """

        assert self.backbone, "Backbone population can only be calculated for Fishbone models"
        self._back_to_unitary()

        upper_backbone_sites, lower_backbone_sites = [], []
        for i in range(self.num_sites_per_strand):
            upper_backbone_sites.append(f"(0, {i})")
            lower_backbone_sites.append(f"({self.num_channels-1}, {i})")
        backbone_sites = upper_backbone_sites + lower_backbone_sites
        return self._calc_average_transfer(backbone_sites, average=average)

    def calc_exciton_transfer(self, average=True):
        """
        Calculate exciton transfer for a double-stranded, non-backbone model.
        This method computes the average exciton transfer between the upper and
        lower strands of the model.

        Returns
        -------
        dict
            A dictionary containing the average exciton transfer populations for
            the upper and lower strands:

            - `upper_strand_pop`: Population transfer for the upper strand.
            - `lower_strand_pop`: Population transfer for the lower strand.

        Raises
        ------
        AssertionError
            If the model is not double-stranded or if it includes a backbone.

        Examples
        --------
        >>> eva = Evaluation([list('GCG'), list('CGC')], relax_rate=0.)
        >>> eva.calc_exciton_transfer()
        ({'electron': 0.6867427675114343,
        'hole': 0.9943813264087192,
        'exciton': 0.45001514054414693},
        {'electron': 0.31325723248857257,
        'hole': 0.005618673591287626,
        'exciton': 0.0001960836601784245})
        """

        assert (
            self.double_stranded
        ), "Exciton transfer can only be calculated for double-stranded models"
        assert not self.backbone, "Exciton transfer can only be calculated for non-backbone models"
        self._back_to_unitary()

        upper_sites, lower_sites = [], []
        for i in range(self.num_sites_per_strand):
            upper_sites.append(f"(0, {i})")
            lower_sites.append(f"(1, {i})")

        upper_strand_pop = self._calc_average_transfer(upper_sites, average=average)
        lower_strand_pop = self._calc_average_transfer(lower_sites, average=average)
        return {
            "upper_strand_pop": upper_strand_pop,
            "lower_strand_pop": lower_strand_pop,
        }


# ----------------------------------------------------------------------


def evaluate(args):
    """Helper function for parallel evaluation of observables."""

    i, observables, evaluation_list = args
    evaluation = evaluation_list[i]
    result = []
    if "lifetime" in observables:
        result.append(evaluation.calc_lifetime())
    if "charge_separation" in observables:
        result.append(evaluation.calc_charge_separation())
    if "dipole_moment" in observables:
        result.append(evaluation.calc_dipole_moment())
    if "exciton_transfer" in observables:
        result.append(evaluation.calc_exciton_transfer())
    return tuple(result)


# ----------------------------------------------------------------------


class EvaluationParallel:
    """
    A class for parallel evaluation of observables across a list of sequences using multiprocessing.

    Parameters
    ----------
    evaluation_list : list
        A list of evaluation objects, each containing a sequence ID and associated data.

    Attributes
    ----------
    num_sequences : int
        Total number of sequences to evaluate.
    available_cpus : int
        Number of CPU cores available on the system.
    num_cpu : int
        Number of CPU cores allocated for multiprocessing.
    evaluation_list : list
        List of evaluation objects to process.
    sequence_id_list : list
        List of sequence IDs extracted from the evaluation objects.
    observables : list of str
        Observables to calculate for each sequence.
    args : list of tuple
        Arguments prepared for multiprocessing, containing sequence index, observables, and evaluation list.

    Methods
    -------
    calc_results(filepath=None, save=True)
        Calculate results for the sequences using multiprocessing and optionally save them to a file.

    """

    def __init__(self, evaluation_list, **kwargs):
        self.num_sequences = len(evaluation_list)
        self.available_cpus = multiprocessing.cpu_count()
        self.num_cpu = kwargs.get("num_cpu", self.available_cpus - 1)
        self.evaluation_list = evaluation_list
        self.sequence_id_list = [ev.sequence_id for ev in self.evaluation_list]
        self.observables = kwargs.get(
            "observables", ["lifetime", "charge_separation", "dipole_moment"]
        )

        self.args = [(i, self.observables, self.evaluation_list) for i in range(self.num_sequences)]

    def calc_results(self, filepath=None, save=True):
        """
        Calculate results for a set of sequences using multiprocessing.

        Parameters
        ----------
        filepath : str, optional
            Path to save the results as a JSON file. If None, results are not saved.
        save : bool, default=True
            Whether to save the results to a file.

        Returns
        -------
        dict
            A dictionary containing the calculated results for each sequence and metadata
            including observables and calculation time.
        """

        t_start = time.time()
        with multiprocessing.Pool(processes=self.num_cpu) as pool:
            result_list = list(
                tqdm(
                    pool.imap(evaluate, self.args),
                    total=self.num_sequences,
                    desc=f"Calculating observables: {self.observables}",
                    leave=True,
                    disable=False,
                )
            )
        t_end = time.time()

        result_dict = dict(zip(self.sequence_id_list, result_list))
        result_dict["metadata"] = {
            "observables": self.observables,
            "calculation_time": t_end - t_start,
        }
        if save:
            save_json(result_dict, filepath)
        return result_dict


# ----------------------------------------------------------------------


def evaluate_pdb(directory_pdb, start_site_idx, end_site_idx, tb_model_name="ELM"):
    """
    Evaluates PDB files within a specified directory and calculates results.

    Parameters
    ----------
    directory_pdb : str
        Path to the directory containing PDB files.
    start_site_idx : int
        Starting index for slicing the sites.
    end_site_idx : int
        Ending index for slicing the sites.
    tb_model_name : str, optional
        Name of the tight-binding model (default is 'ELM').

    Returns
    -------
    None
        Results are saved to a JSON file in the specified directory.
    """

    filenames = find_pdb(directory_pdb)
    evaluation_list = []
    for filename in filenames:
        oligomer = Oligomer(os.path.join(directory_pdb, filename + ".pdb"))
        oligomer.save_tb_params()
        tb_sites = oligomer.sites_id
        idx_slice = slice(start_site_idx, end_site_idx)

        evaluation = Evaluation(tb_sites[:, idx_slice], source=filename, relax_rate=3.0)
        evaluation_list.append(evaluation)

    parallel = EvaluationParallel(evaluation_list)
    parallel.calc_results(os.path.join(directory_pdb, "results.json"))

    for filename in filenames:
        delete_tb_params(filename, tb_model_name)


# ----------------------------------------------------------------------
