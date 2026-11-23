#!/usr/bin/env python3
"""
Fitness Evaluation for SANS Deuteration Profiles
================================================

Evaluates fitness of simulated SANS curves by comparing them to ALL reference
curves found in the ref/ subdirectory.

Fitness = product(area_i for each reference) * ratio * 10 000

"""

import argparse
import csv
import glob
import logging
import os
import re
from typing import Dict, List, Tuple

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import griddata
from scipy.optimize import curve_fit

# ============================================================================
#                           LOGGING CONFIGURATION
# ============================================================================

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
#                           ARGUMENT PARSING
# ============================================================================


def parse_arguments():
    """
    Parses command-line arguments for fitness evaluation.

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Evaluate fitness of simulated SANS curves against reference data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Directory with simulated .dat files and a ref/ subfolder",
    )
    ref_group = parser.add_argument_group(
        "Reference files (deprecated - all .dat in ref/ are used)"
    )
    ref_group.add_argument(
        "--deut-ref",
        type=str,
        default=None,
        help="Deprecated: ignored. All .dat files in ref/ are used.",
    )
    ref_group.add_argument(
        "--prot-ref",
        type=str,
        default=None,
        help="Deprecated: ignored. All .dat files in ref/ are used.",
    )
    param_group = parser.add_argument_group("Evaluation parameters")
    param_group.add_argument(
        "--q-max",
        type=float,
        default=0.3,
        help="Maximum q value for truncation (default: 0.3 Å⁻¹)",
    )
    param_group.add_argument(
        "--ratio-threshold",
        type=float,
        default=0.01,
        help="Minimum Imax/background ratio to accept a curve (default: 0.01)",
    )
    param_group.add_argument(
        "--gamma",
        type=float,
        default=2,
        help="Exponent for Imax/background ratio in fitness "
        "formula: fitness = product(areas) * ratio^gamma. "
        "Default: 2 (quadratic). Use 0 to ignore ratio, "
        "1 for linear weighting.",
    )
    param_group.add_argument(
        "--i0-threshold", type=float, default=None, help=argparse.SUPPRESS
    )
    output_group = parser.add_argument_group("Output control")
    output_group.add_argument("-v", "--verbose", action="store_true")
    output_group.add_argument("-q", "--quiet", action="store_true")
    output_group.add_argument(
        "--csv",
        type=str,
        default=None,
        metavar="FILE",
        help="Write filename, fitness_score and ratio to a CSV file",
    )
    args = parser.parse_args()
    if args.quiet:
        logger.setLevel(logging.WARNING)
    elif args.verbose:
        logger.setLevel(logging.DEBUG)
    if args.deut_ref or args.prot_ref:
        logger.warning(
            "--deut-ref / --prot-ref are deprecated: all .dat files in ref/ are used."
        )
    return args


# ============================================================================
#                           FILE DISCOVERY
# ============================================================================


def find_dat_files(directory):
    """
    Find all .dat files in the specified directory.

    Args:
        directory (str): Path to the directory

    Returns:
        List[str]: Sorted list of paths to .dat files

    Raises:
        NotADirectoryError: If the path is not a directory
        FileNotFoundError: If no .dat files are found
    """
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"'{directory}' is not a valid directory")
    dat_files = [
        f for f in glob.glob(os.path.join(directory, "*.dat")) if os.path.isfile(f)
    ]
    if not dat_files:
        raise FileNotFoundError(f"No .dat files found in '{directory}'")
    return sorted(dat_files)


def find_all_reference_files(ref_dir):
    """Return all .dat files found in ref_dir."""
    if not os.path.isdir(ref_dir):
        raise NotADirectoryError(f"Reference directory not found: {ref_dir}")
    ref_files = [
        f for f in glob.glob(os.path.join(ref_dir, "*.dat")) if os.path.isfile(f)
    ]
    if not ref_files:
        raise FileNotFoundError(
            f"No .dat reference files found in {ref_dir}. "
            "Ensure the ref/ subdirectory contains at least one .dat file."
        )
    return sorted(ref_files)


# Kept for backward compatibility with external callers
def find_reference_files(ref_dir, deut_ref=None, prot_ref=None):
    """Deprecated: use find_all_reference_files instead."""
    logger.warning("find_reference_files is deprecated; use find_all_reference_files.")
    paths = find_all_reference_files(ref_dir)
    return paths[0], paths[1] if len(paths) > 1 else paths[0]


# ============================================================================
#                           PARSING & UTILITIES
# ============================================================================


def parse_sans_file(filepath):
    q_values, I_values = [], []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    q_values.append(float(parts[0]))
                    I_values.append(float(parts[1]))
                except ValueError:
                    continue
    q, I = np.array(q_values), np.array(I_values)
    if len(q) == 0:
        raise ValueError(f"No valid data points in file: {filepath}")

    return q, I


def truncate_to_q_max(q, I, q_max):
    """
    Truncate arrays to keep only values where q <= q_max.

    Args:
        q (np.array): q values
        I (np.array): I values
        q_max (float): Maximum q value

    Returns:
        tuple: (q_truncated, I_truncated)
    """
    mask = q <= q_max
    return q[mask], I[mask]


def regrid_curve(q_source, I_source, q_target):
    """
    Interpolate a curve onto a target q grid.

    Args:
        q_source (np.array): Original q values
        I_source (np.array): Original I values
        q_target (np.array): Target q values

    Returns:
        np.array: Interpolated I values
    """
    # Use cubic interpolation to regrid
    I_regridded = griddata(q_source, I_source, q_target, method="cubic")

    # Fallback to linear interpolation for NaN values
    if np.any(np.isnan(I_regridded)):
        # Fall back to linear interpolation for NaN values
        I_linear = griddata(q_source, I_source, q_target, method="linear")
        I_regridded = np.where(np.isnan(I_regridded), I_linear, I_regridded)

    return I_regridded


def scale_curves(q, I_reference, I_to_scale):
    """
    Scale I_to_scale to match I_reference using linear model: I_scaled = a * I_to_scale + b.

    Args:
        q (np.array): q values (same for both curves)
        I_reference (np.array): Reference intensity
        I_to_scale (np.array): Intensity to scale

    Returns:
        tuple: (I_scaled, a, b, chi2)
    """

    def linear_model(I, a, b):
        return a * I + b

    try:
        params, _ = curve_fit(linear_model, I_to_scale, I_reference)
        a, b = params
        I_scaled = a * I_to_scale + b

        # Calculate chi^2 (area difference metric)
        chi2 = np.sum((I_reference - I_scaled) ** 2)

        return I_scaled, a, b, chi2

    except Exception as e:
        logger.warning(f"Scaling failed: {e}")
        return I_to_scale, 1.0, 0.0, np.inf


def calculate_area_difference(q, I1, I2):
    """
    Compute the area between two curves using Simpson integration.
    """
    return simpson(np.abs(I1 - I2), x=q)


# ============================================================================
#                     RATIO CHECK
# ============================================================================


def extract_d2o_from_filename(filename):
    """
    Extract the D2O percentage from a curve filename.

    Expected filename pattern: genXX_chrXXX_d2oYY_deutAAZZ.dat (or .out)
    Example: gen01_chr003_d2o42_deutAA05.dat -> returns 42

    Args:
        filename (str): Basename or full path of the curve file

    Returns:
        int or None: D2O percentage, or None if not found in filename
    """
    basename = os.path.basename(filename)
    match = re.search(r"_d2o(\d+)", basename)
    return int(match.group(1)) if match else None


def compute_incoherent_background(d2o_percent):
    """
    Compute the incoherent scattering background for a given D2O percentage.

    Formula: background = -0.0117 * d2o_percent + 1.25

    Args:
        d2o_percent (int or float): D2O percentage (0-100)

    Returns:
        float: Background value
    """
    return -0.0117 * d2o_percent + 1.25


def ratio_check(I, d2o_percent, ratio_threshold=0.01):
    """
    Check if a curve is valid based on the Imax / background ratio.

    A curve is accepted if:
        ratio = Imax / background > threshold
    where:
        Imax       = maximum value of I (not necessarily I[0])
        background = -0.0117 * d2o_percent + 1.25

    Args:
        I (np.array): Simulated intensity array
        d2o_percent (int or float): D2O percentage extracted from the filename
        threshold (float): Minimum ratio to accept the curve (default: 0.01)

    Returns:
        bool: True if the curve passes the ratio check, False otherwise
    """
    background = compute_incoherent_background(d2o_percent)
    if background <= 0:
        logger.warning(
            f"Non-positive background ({background:.4f}) for D2O={d2o_percent}%; curve rejected"
        )
        return 0.0
    if np.max(I) < ratio_threshold:  # discard courves with Imax<0.01
        logger.debug(f"Imax too low (Imax={np.max(I):.4e}); curve rejected")
        return 0.0
    ratio = np.max(I) / background
    logger.debug(
        f"Ratio check: Imax={np.max(I):.4e}, background={background:.4f}, ratio={ratio:.4f}"
    )
    return ratio


# ============================================================================
#  SCALING AND FITNESS
# ============================================================================


def scaling_and_compare(q, I, references, q_max, ratio, gamma=2):
    """
    Compute fitness as the product of area differences against all references,
    multiplied by ratio^gamma.

    Args:
        q:          q values of the simulated curve
        I:          Intensity of the simulated curve
        references: List of (q_ref, I_ref) tuples for each reference curve
        q_max:      Truncation limit
        ratio:      Imax / background ratio for the simulated curve
        gamma:      Exponent applied to ratio (default: 2).
                    0 = no ratio weighting, 1 = linear, 2 = quadratic.

    Returns:
        float: product(area_i) * ratio^gamma
    """
    q_trunc, I_trunc = truncate_to_q_max(q, I, q_max)

    fitness_product = 1.0
    for q_ref, I_ref in references:
        q_ref_trunc, I_ref_trunc = truncate_to_q_max(q_ref, I_ref, q_max)
        I_ref_regrid = regrid_curve(q_ref_trunc, I_ref_trunc, q_trunc)
        I_scaled, _, _, _ = scale_curves(q_trunc, I_ref_regrid, I_trunc)
        area = calculate_area_difference(q_trunc, I_scaled, I_ref_regrid)
        fitness_product *= area

    return fitness_product * (ratio**gamma)


def fitness(q, I, references, q_max, file_path, ratio_threshold=0.01, gamma=2):
    """
    Compute raw fitness score for a single curve against all references.

    Args:
        q, I:            Simulated curve data
        references:      List of (q_ref, I_ref) tuples
        q_max:           q truncation limit
        file_path:       Used to extract D2O percentage for ratio check
        ratio_threshold: Minimum ratio to accept a curve
        gamma:           Exponent applied to ratio in fitness formula

    Returns:
        (float, float): (fitness_score, ratio)
    """
    d2o_percent = extract_d2o_from_filename(file_path)

    if d2o_percent is None:
        logger.debug(
            f"Could not extract D2O from '{os.path.basename(file_path)}'; ratio check skipped"
        )
        ratio = 0.0
    else:
        ratio = ratio_check(I, d2o_percent, ratio_threshold)
        if ratio < ratio_threshold:
            logger.debug(
                f"Ratio check failed for '{os.path.basename(file_path)}': ratio={ratio:.4f}"
            )
            return 0.0, ratio

    return scaling_and_compare(q, I, references, q_max, ratio, gamma=gamma), ratio


def normalize_fitness(fitness_values):
    """
    Normalize fitness scores to [0, 1] (higher area → higher fitness).

    Args:
        fitness_values (list or np.array): Raw fitness scores

    Returns:
        np.array: Normalized scores
    """
    arr = np.array(fitness_values, dtype=float)

    if len(arr) == 0:
        return arr

    # Identify zero scores (failed ratio check)
    zero_mask = arr == 0.0
    non_zero = arr[~zero_mask]

    if len(non_zero) == 0:
        return arr
    min_val, max_val = np.min(non_zero), np.max(non_zero)
    if max_val == min_val:
        # All non-zero equal → assign 0.5
        arr[~zero_mask] = 0.5
    else:
        arr[~zero_mask] = (non_zero - min_val) / (max_val - min_val)

    return arr


# ============================================================================
#                           MAIN EVALUATION
# ============================================================================


def evaluate_population_fitness(
    directory,
    deut_ref=None,
    prot_ref=None,
    q_max=0.3,
    i0_threshold=None,
    ratio_threshold=0.01,
    gamma=2,
):
    """
    Evaluate fitness for all .dat files in directory against all references in ref/.

    Args:
        directory:       Directory with simulated .dat files and a ref/ subfolder
        deut_ref:        Deprecated, ignored
        prot_ref:        Deprecated, ignored
        q_max:           q truncation limit
        i0_threshold:    Deprecated, ignored
        ratio_threshold: Minimum Imax/background ratio
        gamma:           Exponent applied to ratio in fitness formula (default: 2)

    Returns:
        (raw_fitness_scores, file_paths, ratios)
    """
    if deut_ref or prot_ref:
        logger.debug(
            "deut_ref/prot_ref are deprecated; all .dat files in ref/ are used."
        )

    sim_files = find_dat_files(directory)
    logger.info(f"Found {len(sim_files)} simulated .dat files")

    ref_dir = os.path.join(directory, "ref")
    ref_paths = find_all_reference_files(ref_dir)
    logger.info(f"Found {len(ref_paths)} reference file(s) in ref/")

    references = []
    for ref_path in ref_paths:
        q_ref, I_ref = parse_sans_file(ref_path)
        references.append((q_ref, I_ref))
        logger.info(
            f"  Loaded reference: {os.path.basename(ref_path)} [{len(q_ref)} points]"
        )

    raw_scores, valid_files, ratios = [], [], []

    for file_path in sim_files:
        try:
            q, I = parse_sans_file(file_path)
            score, ratio = fitness(
                q, I, references, q_max, file_path, ratio_threshold, gamma=gamma
            )
            raw_scores.append(score)
            valid_files.append(file_path)
            ratios.append(ratio)
            logger.debug(
                f"{os.path.basename(file_path)}: ratio={ratio:.4f} score={score:.6e}"
            )
        except Exception as e:
            logger.warning(f"Skipping {os.path.basename(file_path)}: {e}")
            raw_scores.append(0.0)
            valid_files.append(file_path)
            ratios.append(0.0)

    return np.array(raw_scores, dtype=float), valid_files, np.array(ratios, dtype=float)


# ============================================================================
#  ENSEMBLE (MULTI-PDB) HELPERS
# ============================================================================


def mean_curve(file_paths: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse all curves, regrid onto the first file's q-grid, average.

    Args:
        file_paths: List of .dat file paths (one per conformation).

    Returns:
        (q, I_mean, I_std): q-grid of the first file, mean and std of the
        regridded intensities along axis 0.
    """
    curves = [parse_sans_file(f) for f in file_paths]
    q0 = curves[0][0]
    I_stack = np.stack([regrid_curve(q, I, q0) for q, I in curves])
    return q0, I_stack.mean(axis=0), I_stack.std(axis=0)


def group_dat_files_by_chromosome(directory) -> Dict[str, List[str]]:
    """Group the .dat files of *directory* by chromosome-level stem.

    The ``_confNNNN`` ensemble suffix is stripped by
    :func:`src.ensemble.chrom_key_from_filename`, so all per-conformation
    curves of one chromosome land in the same group.
    """
    from src.ensemble import chrom_key_from_filename

    groups: Dict[str, List[str]] = {}
    for f in find_dat_files(directory):
        groups.setdefault(chrom_key_from_filename(f), []).append(f)
    return groups


def evaluate_population_fitness_mean(
    directory,
    q_max=0.3,
    ratio_threshold=0.01,
    gamma=2,
):
    """Ensemble variant of :func:`evaluate_population_fitness`.

    The fitness of each chromosome is computed on the mean of its
    per-conformation .dat files.  The D2O ratio check still works:
    ``_d2o(\\d+)`` matches inside conf-suffixed filenames, so the first
    file of each group is passed to :func:`fitness` unchanged.

    Args:
        directory:       Directory with simulated .dat files and a ref/ subfolder
        q_max:           q truncation limit
        ratio_threshold: Minimum Imax/background ratio
        gamma:           Exponent applied to ratio in fitness formula

    Returns:
        (raw_scores, chrom_keys, ratios, groups)
    """
    groups = group_dat_files_by_chromosome(directory)
    logger.info(
        f"Found {sum(len(v) for v in groups.values())} simulated .dat files "
        f"in {len(groups)} chromosome group(s)"
    )

    ref_dir = os.path.join(directory, "ref")
    ref_paths = find_all_reference_files(ref_dir)
    logger.info(f"Found {len(ref_paths)} reference file(s) in ref/")
    references = [parse_sans_file(p) for p in ref_paths]

    scores, keys, ratios = [], [], []
    for key in sorted(groups):
        files = groups[key]
        q, I_mean, _ = mean_curve(files)
        try:
            score, ratio = fitness(
                q, I_mean, references, q_max, files[0],
                ratio_threshold, gamma=gamma,
            )
        except Exception as e:
            logger.warning(f"Skipping chromosome group {key}: {e}")
            score, ratio = 0.0, 0.0
        scores.append(score)
        keys.append(key)
        ratios.append(ratio)
        logger.debug(f"{key}: ratio={ratio:.4f} score={score:.6e} ({len(files)} curve(s))")

    return (np.array(scores, dtype=float), keys,
            np.array(ratios, dtype=float), groups)


def average_ref_curves(primus_dir,
                       kinds=("total_deuteration", "total_protonation")) -> None:
    """Replace per-conformation reference curves in ``primus_dir/ref/`` by one
    mean curve per kind (``mean_<kind>.dat``).

    Called once in ensemble mode, after the generation-0 full-mode batch run
    and before the first fitness evaluation, so that ``ref/`` contains exactly
    one mean curve per reference kind and the "all .dat in ref/" contract of
    :func:`evaluate_population_fitness` keeps holding.

    The 3-column mean file reloads fine through :func:`parse_sans_file`
    (only columns 0 and 1 are read).
    """
    ref_dir = os.path.join(primus_dir, "ref")
    for kind in kinds:
        files = sorted(glob.glob(
            os.path.join(ref_dir, f"*_conf[0-9][0-9][0-9][0-9]_{kind}.dat")))
        if not files:
            continue
        q, I_mean, I_std = mean_curve(files)
        out = os.path.join(ref_dir, f"mean_{kind}.dat")
        with open(out, "w") as fh:
            fh.write(f"# mean over {len(files)} conformation(s)\n")
            fh.write("# q  I_mean  I_std\n")
            for qi, mi, si in zip(q, I_mean, I_std):
                fh.write(f"{qi:.6f}\t{mi:.6e}\t{si:.6e}\n")
        for f in files:
            os.unlink(f)
        logger.info(
            f"Mean reference: {os.path.basename(out)} "
            f"({len(files)} curve(s) averaged)"
        )


def main():
    args = parse_arguments()
    try:
        fitness_scores, sim_files, ratios = evaluate_population_fitness(
            directory=args.directory,
            q_max=args.q_max,
            ratio_threshold=args.ratio_threshold,
            gamma=args.gamma,
        )
    except (NotADirectoryError, FileNotFoundError, ValueError) as e:
        logger.error(f"Evaluation failed: {e}")
        import sys

        sys.exit(1)

    if not args.quiet:
        print("\n" + "=" * 60)
        print("FITNESS EVALUATION SUMMARY")
        print("=" * 60)
        print(f"Directory: {args.directory}")
        print(f"Q max: {args.q_max:.3f} Å⁻¹")
        print(f"Ratio threshold: {args.ratio_threshold}")
        print(f"Files evaluated: {len(fitness_scores)}")
        print(
            f"Files passing ratio check: {np.sum(fitness_scores > 0)}/{len(fitness_scores)}"
        )
        print(f"Best fitness: {np.max(fitness_scores):.4e}")
        if len(fitness_scores) > 0 and np.max(fitness_scores) > 0:
            best_idx = np.argmax(fitness_scores)
            print(f"Best file: {os.path.basename(sim_files[best_idx])}")
        print(f"Average fitness: {np.mean(fitness_scores):.4e}")
        print("=" * 60 + "\n")

    for score in fitness_scores:
        print(f"{score:.4e}")

    if args.csv:
        csv_path = args.csv
        try:
            with open(csv_path, "w", newline="") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["filename", "fitness_score", "ratio"])
                for file_path, score, ratio in zip(sim_files, fitness_scores, ratios):
                    writer.writerow(
                        [os.path.basename(file_path), f"{score:.4e}", f"{ratio:.4e}"]
                    )
            if not args.quiet:
                logger.info(
                    f"Results saved to '{csv_path}' ({len(fitness_scores)} rows)"
                )
        except OSError as e:
            logger.error(f"Could not write CSV file '{csv_path}': {e}")


if __name__ == "__main__":
    import sys

    main()
