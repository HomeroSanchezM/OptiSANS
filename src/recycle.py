#!/usr/bin/env python3
"""
Recycle workflow: D2O contrast variation scan with a fixed deuteration pattern.

For a given PDB file and amino-acid deuteration pattern, generates SANS
curves for every D2O percentage from 0 to 100 %, evaluates their fitness
against three reference curves, and produces two diagnostic plots.
"""

import csv
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

import numpy as np

# Configure module-level logger (same style as generate_deuterated_pdbs.py)
logger = logging.getLogger(__name__)


# ============================================================================
#                           PUBLIC ENTRY POINT
# ============================================================================

def recycle_workflow(
    pdb_file: str,
    d2o_ref: int,
    amino_acids: List[str],
    output_dir: Optional[str] = "./",
    batch_script: str = "./parallel_process_pdb.sh",
    d2o_step: int = 1,
    q_max: float = 0.3,
    ratio_threshold: float = 0.01,
    n_jobs: int = 150,
    no_default_ref: bool = False,
    gamma: float = 2,
    conc: float = 2.5,
) -> None:
    """
    Full contrast-variation scan pipeline.

    Args:
        pdb_file:        Path to the fully protonated input PDB file.
        d2o_ref:         D2O percentage (0–100) used as the pattern reference.
                         The SANS curve for this exact D2O value (with the given
                         AA deuteration pattern) is copied into ref/ as a reference.
        amino_acids:     List of 3-letter AA codes to deuterate (e.g. ["LEU","LYS"]).
                         Pass an empty list for no AA deuteration.
        output_dir:      Base output directory. Default: "{pdb_stem}_recycle/".
        batch_script:    Path to parallel_process_pdb.sh.
        d2o_step:        Step size for the D2O scan (default 1 = every percent).
        q_max:           q truncation limit for fitness evaluation (A^-1).
        ratio_threshold: Minimum Imax/background ratio to accept a curve.
        n_jobs:          Number of parallel Pepsi-SANS jobs.
        conc:            Pepsi-SANS concentration passed as --conc.

    Raises:
        FileNotFoundError: If pdb_file or batch_script is not found.
        ValueError:        If d2o_ref is not reachable given d2o_step, or if
                           d2o_ref is outside [0, 100].
    """
    # ---- 0. Validate inputs ------------------------------------------------
    _validate_inputs(pdb_file, d2o_ref, amino_acids, batch_script, d2o_step)

    pdb_stem = Path(pdb_file).stem

    # ---- 1. Create directory structure -------------------------------------
    if output_dir is None:
        base_dir = Path(f"{pdb_stem}_recycle")
    else:
        base_dir = Path(output_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    # PDB output folder — naming is CRITICAL for parallel_process_pdb.sh:
    # The script strips the "_deuterated_pdbs" suffix to get the prefix, then
    # creates "{prefix}_primus_out" as the SANS output folder.
    pdb_dir = base_dir / f"{pdb_stem}_recycle_deuterated_pdbs"
    pdb_dir.mkdir(exist_ok=True)
    ref_pdb_dir = pdb_dir / "ref"
    ref_pdb_dir.mkdir(exist_ok=True)

    # Expected SANS output folder (created by the batch script)
    sans_dir = base_dir / f"{pdb_stem}_recycle_primus_out"

    logger.info("=" * 70)
    logger.info("              OPTISANS RECYCLE — Contrast Variation Scan")
    logger.info("=" * 70)
    logger.info(f"Input PDB        : {pdb_file}")
    logger.info(f"Deuterated AAs   : {', '.join(amino_acids) if amino_acids else 'none'}")
    logger.info(f"D2O pattern ref  : {d2o_ref}%")
    logger.info(f"D2O scan range   : 0–100% (step {d2o_step})")
    logger.info(f"PDB directory    : {pdb_dir.absolute()}")
    logger.info(f"SANS directory   : {sans_dir.absolute()}")
    logger.info("=" * 70)

    # ---- 2. Build 20-element deuteration vector ----------------------------
    deuteration_vector = _build_deuteration_vector(amino_acids)

    # ---- 3. Generate main PDB files (_d2o0.pdb … _d2o100.pdb) -------------
    _generate_main_pdbs(pdb_file, deuteration_vector, pdb_dir, d2o_step)

    # ---- 4. Generate reference PDB files -----------------------------------
    if not no_default_ref:
        _generate_reference_pdbs(pdb_file, pdb_stem, ref_pdb_dir)
    else:
        logger.info("Skipping default reference PDB generation (--no-default-ref).")

    # ---- 5. Run Pepsi-SANS via parallel_process_pdb.sh --------------------
    _run_pepsi_sans(pdb_dir, batch_script, n_jobs, conc=conc)

    # ---- 6. Copy pattern reference curve to ref/ --------------------------
    _copy_pattern_reference(pdb_stem, d2o_ref, sans_dir)

    # ---- 7. Fitness evaluation → result.csv --------------------------------
    csv_path = _run_fitness_evaluation(sans_dir, q_max, ratio_threshold,gamma)

    # ---- 8. Plot 1: fitness vs D2O% ----------------------------------------
    _plot_fitness(csv_path, sans_dir)

    # ---- 9. Plot 2: I(0) vs D2O% -------------------------------------------
    _plot_i0(sans_dir)

    logger.info("=" * 70)
    logger.info("RECYCLE PIPELINE COMPLETE")
    logger.info(f"  PDB files      : {pdb_dir.absolute()}")
    logger.info(f"  SANS curves    : {sans_dir.absolute()}")
    logger.info(f"  Fitness CSV    : {csv_path.absolute()}")
    logger.info(f"  Fitness plot   : {(sans_dir / 'fitness_vs_d2o.png').absolute()}")
    logger.info(f"  I(0) plots     : {sans_dir.absolute()}/I0_vs_d2o_*.png")
    logger.info("=" * 70)


# ============================================================================
#                           PRIVATE HELPERS
# ============================================================================

def _validate_inputs(
    pdb_file: str,
    d2o_ref: int,
    amino_acids: List[str],
    batch_script: str,
    d2o_step: int,
) -> None:
    """Raise clear errors for bad inputs before any expensive work starts."""
    from pdb_deuteration import AA_INDEX   # local import to keep module lean

    if not Path(pdb_file).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    if not Path(batch_script).exists():
        raise FileNotFoundError(
            f"Batch script not found: {batch_script}. "
            "Make sure you are running from the project root."
        )
    if not (0 <= d2o_ref <= 100):
        raise ValueError(f"d2o_ref must be in [0, 100], got {d2o_ref}")
    if d2o_step < 1 or d2o_step > 100:
        raise ValueError(f"d2o_step must be in [1, 100], got {d2o_step}")
    if d2o_ref % d2o_step != 0:
        raise ValueError(
            f"d2o_ref={d2o_ref} is not reachable with step={d2o_step}. "
            f"Use a value that is a multiple of {d2o_step}, "
            f"e.g. {(d2o_ref // d2o_step) * d2o_step} or "
            f"{((d2o_ref // d2o_step) + 1) * d2o_step}."
        )
    invalid = [aa for aa in amino_acids if aa.upper() not in AA_INDEX]
    if invalid:
        raise ValueError(
            f"Unknown amino acid code(s): {', '.join(invalid)}. "
            "Use standard 3-letter codes (ALA, ARG, ASN, …)."
        )


def _build_deuteration_vector(amino_acids: List[str]) -> List[bool]:
    """Build a 20-element deuteration vector from a list of 3-letter AA codes."""
    from pdb_deuteration import AMINO_ACIDS, AA_INDEX
    vector = [False] * len(AMINO_ACIDS)
    for aa_code in amino_acids:
        idx = AA_INDEX.get(aa_code.upper())
        if idx is not None:
            vector[idx] = True
    return vector


def _generate_main_pdbs(
    pdb_file: str,
    deuteration_vector: List[bool],
    pdb_dir: Path,
    d2o_step: int,
) -> None:
    """
    Generate one deuterated PDB per D2O value in range(0, 101, d2o_step).

    CRITICAL — Filename format: `_d2o{XX}.pdb`
    `parallel_process_pdb.sh` extracts the D2O integer from the filename with the
    regex `_d2o([0-9]+)` and sets Pepsi-SANS `--d2o XX/100`.
    """
    from pdb_deuteration import PdbDeuteration

    d2o_values = list(range(0, 101, d2o_step))
    logger.info(f"Generating {len(d2o_values)} deuterated PDB files in {pdb_dir.name} …")

    for d2o in d2o_values:
        out_path = pdb_dir / f"_d2o{d2o}.pdb"
        try:
            deut = PdbDeuteration(pdb_file)
            deut.apply_deuteration(deuteration_vector, d2o)
            deut.save(str(out_path))
            logger.debug(f"  Written: {out_path.name}")
        except Exception as exc:
            logger.error(f"  FAILED to generate {out_path.name}: {exc}")
            raise

    logger.info(f"  {len(d2o_values)} PDB files written.")


def _generate_reference_pdbs(
    pdb_file: str,
    pdb_stem: str,
    ref_pdb_dir: Path,
) -> None:
    """
    Generate the 2 reference PDB files with NO amino acid deuteration:
      - {stem}_total_protonation.pdb   (D2O = 0 %,  all-H, no exchange)
      - {stem}_total_deuteration.pdb   (D2O = 100 %, all labile H exchanged)

    CRITICAL — Filename suffixes must end in `_total_protonation` and
    `_total_deuteration` exactly.  `parallel_process_pdb.sh` matches these
    suffixes to set Pepsi-SANS `--d2o 0` and `--d2o 1` respectively.
    """
    from pdb_deuteration import PdbDeuteration, AMINO_ACIDS

    no_deut = [False] * len(AMINO_ACIDS)

    logger.info("Generating 2 reference PDB files in ref/ …")

    # Protonated in H2O (D2O = 0)
    prot = PdbDeuteration(pdb_file)
    prot.apply_deuteration(no_deut, 0)
    prot_path = ref_pdb_dir / f"{pdb_stem}_total_protonation.pdb"
    prot.save(str(prot_path))
    logger.info(f"  Protonation ref : {prot_path.name}")

    # Protonated in D2O (D2O = 100, all labile H exchanged)
    deut = PdbDeuteration(pdb_file)
    deut.apply_deuteration(no_deut, 100)
    deut_path = ref_pdb_dir / f"{pdb_stem}_total_deuteration.pdb"
    deut.save(str(deut_path))
    logger.info(f"  Deuteration ref : {deut_path.name}")


def _run_pepsi_sans(
    pdb_dir: Path,
    batch_script: str,
    n_jobs: int,
    conc: float = 2.5,
) -> None:
    """
    Run Pepsi-SANS in parallel via `parallel_process_pdb.sh`.

    The script is called as:
        parallel_process_pdb.sh <abs_pdb_dir> <n_jobs> --conc <conc>
    and is executed with `cwd` set to the directory that contains the script
    (typically the project root, where `./Pepsi-SANS-Linux/Pepsi-SANS` lives).
    """
    script_path = Path(batch_script).resolve()
    script_cwd = script_path.parent

    cmd = [str(script_path), str(pdb_dir.resolve()), str(n_jobs), "--conc", str(conc)]

    logger.info("=" * 70)
    logger.info("Running Pepsi-SANS simulation …")
    logger.info(f"  Script : {script_path}")
    logger.info(f"  PDB dir: {pdb_dir.resolve()}")
    logger.info(f"  Jobs   : {n_jobs}")
    logger.info(f"  Conc   : {conc}")
    logger.info("=" * 70)

    try:
        result = subprocess.run(
            cmd,
            check=True,
            text=True,
            capture_output=True,
            cwd=str(script_cwd),
        )
        logger.debug(result.stdout)
    except subprocess.CalledProcessError as exc:
        logger.error(f"Pepsi-SANS simulation failed (exit {exc.returncode})")
        logger.error(exc.stderr)
        raise RuntimeError(f"parallel_process_pdb.sh failed: {exc}") from exc
    except FileNotFoundError:
        raise RuntimeError(
            f"Batch script not executable: {script_path}. "
            "Make sure it is executable (`chmod +x parallel_process_pdb.sh`) "
            "and that GNU parallel is installed."
        )

    logger.info("Pepsi-SANS simulation complete.")


def _copy_pattern_reference(
    pdb_stem: str,
    d2o_ref: int,
    sans_dir: Path,
) -> None:
    """
    Copy the SANS curve for the given deuteration pattern at d2o_ref% from the
    main output folder to ref/ as the 3rd reference curve.

    Source : {sans_dir}/_d2o{d2o_ref}.dat
    Dest   : {sans_dir}/ref/{pdb_stem}_pattern_d2o{d2o_ref}.dat
    """
    src = sans_dir / f"_d2o{d2o_ref}.dat"
    if not src.exists():
        raise FileNotFoundError(
            f"Expected SANS curve not found: {src}. "
            f"Pepsi-SANS may have failed for D2O={d2o_ref}%."
        )

    ref_dir = sans_dir / "ref"
    ref_dir.mkdir(exist_ok=True)   # should already exist from the script, but be safe

    dst = ref_dir / f"{pdb_stem}_pattern_d2o{d2o_ref}.dat"
    shutil.copy2(src, dst)
    logger.info(f"Copied pattern reference: {dst.name}  ->  ref/")


def _run_fitness_evaluation(
    sans_dir: Path,
    q_max: float,
    ratio_threshold: float,
    gamma: float,
) -> Path:
    """
    Evaluate fitness for all .dat files in sans_dir against the 3 reference
    curves in sans_dir/ref/.  Write results to sans_dir/result.csv.

    Returns the Path to the written CSV file.
    """
    from fitness_evaluation import evaluate_population_fitness

    logger.info("=" * 70)
    logger.info("Running fitness evaluation …")

    fitness_scores, dat_files, ratios = evaluate_population_fitness(
        directory=str(sans_dir),
        q_max=q_max,
        ratio_threshold=ratio_threshold,
        gamma=gamma,
    )

    # Summary
    n_pass = int(np.sum(np.array(fitness_scores) > 0))
    logger.info(
        f"  Evaluated {len(fitness_scores)} curves — "
        f"{n_pass} passed ratio check — "
        f"best fitness: {float(np.max(fitness_scores)):.4e}"
    )

    # Write result.csv  (format expected by plot_contrast_fitness.py)
    csv_path = sans_dir / "result.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["filename", "fitness_score", "ratio"])
        for file_path, score, ratio in zip(dat_files, fitness_scores, ratios):
            writer.writerow([
                Path(file_path).name,
                f"{score:.4e}",
                f"{ratio:.4e}",
            ])

    logger.info(f"  Fitness CSV     : {csv_path}")
    return csv_path


def _plot_fitness(csv_path: Path, sans_dir: Path) -> None:
    """Generate the fitness-vs-D2O% plot and save it next to result.csv."""
    try:
        from plot_contrast_fitness import plot_fitness_from_csv
    except ImportError as exc:
        logger.warning(f"Could not import plot_contrast_fitness: {exc}. Skipping fitness plot.")
        return

    output = sans_dir / "fitness_vs_d2o.png"
    try:
        plot_fitness_from_csv(
            csv_file=str(csv_path),
            output=str(output),
            no_labels=False,
            label_step=5,
        )
        logger.info(f"  Fitness plot    : {output}")
    except Exception as exc:
        logger.warning(f"Fitness plot failed: {exc}")


def _plot_i0(sans_dir: Path) -> None:
    """Generate I(0) vs D2O% plots and save them in sans_dir."""
    try:
        from plot_i0_d2o import plot_i0_vs_d2o
    except ImportError as exc:
        logger.warning(f"Could not import plot_i0_d2o: {exc}. Skipping I(0) plots.")
        return

    output_prefix = str(sans_dir / "I0_vs_d2o")
    try:
        plot_i0_vs_d2o(str(sans_dir), output_prefix)
        logger.info(f"  I(0) plots      : {output_prefix}_*.png")
    except Exception as exc:
        logger.warning(f"I(0) plot failed: {exc}")
