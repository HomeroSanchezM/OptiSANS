#!/usr/bin/env python3
"""OptiSANS — Unified CLI for protein deuteration optimization for SANS."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import typer

# Ensure sibling modules in src/ are importable when running as an entry point.
_src_dir = str(Path(__file__).resolve().parent)
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)


app = typer.Typer(
    name="optisans",
    help="Protein deuteration optimization for SANS experiments.",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)

BANNER = r"""

 ________  ________  _________  ___  ________  ________  ________   ________
|\   __  \|\   __  \|\___   ___\\  \|\   ____\|\   __  \|\   ___  \|\   ____\
\ \  \|\  \ \  \|\  \|___ \  \_\ \  \ \  \___|\ \  \|\  \ \  \\ \  \ \  \___|_
 \ \  \\\  \ \   ____\   \ \  \ \ \  \ \_____  \ \   __  \ \  \\ \  \ \_____  \
  \ \  \\\  \ \  \___|    \ \  \ \ \  \|____|\  \ \  \ \  \ \  \\ \  \|____|\  \
   \ \_______\ \__\        \ \__\ \ \__\____\_\  \ \__\ \__\ \__\\ \__\____\_\  \
    \|_______|\|__|         \|__|  \|__|\_________\|__|\|__|\|__| \|__|\_________\
                                       \|_________|                   \|_________|

"""


@app.callback(invoke_without_command=True)
def _banner(ctx: typer.Context) -> None:
    """Print the OptiSANS banner, then run the requested subcommand."""
    typer.echo(BANNER)
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


VALID_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


@app.command()
def aa():
    """Affiche la liste des 20 acides aminés standards avec leurs codes."""
    AA_TABLE = [
        ("ALA", "A", "Alanine"),
        ("ARG", "R", "Arginine"),
        ("ASN", "N", "Asparagine"),
        ("ASP", "D", "Aspartic acid"),
        ("CYS", "C", "Cysteine"),
        ("GLU", "E", "Glutamic acid"),
        ("GLN", "Q", "Glutamine"),
        ("GLY", "G", "Glycine"),
        ("HIS", "H", "Histidine"),
        ("ILE", "I", "Isoleucine"),
        ("LEU", "L", "Leucine"),
        ("LYS", "K", "Lysine"),
        ("MET", "M", "Methionine"),
        ("PHE", "F", "Phenylalanine"),
        ("PRO", "P", "Proline"),
        ("SER", "S", "Serine"),
        ("THR", "T", "Threonine"),
        ("TRP", "W", "Tryptophan"),
        ("TYR", "Y", "Tyrosine"),
        ("VAL", "V", "Valine"),
    ]
    typer.echo("\n Amino acids available for deuteration\n")
    typer.echo(f"  {'Code 3':<10} {'Code 1':<10} {'Name'}")
    typer.echo("  " + "─" * 38)
    for code3, code1, name in AA_TABLE:
        typer.echo(f"  {code3:<10} {code1:<10} {name}")
    typer.echo("")
    typer.echo(
        "  Usage : --aa 'LEU LYS PRO'  or  --aa LEU,LYS,PRO\n"
        "  Note  : ASN+ASP and GLU+GLN are always deuterated together (linked pairs).\n"
    )


@app.command()
def run(
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Source PDB file (all hydrogens must be explicit and protonated).",
    ),
    population_size: Optional[int] = typer.Option(
        None,
        "-p",
        "--population-size",
        help="Population size (must be a multiple of 3).",
    ),
    generations: Optional[int] = typer.Option(
        None,
        "-g",
        "--generations",
        help="Maximum number of generations to run.",
    ),
    elitism: Optional[int] = typer.Option(
        None,
        "-e",
        "--elitism",
        help="Number of elite individuals preserved (must be <= population_size / 3).",
    ),
    d2o_var: Optional[int] = typer.Option(
        None,
        "--d2o-var",
        help="Maximum D2O variation per mutation (0-100).",
    ),
    seed: Optional[int] = typer.Option(
        None,
        "--seed",
        help="Random seed for reproducibility.",
    ),
    patience: Optional[int] = typer.Option(
        None,
        "--patience",
        help=(
            "Early stopping: number of consecutive generations without fitness "
            "improvement before stopping. Default: 50. Set to 0 to disable."
        ),
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        exists=True,
        help="config.ini file (CLI arguments override INI values).",
    ),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        help=(
            "Output directory for generated PDB files. "
            "Default: '<pdb_basename>_deuterated_pdbs/' in the current folder."
        ),
    ),
    batch_script: Optional[Path] = typer.Option(
        None,
        "--batch-script",
        help="Path to parallel_process_pdb.sh.",
    ),
    conc: float = typer.Option(
        2.5,
        "--conc",
        help="Pepsi-SANS concentration passed as --conc (default: 2.5).",
    ),
    q_max: Optional[float] = typer.Option(
        None,
        "--q-max",
        help="Maximum q value for fitness evaluation (A^-1).",
    ),
    ratio_threshold: Optional[float] = typer.Option(
        None,
        "--ratio-threshold",
        help="Minimum Imax/background ratio to accept a curve.",
    ),
    gamma: Optional[float] = typer.Option(
        None,
        "--gamma",
        help=(
            "Exponent for Imax/background ratio in fitness formula: "
            "fitness = product(areas) * ratio^gamma. "
            "Default: 2 (quadratic). 0 = ignore ratio, 1 = linear."
        ),
    ),
    d2o_values: Optional[List[int]] = typer.Option(
        None,
        "--d2o",
        help="Lock D2O to fixed values (repeat flag, e.g. --d2o 0 --d2o 42 --d2o 100).",
    ),
    no_default_ref: bool = typer.Option(
        False,
        "--no-default-ref",
        help="Do not create the default protonated-in-D2O / H2O reference PDBs.",
    ),
    ref: Optional[List[Path]] = typer.Option(
        None,
        "--ref",
        help=(
            "Additional reference PDB file(s) to copy into ref/ and use for "
            "fitness evaluation. Repeat the flag for multiple files: "
            "--ref ref1.pdb --ref ref2.pdb. "
            "Can be combined with or without --no-default-ref."
        ),
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Run the full genetic algorithm on a PDB file.

    Outputs are written to the following directories by default:
      '<pdb_basename>_deuterated_pdbs/'   — generated deuterated PDB files
      '<pdb_basename>_primus_out/'        — Pepsi-SANS simulation results
      '<pdb_basename>_final_results/'     — best solution summary
    """
    argv = ["generate_deuterated_pdbs", str(pdb_file)]
    if population_size is not None:
        argv += ["-p", str(population_size)]
    if generations is not None:
        argv += ["-g", str(generations)]
    if elitism is not None:
        argv += ["-e", str(elitism)]
    if d2o_var is not None:
        argv += ["--d2o-var", str(d2o_var)]
    if seed is not None:
        argv += ["--seed", str(seed)]
    if patience is not None:
        argv += ["--patience", str(patience)]
    if config is not None:
        argv += ["--config", str(config)]
    if output_dir is not None:
        argv += ["--output_dir", str(output_dir)]
    if batch_script is not None:
        argv += ["--batch_script", str(batch_script)]
    argv += ["--conc", str(conc)]
    if q_max is not None:
        argv += ["--q-max", str(q_max)]
    if ratio_threshold is not None:
        argv += ["--ratio-threshold", str(ratio_threshold)]
    if gamma is not None:
        argv += ["--gamma", str(gamma)]
    if d2o_values:
        argv += ["--d2o"] + [str(v) for v in d2o_values]
    if no_default_ref:
        argv += ["--no_default_ref"]
    if ref:
        argv += ["--ref"] + [str(r) for r in ref]
    if verbose:
        argv += ["--verbose"]

    sys.argv = argv
    try:
        from generate_deuterated_pdbs import main as ga_main
    except ImportError as exc:
        typer.echo(f"Import error: {exc}", err=True)
        typer.echo(
            "Make sure you are running inside the pixi environment: "
            "pixi run optisans run ...",
            err=True,
        )
        raise typer.Exit(1)
    ga_main()


@app.command()
def deuterate(
    source: Path = typer.Argument(
        ...,
        exists=True,
        help="Input PDB file or directory containing PDB files.",
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output PDB file or output directory when processing a folder.",
    ),
    d2o: float = typer.Option(
        0.0,
        "--d2o",
        help="D2O percentage for labile hydrogen exchange (0-100).",
    ),
    amino_acids: Optional[str] = typer.Option(
        None,
        "-a",
        "--aa",
        help=(
            "Amino acid types to deuterate. Separate multiple codes with spaces or "
            "commas (e.g. --aa 'LEU LYS PRO' or --aa LEU,LYS,PRO). "
            "Run 'optisans aa' to see all available codes."
        ),
    ),
    all_aa: bool = typer.Option(
        False,
        "--all",
        help="Deuterate all amino acid types.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging.",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        help=(
            "INI configuration file for deuteration parameters "
            "(sections [DEUTERATION] and [AMINO_ACIDS]). "
            "CLI arguments always override values from the config file. "
            "See pdb_config.ini in the project root for an example."
        ),
    ),
):
    """Deuterate one or more PDB files according to a given specification."""
    import re as _re

    aa_list: List[str] = []
    if amino_acids:
        aa_list = [
            x.strip().upper() for x in _re.split(r"[,\s]+", amino_acids) if x.strip()
        ]
    if all_aa:
        aa_list = [aa.code_3 for aa in VALID_AA]

    for aa in aa_list:
        if aa not in VALID_AA:
            typer.echo(f"Error: invalid amino acid code: {aa}", err=True)
            typer.echo(
                f"Valid codes: {', '.join(sorted(VALID_AA))}",
                err=True,
            )
            raise typer.Exit(1)

    if config is not None:
        try:
            from pdb_deuteration import load_config_ini, merge_config, validate_config
            ini_cfg = load_config_ini(str(config))
        except (FileNotFoundError, ValueError) as exc:
            typer.echo(f"Error loading config file: {exc}", err=True)
            raise typer.Exit(1)
    else:
        ini_cfg = {}

    if source.is_dir():
        pdb_files = sorted(p.name for p in source.iterdir() if p.is_file() and p.suffix.lower() == ".pdb")
        if not pdb_files:
            typer.echo(f"No .pdb files found in directory: {source}", err=True)
            raise typer.Exit(1)
        if not output.exists():
            output.mkdir(parents=True, exist_ok=True)
        if not output.is_dir():
            typer.echo(
                "When SOURCE is a directory, --output must be a directory.",
                err=True,
            )
            raise typer.Exit(1)

        typer.echo(f"Processing {len(pdb_files)} PDB files from: {source}")
        typer.echo(f"Output directory: {output}")

        ok = 0
        ko = 0
        for idx, name in enumerate(pdb_files, start=1):
            in_path = source / name
            out_path = output / name

            try:
                from pdb_deuteration import build_batch_config, deuterate_file
            except ImportError as exc:
                typer.echo(f"Import error: {exc}", err=True)
                raise typer.Exit(1)

            try:
                cfg = build_batch_config(
                    input_pdb=str(in_path),
                    output_pdb=str(out_path),
                    d2o_percent=d2o,
                    all_aa=all_aa,
                    amino_acids=aa_list or None,
                    config_path=str(config) if config else None,
                    verbose=verbose,
                )
            except ValueError as exc:
                typer.echo(f"[{idx}/{len(pdb_files)}] Invalid config for {name}: {exc}", err=True)
                ko += 1
                continue

            typer.echo(f"[{idx}/{len(pdb_files)}] Deuterating {name} ...")
            success = deuterate_file(
                input_pdb=cfg["input_pdb"],
                output_pdb=cfg["output_pdb"],
                d2o_percent=cfg["d2o_percent"],
                deuteration_vector=cfg["deuteration_vector"],
                verbose=cfg["verbose"],
            )
            if success:
                ok += 1
            else:
                ko += 1

        typer.echo(f"Done. succeeded={ok}, failed={ko} out of {len(pdb_files)} files.")
        if ko:
            raise typer.Exit(1)
        return

    single_input = source
    if output.exists() and output.is_dir():
        single_output = output / source.name
    else:
        single_output = output

    argv = ["pdb_deuteration"]
    argv += [
        "-i",
        str(single_input),
        "-o",
        str(single_output),
        "--d2o",
        str(d2o),
    ]
    if config is not None:
        argv += [str(config)]
    if all_aa:
        argv.append("--all")
    elif aa_list:
        for aa in aa_list:
            argv.append(f"--{aa}")
    if verbose:
        argv.append("--verbose")

    sys.argv = argv
    try:
        from pdb_deuteration import main as deut_main
    except ImportError as exc:
        typer.echo(f"Import error: {exc}", err=True)
        typer.echo(
            "Make sure you are running inside the pixi environment: "
            "pixi run optisans deuterate ...",
            err=True,
        )
        raise typer.Exit(1)
    deut_main()


@app.command()
def evaluate(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=True,
        help="Directory containing .dat files and a ref/ subfolder.",
    ),
    q_max: float = typer.Option(
        0.3,
        "--q-max",
        help="Maximum q value for truncation (A^-1).",
    ),
    ratio_threshold: float = typer.Option(
        0.01,
        "--ratio-threshold",
        help="Minimum Imax/background ratio to accept a curve.",
    ),
    gamma: float = typer.Option(
        2,
        "--gamma",
        help=(
            "Exponent for Imax/background ratio in fitness formula: "
            "fitness = product(areas) * ratio^gamma. "
            "Default: 2 (quadratic). 0 = ignore ratio, 1 = linear."
        ),
    ),
    csv_output: Optional[Path] = typer.Option(
        None,
        "--csv",
        help="CSV output file for fitness scores.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Re-evaluate fitness of existing .dat files without re-running Pepsi-SANS."""
    argv = [
        "fitness_evaluation",
        str(directory),
        "--q-max",
        str(q_max),
        "--ratio-threshold",
        str(ratio_threshold),
        "--gamma",
        str(gamma),
    ]
    if csv_output is not None:
        argv += ["--csv", str(csv_output)]
    if verbose:
        argv.append("--verbose")

    sys.argv = argv
    try:
        from fitness_evaluation import main as eval_main
    except ImportError as exc:
        typer.echo(f"Import error: {exc}", err=True)
        typer.echo(
            "Make sure you are running inside the pixi environment: "
            "pixi run optisans evaluate ...",
            err=True,
        )
        raise typer.Exit(1)
    eval_main()


@app.command()
def recycle(
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Fully protonated input PDB file.",
    ),
    d2o: int = typer.Option(
        ...,
        "--d2o",
        help=(
            "D2O percentage (0–100) for the pattern reference curve. "
            "The SANS curve at this exact D2O value (with the requested AA pattern) "
            "is placed in ref/ as the third fitness reference."
        ),
    ),
    amino_acids: Optional[str] = typer.Option(
        None,
        "-a",
        "--aa",
        help=(
            "Amino acid types to deuterate. Separate multiple codes with spaces or "
            "commas (e.g. --aa 'LEU LYS PRO' or --aa LEU,LYS,PRO). "
            "Run 'optisans aa' to see all available codes."
        ),
    ),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        help=(
            "Base output directory. "
            "Default: '<pdb_stem>_recycle/' in the current folder. "
            "Contains '<stem>_recycle_deuterated_pdbs/' (PDB files) and "
            "'<stem>_recycle_primus_out/' (SANS results and plots)."
        ),
    ),
    step: int = typer.Option(
        1,
        "--step",
        help="Step between D2O values in the scan (1 = every percent point, default).",
    ),
    batch_script: Path = typer.Option(
        Path("./parallel_process_pdb.sh"),
        "--batch-script",
        help="Path to parallel_process_pdb.sh.",
    ),
    conc: float = typer.Option(
        2.5,
        "--conc",
        help="Pepsi-SANS concentration passed as --conc (default: 2.5).",
    ),
    q_max: float = typer.Option(
        0.3,
        "--q-max",
        help="Maximum q value for fitness evaluation (Å⁻¹, default 0.3).",
    ),
    ratio_threshold: float = typer.Option(
        0.01,
        "--ratio-threshold",
        help="Minimum Imax/background ratio to accept a curve (default 0.01).",
    ),
    gamma: float = typer.Option(
        2,
        "--gamma",
        help=(
            "Exponent for Imax/background ratio in fitness formula: "
            "fitness = product(areas) * ratio^gamma. "
            "Default: 2 (quadratic). 0 = ignore ratio, 1 = linear."
        ),
    ),
    n_jobs: int = typer.Option(
        150,
        "--jobs",
        "-j",
        help="Number of parallel Pepsi-SANS jobs (default 150).",
    ),
    no_default_ref: bool = typer.Option(
        False,
        "--no-default-ref",
        help=(
            "Do not generate the default protonated-in-D2O / H2O reference PDBs. "
            "Use this only if reference .dat files are already present in the "
            "ref/ subfolder of the SANS output directory."
        ),
    ),
):
    """
    Contrast variation scan: generate all D2O PDBs (0–100%) with a fixed
    deuteration pattern, simulate SANS, evaluate fitness, and produce plots.

    Outputs (inside {pdb_stem}_recycle/ unless --output-dir is set):
      {stem}_recycle_deuterated_pdbs/   — all deuterated PDB files
      {stem}_recycle_primus_out/        — Pepsi-SANS .dat files + result.csv + plots
      {stem}_recycle_primus_out/ref/    — 3 reference SANS curves used for fitness
    """
    import re as _re

    aa_list: List[str] = []
    if amino_acids:
        aa_list = [
            x.strip().upper() for x in _re.split(r"[,\s]+", amino_acids) if x.strip()
        ]

    # Validate AA codes early, before any expensive operations
    invalid = [aa for aa in aa_list if aa not in VALID_AA]
    if invalid:
        typer.echo(f"Error: invalid amino acid code(s): {', '.join(invalid)}", err=True)
        typer.echo(f"Valid codes: {', '.join(sorted(VALID_AA))}", err=True)
        raise typer.Exit(1)

    try:
        from recycle import recycle_workflow
    except ImportError as exc:
        typer.echo(f"Import error: {exc}", err=True)
        typer.echo(
            "Make sure you are running inside the pixi environment: "
            "pixi run optisans recycle ...",
            err=True,
        )
        raise typer.Exit(1)

    try:
        recycle_workflow(
            pdb_file=str(pdb_file),
            d2o_ref=d2o,
            amino_acids=aa_list,
            output_dir=str(output_dir) if output_dir is not None else None,
            batch_script=str(batch_script),
            d2o_step=step,
            q_max=q_max,
            ratio_threshold=ratio_threshold,
            n_jobs=n_jobs,
            no_default_ref=no_default_ref,
            gamma=gamma,
            conc=conc,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1)


@app.command()
def simulate(
    source: Path = typer.Argument(
        ...,
        exists=True,
        help="Input PDB file or directory containing PDB files.",
    ),
    batch_script: Path = typer.Option(
        Path("./parallel_process_pdb.sh"),
        "--batch-script",
        help="Path to parallel_process_pdb.sh.",
    ),
    jobs: int = typer.Option(
        150,
        "--jobs",
        "-j",
        help="Number of parallel Pepsi-SANS jobs (default 150).",
    ),
    conc: float = typer.Option(
        2.5,
        "--conc",
        help="Pepsi-SANS concentration passed as --conc (default: 2.5).",
    ),
    d2o: Optional[float] = typer.Option(
        None,
        "--d2o",
        help="Explicit D2O value for all files. If omitted, extracted from filenames.",
    ),
):
    """Run Pepsi-SANS simulations on one or more PDB files.

    D2O values are extracted from filenames by default:
      *_d2oXX.pdb -> --d2o XX/100
      *_total_deuteration.pdb -> --d2o 1
      *_total_protonation.pdb -> --d2o 0
    """
    if not batch_script.exists():
        typer.echo(f"Error: batch script not found: {batch_script}", err=True)
        raise typer.Exit(1)

    tmp_list_file = None
    try:
        if source.is_dir():
            pdb_dir = source.resolve()
            cmd = [str(batch_script.resolve()), str(pdb_dir), str(jobs), "--conc", str(conc)]
        else:
            pdb_dir = source.parent.resolve()
            fd, tmp_list_file = tempfile.mkstemp(suffix=".txt", prefix="pdb_list_")
            with os.fdopen(fd, "w") as fh:
                fh.write(str(source.resolve()) + "\n")
            cmd = [str(batch_script.resolve()), str(pdb_dir), tmp_list_file, str(jobs), "--conc", str(conc)]

        if d2o is not None:
            cmd += ["--d2o", str(d2o)]

        typer.echo(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=str(batch_script.parent.resolve()))
        if result.returncode != 0:
            raise typer.Exit(result.returncode)
    finally:
        if tmp_list_file and os.path.exists(tmp_list_file):
            os.unlink(tmp_list_file)


def main():
    app()


if __name__ == "__main__":
    main()
