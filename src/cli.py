#!/usr/bin/env python3
"""OptiSANS — Unified CLI for protein deuteration optimization for SANS."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import typer

# Single source for the amino-acid tables (see src/aa.py).
from src.aa import AA_DICT, AMINO_ACIDS


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


VALID_AA = {aa.code_3 for aa in AMINO_ACIDS}


@app.command()
def aa():
    """Affiche la liste des 20 acides aminés standards avec leurs codes."""
    typer.echo("\n Amino acids available for deuteration\n")
    typer.echo(f"  {'Code 3':<10} {'Code 1':<10} {'Name'}")
    typer.echo("  " + "─" * 38)
    for aa in AMINO_ACIDS:
        typer.echo(f"  {aa.code_3:<10} {aa.code_1:<10} {aa.name}")
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
    restriction: Optional[List[str]] = typer.Option(
        None,
        "--restriction",
        help=(
            "Restrict (forbid) the deuteration of one or more amino-acid types, "
            "given as 3-letter codes (repeatable). Linked pairs always stay "
            "linked: --restriction ASN also restricts ASP, and a warning is "
            "printed."
        ),
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
    if restriction:
        for code in restriction:
            argv += ["--restriction", code]
    if no_default_ref:
        argv += ["--no_default_ref"]
    if ref:
        argv += ["--ref"] + [str(r) for r in ref]
    if verbose:
        argv += ["--verbose"]

    sys.argv = argv
    try:
        from src.generate_deuterated_pdbs import main as ga_main
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
    source: Optional[Path] = typer.Argument(
        None,
        help="Input PDB file or directory containing PDB files.",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "-o",
        "--output",
        help="Output PDB file, or output directory when SOURCE is a folder.",
    ),
    d2o: Optional[float] = typer.Option(
        None,
        "--d2o",
        help="D2O percentage for labile hydrogen exchange (0-100).",
    ),
    amino_acids: Optional[str] = typer.Option(
        None,
        "-a",
        "--aa",
        help=(
            "Amino acid types to deuterate (3-letter codes, space- or "
            "comma-separated). Linked pairs always stay linked: --aa ASN also "
            "deuterates ASP and a warning is printed. Run 'optisans aa' to see "
            "all available codes."
        ),
    ),
    all_aa: bool = typer.Option(
        False,
        "--all",
        help="Deuterate all amino acid types.",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        help=(
            "INI file with the historical [DEUTERATION] and [AMINO_ACIDS] "
            "sections, used as defaults when the CLI flags are omitted."
        ),
    ),
    config_deuterate: Optional[Path] = typer.Option(
        None,
        "--config-deuterate",
        help=(
            "INI file with the [DEUTERPDB.*] engine options. "
            "Default: ./config_deuterate.ini when present."
        ),
    ),
    policy_file: Optional[Path] = typer.Option(
        None,
        "--policy",
        help=(
            "Explicit DeuterPDB policy YAML (escape hatch; takes precedence "
            "over the [DEUTERPDB.*] sections)."
        ),
    ),
    mode: Optional[str] = typer.Option(
        None, "--mode", help="Labile exchange mode: simple | kinetic | dual."
    ),
    sampling: Optional[str] = typer.Option(
        None, "--sampling", help="Sampling of the labile exchange: exact | bernoulli."
    ),
    rounding: Optional[str] = typer.Option(
        None, "--rounding", help="Exact-count rounding: floor | round | ceil."
    ),
    seed: Optional[int] = typer.Option(
        None, "--seed", help="RNG seed (default: the policy seed, 42)."
    ),
    lability_cutoff: Optional[float] = typer.Option(
        None, "--lability-cutoff", help="Bond cutoff for labile H classification (angstrom)."
    ),
    quarantine: Optional[str] = typer.Option(
        None, "--quarantine", help="Ambiguous labile H handling: warn | error."
    ),
    probability_burried: Optional[float] = typer.Option(
        None, "--probability-burried", help="Buried-H probability (dual mode only)."
    ),
    kinetic_rate: Optional[float] = typer.Option(
        None, "--kinetic-rate", help="Accessible exchange rate in 1/min (kinetic mode)."
    ),
    buried_rate: Optional[float] = typer.Option(
        None, "--buried-rate", help="Buried exchange rate in 1/min (kinetic mode)."
    ),
    time_minutes: Optional[float] = typer.Option(
        None, "--time-minutes", help="Incubation time in minutes (kinetic mode)."
    ),
    accessibility: bool = typer.Option(
        False, "--accessibility", help="Enable SASA accessibility (required by kinetic/dual)."
    ),
    report: Optional[bool] = typer.Option(
        None,
        "--report/--no-report",
        help="Write the per-atom trace + provenance report (default: on).",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Deuterate PDB files with the integrated DeuterPDB engine.

    Non-labile hydrogens are deuterated deterministically for the selected
    residue types (linked pairs ASN+ASP and GLU+GLN always deuterate together);
    labile hydrogens exchange with a probability set by the mode and the D2O
    percentage (simple mode: d2o/100), by exact-count sampling with a seeded RNG.
    """
    import logging as _logging
    import re as _re

    from src.deuterpdb import Policy
    from src.deuteration import build_template, deuterate_codes, load_policy_and_options
    from src.deuteration.policy_ini import DeuterPDBConfigError, read_deuteration_ini

    _logging.basicConfig(
        level=_logging.DEBUG if verbose else _logging.INFO,
        format="%(message)s",
    )

    # ---- configuration files ------------------------------------------------
    config_deuterate_path: Optional[str] = None
    if config_deuterate is not None:
        if not config_deuterate.exists():
            typer.echo(f"Error: config file not found: {config_deuterate}", err=True)
            raise typer.Exit(1)
        config_deuterate_path = str(config_deuterate)
    elif Path("config_deuterate.ini").exists():
        config_deuterate_path = "config_deuterate.ini"

    legacy_path = str(config) if config is not None else config_deuterate_path
    try:
        legacy = read_deuteration_ini(legacy_path)
    except Exception as exc:
        typer.echo(f"Error reading {legacy_path}: {exc}", err=True)
        raise typer.Exit(1)

    # ---- policy: CLI flags > INI > defaults ---------------------------------
    overrides: Dict = {}

    def _put(section: str, key: str, value) -> None:
        if value is None:
            return
        node = overrides.setdefault(section, {})
        parts = key.split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = value

    _put("labile_exchange", "mode", mode)
    _put("labile_exchange", "sampling", sampling)
    _put("labile_exchange", "rounding", rounding)
    _put("labile_exchange", "seed", seed)
    _put("labile_exchange", "probability_burried", probability_burried)
    _put("labile_exchange", "kinetic.accessible_rate", kinetic_rate)
    _put("labile_exchange", "kinetic.buried_rate", buried_rate)
    _put("labile_exchange", "kinetic.time_minutes", time_minutes)
    _put("lability", "cutoff_angstrom", lability_cutoff)
    _put("lability", "quarantine", quarantine)
    if accessibility:
        _put("accessibility", "enabled", True)

    if policy_file is not None:
        try:
            policy = Policy.from_yaml(str(policy_file))
            options = {"write_report": True, "write_pdb": True}
        except Exception as exc:
            typer.echo(f"Error loading policy file {policy_file}: {exc}", err=True)
            raise typer.Exit(1)
    else:
        try:
            policy, options = load_policy_and_options(
                config_deuterate_path,
                overrides=overrides,
                default_write_report=True,
            )
        except (DeuterPDBConfigError, FileNotFoundError, ValueError) as exc:
            typer.echo(f"Error: {exc}", err=True)
            raise typer.Exit(1)

    if report is not None:
        options["write_report"] = report

    # ---- resolve input / output ---------------------------------------------
    src_path: Optional[Path] = Path(source) if source is not None else None
    if src_path is None and legacy.get("input_pdb"):
        src_path = Path(str(legacy["input_pdb"]))
    if src_path is None:
        typer.echo(
            "Error: no input given (SOURCE argument or [DEUTERATION] input_pdb).",
            err=True,
        )
        raise typer.Exit(1)
    if not src_path.exists():
        typer.echo(f"Error: input not found: {src_path}", err=True)
        raise typer.Exit(1)

    d2o_value = float(d2o) if d2o is not None else float(legacy.get("d2o_percent", 0.0))
    if not 0.0 <= d2o_value <= 100.0:
        typer.echo(f"Error: --d2o must be in [0, 100], got {d2o_value}", err=True)
        raise typer.Exit(1)

    # ---- amino-acid selection (linked pairs stay linked) --------------------
    codes: List[str] = []
    all_residues = False
    if all_aa:
        all_residues = True
    else:
        if amino_acids:
            codes = [c.upper() for c in _re.split(r"[,\s]+", amino_acids) if c]
        else:
            codes = [c.upper() for c in legacy.get("amino_acids", [])]
        invalid = [c for c in codes if c not in VALID_AA]
        if invalid:
            typer.echo(
                f"Error: invalid amino acid code(s): {', '.join(invalid)}", err=True
            )
            typer.echo(f"Valid codes: {', '.join(sorted(VALID_AA))}", err=True)
            raise typer.Exit(1)

    policy_seed = policy.get("labile_exchange", "seed")
    run_seed = 42 if policy_seed is None else int(policy_seed)
    if seed is not None:
        run_seed = int(seed)

    typer.echo(
        f"DeuterPDB policy {policy.hash} | mode={policy.get('labile_exchange', 'mode')} "
        f"| sampling={policy.get('labile_exchange', 'sampling')} "
        f"| d2o={d2o_value}% | seed={run_seed} "
        f"| report={'on' if options['write_report'] else 'off'}"
    )

    # ---- processing ---------------------------------------------------------
    def _process_one(in_path: Path, out_path: Path) -> None:
        template = build_template(str(in_path), policy)
        result = deuterate_codes(
            template,
            codes=codes,
            all_residues=all_residues,
            d2o=d2o_value,
            seed=run_seed,
            collect_trace=True,
        )
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        result.write_pdb(str(out_path))
        if options["write_report"]:
            result.write_report(str(out_path.with_suffix("")))
        stats = result.stats
        typer.echo(
            f"  {in_path.name} -> {out_path}  |  D={stats['deuterium_atoms']} "
            f"({stats['deuterium_pct']:.2f}% of H/D)  "
            f"|  non-labile D={stats['non_labile_D']}"
        )

    if src_path.is_dir():
        pdb_files = sorted(
            p for p in src_path.iterdir() if p.is_file() and p.suffix.lower() == ".pdb"
        )
        if not pdb_files:
            typer.echo(f"No .pdb files found in directory: {src_path}", err=True)
            raise typer.Exit(1)
        out_dir = Path(output) if output is not None else src_path
        if out_dir.exists() and not out_dir.is_dir():
            typer.echo(
                "When SOURCE is a directory, --output must be a directory.",
                err=True,
            )
            raise typer.Exit(1)
        out_dir.mkdir(parents=True, exist_ok=True)
        typer.echo(f"Processing {len(pdb_files)} PDB file(s) from {src_path} -> {out_dir}")
        failures = 0
        for in_path in pdb_files:
            try:
                _process_one(in_path, out_dir / in_path.name)
            except Exception as exc:
                typer.echo(f"  FAILED {in_path.name}: {exc}", err=True)
                failures += 1
        if failures:
            raise typer.Exit(1)
        return

    # single file
    if output is None and legacy.get("output_pdb"):
        out_path = Path(str(legacy["output_pdb"]))
    elif output is None:
        out_path = Path("deuterated.pdb")
    else:
        out_path = Path(output)
        if out_path.exists() and out_path.is_dir():
            out_path = out_path / src_path.name

    try:
        _process_one(src_path, out_path)
    except Exception as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1)


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
        from src.fitness_evaluation import main as eval_main
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
    config_deuterate: Optional[Path] = typer.Option(
        None,
        "--config-deuterate",
        help=(
            "INI file with the [DEUTERPDB.*] engine options. "
            "Default: ./config_deuterate.ini when present."
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
        from src.recycle import recycle_workflow
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
            policy_ini=(
                str(config_deuterate)
                if config_deuterate is not None
                else ("config_deuterate.ini" if Path("config_deuterate.ini").exists() else None)
            ),
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
