#!/usr/bin/env python3
"""OptiSANS — Unified CLI for protein deuteration optimization for SANS."""

import sys
import subprocess
from pathlib import Path
from typing import Optional, List

# Ensure sibling modules in src/ are importable when running as an entry point.
_src_dir = str(Path(__file__).resolve().parent)
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

import typer

app = typer.Typer(
    name="optisans",
    help="Protein deuteration optimization for SANS experiments.",
    
    add_completion=False,
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
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@app.command()
def run(
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Source PDB file (all hydrogens must be explicit and protonated).",
    ),
    population_size: Optional[int] = typer.Option(
        None, "-p", "--population-size",
        help="Population size (must be a multiple of 3).",
    ),
    generations: Optional[int] = typer.Option(
        None, "-g", "--generations",
        help="Maximum number of generations to run.",
    ),
    elitism: Optional[int] = typer.Option(
        None, "-e", "--elitism",
        help="Number of elite individuals preserved (must be <= population_size / 3).",
    ),
    d2o_var: Optional[int] = typer.Option(
        None, "--d2o-var",
        help="Maximum D2O variation per mutation (0-100).",
    ),
    seed: Optional[int] = typer.Option(
        None, "--seed",
        help="Random seed for reproducibility.",
    ),
    patience: Optional[int] = typer.Option(
        None, "--patience",
        help=(
            "Early stopping: number of consecutive generations without fitness "
            "improvement before stopping. Default: 50. Set to 0 to disable."
        ),
    ),
    config: Optional[Path] = typer.Option(
        None, "--config",
        exists=True,
        help="config.ini file (CLI arguments override INI values).",
    ),
    output_dir: Optional[Path] = typer.Option(
        None, "--output-dir",
        help="Output directory.",
    ),
    batch_script: Optional[Path] = typer.Option(
        None, "--batch-script",
        help="Path to parallel_process_pdb.sh.",
    ),
    q_max: Optional[float] = typer.Option(
        None, "--q-max",
        help="Maximum q value for fitness evaluation (A^-1).",
    ),
    ratio_threshold: Optional[float] = typer.Option(
        None, "--ratio-threshold",
        help="Minimum Imax/background ratio to accept a curve.",
    ),
    d2o_values: Optional[List[int]] = typer.Option(
        None, "--d2o",
        help="Lock D2O to fixed values (repeat flag, e.g. --d2o 0 --d2o 42 --d2o 100).",
    ),
    no_default_ref: bool = typer.Option(
        False, "--no-default-ref",
        help="Do not create the default protonated-in-D2O / H2O reference PDBs.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Run the full genetic algorithm on a PDB file."""
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
    if q_max is not None:
        argv += ["--q-max", str(q_max)]
    if ratio_threshold is not None:
        argv += ["--ratio-threshold", str(ratio_threshold)]
    if d2o_values:
        argv += ["--d2o"] + [str(v) for v in d2o_values]
    if no_default_ref:
        argv += ["--no_default_ref"]
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
    pdb_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Input PDB file.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output",
        help="Output PDB file.",
    ),
    d2o: float = typer.Option(
        0.0, "--d2o",
        help="D2O percentage for labile hydrogen exchange (0-100).",
    ),
    amino_acids: Optional[List[str]] = typer.Option(
        None, "-a", "--aa",
        help="Amino acid types to deuterate (3-letter codes, e.g. --aa ALA --aa GLY).",
    ),
    all_aa: bool = typer.Option(
        False, "--all",
        help="Deuterate all amino acid types.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Deuterate a single PDB file according to a given specification."""
    argv = [
        "pdb_deuteration",
        "-i", str(pdb_file),
        "-o", str(output),
        "--d2o", str(d2o),
    ]
    if all_aa:
        argv.append("--all")
    elif amino_acids:
        for aa in amino_acids:
            aa_upper = aa.upper()
            if aa_upper not in VALID_AA:
                typer.echo(f"Error: invalid amino acid code: {aa}", err=True)
                typer.echo(
                    f"Valid codes: {', '.join(sorted(VALID_AA))}", err=True,
                )
                raise typer.Exit(1)
            argv.append(f"--{aa_upper}")
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
def batch(
    pdb_files: List[Path] = typer.Argument(
        ...,
        exists=True,
        help="PDB files to process (one or more).",
    ),
    config: Optional[Path] = typer.Option(
        None, "--config",
        exists=True,
        help="config.ini file for GA parameters.",
    ),
    batch_script: Path = typer.Option(
        Path("run_convergence_simulation_multiprotein.sh"),
        "--batch-script",
        help="Path to run_convergence_simulation_multiprotein.sh.",
    ),
):
    """Run GA simulations across multiple proteins (convergence study).

    Seeds are defined inside run_convergence_simulation_multiprotein.sh.
    To use different seeds, edit the SEEDS variable in that script or pass
    --config with an appropriate config.ini file.
    """
    for pdb in pdb_files:
        if not pdb.exists():
            typer.echo(f"Error: PDB file not found: {pdb}", err=True)
            raise typer.Exit(1)

    if not batch_script.exists():
        typer.echo(f"Error: batch script not found: {batch_script}", err=True)
        raise typer.Exit(1)

    cmd = ["bash", str(batch_script)]
    if config is not None:
        cmd += ["--config", str(config)]
    cmd += [str(p) for p in pdb_files]

    typer.echo(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise typer.Exit(result.returncode)


@app.command()
def evaluate(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=True,
        help="Directory containing .dat files and a ref/ subfolder.",
    ),
    q_max: float = typer.Option(
        0.3, "--q-max",
        help="Maximum q value for truncation (A^-1).",
    ),
    ratio_threshold: float = typer.Option(
        0.01, "--ratio-threshold",
        help="Minimum Imax/background ratio to accept a curve.",
    ),
    csv_output: Optional[Path] = typer.Option(
        None, "--csv",
        help="CSV output file for fitness scores.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose",
        help="Enable verbose logging.",
    ),
):
    """Re-evaluate fitness of existing .dat files without re-running Pepsi-SANS."""
    argv = [
        "fitness_evaluation",
        str(directory),
        "--q-max", str(q_max),
        "--ratio-threshold", str(ratio_threshold),
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
def plot(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=True,
        help="Results directory (must contain best_fitness_summary.csv and "
        "generation_XX_summary.txt files).",
    ),
    output: Optional[Path] = typer.Option(
        None, "-o", "--output",
        help="Output path for the fitness plot "
        "(default: Fitness_evolution.png inside the directory).",
    ),
    annotate: bool = typer.Option(
        False, "--annotate",
        help="Add a colour grid of deuterated AAs below the fitness plot.",
    ),
    show_min: bool = typer.Option(
        False, "--min",
        help="With --annotate: show stat values only on change and at last column.",
    ),
    fitness_only: bool = typer.Option(
        False, "--fitness-only",
        help="Generate only the fitness evolution plot (skip D2O scatter plots).",
    ),
    scatter_only: bool = typer.Option(
        False, "--scatter-only",
        help="Generate only the D2O vs %D scatter plots (skip fitness plot).",
    ),
    interactive: bool = typer.Option(
        False, "--interactive",
        help="Display the fitness plot interactively (matplotlib show).",
    ),
):
    """Generate result plots: fitness evolution and D2O vs %D scatter plots."""
    if fitness_only and scatter_only:
        typer.echo(
            "Error: --fitness-only and --scatter-only are mutually exclusive.",
            err=True,
        )
        raise typer.Exit(1)

    # Plot 1: fitness evolution
    if not scatter_only:
        csv_file = directory / "best_fitness_summary.csv"
        if not csv_file.exists():
            typer.echo(
                f"Warning: best_fitness_summary.csv not found in {directory}",
                err=True,
            )
        else:
            argv = ["plot_fitness_evolution", str(csv_file)]
            if output is not None:
                argv += ["-o", str(output)]
            if annotate:
                argv.append("--annotate")
            if show_min:
                argv.append("--min")
            if interactive:
                argv.append("--interactive")
            sys.argv = argv
            try:
                from plot_fitness_evolution import main as plot_main
                plot_main()
            except ImportError as exc:
                typer.echo(f"Import error: {exc}", err=True)
                typer.echo(
                    "Make sure you are running inside the pixi environment: "
                    "pixi run optisans plot ...",
                    err=True,
                )
                raise typer.Exit(1)

    # Plot 2: D2O vs %D scatter
    if not fitness_only:
        sys.argv = ["d2o_vs_d", str(directory)]
        try:
            from d2o_vs_d import main as scatter_main
            scatter_main()
        except ImportError as exc:
            typer.echo(f"Import error: {exc}", err=True)
            typer.echo(
                "Make sure you are running inside the pixi environment: "
                "pixi run optisans plot ...",
                err=True,
            )
            raise typer.Exit(1)


def main():
    app()


if __name__ == "__main__":
    main()
