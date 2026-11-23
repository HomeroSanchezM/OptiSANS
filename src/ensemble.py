"""Ensemble (multi-PDB) support for ``optisans run``.

A run accepts either a single PDB file or a directory of PDB files
(conformations).  In ensemble mode (input directory, possibly truncated by
``--conf N``) every chromosome is deuterated once per conformation
(``_confNNNN`` filename suffix), simulated by Pepsi-SANS, and its fitness
is computed on the mean of the per-conformation SANS curves.

Single-PDB mode is the degenerate case of one conformation: filenames keep
their legacy format (no ``_confNNNN`` suffix) and the existing GA loop,
summaries and cleanup are reused unchanged.
"""

from __future__ import annotations

import logging
import random
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from src.GA import Chromosome
from src.deuteration import (
    build_template,
    composition_stats,
    deuterate_chromosome,
)

logger = logging.getLogger(__name__)

# Suffix appended to chromosome filenames in ensemble mode: _confNNNN
CONF_SUFFIX_RE = re.compile(r"_conf\d{4}$")


# ---------------------------------------------------------------------------
# Input resolution
# ---------------------------------------------------------------------------

def resolve_pdb_inputs(path, n_conf: Optional[int] = None,
                       seed: Optional[int] = None) -> List[Path]:
    """Resolve the positional input to a list of PDB paths.

    file -> ``[file]`` (``n_conf`` ignored with a warning, unless it is 1);
    dir  -> sorted ``*.pdb``, optionally seeded-shuffled, truncated to N.

    Raises ValueError on a missing path, an empty directory, or N out of
    range (the error message reports how many PDBs were found).
    """
    p = Path(path)
    if not p.exists():
        raise ValueError(f"Input path not found: {p}")
    if p.is_file():
        if n_conf not in (None, 1):
            logger.warning("--conf ignored: input is a single PDB file")
        return [p]
    pdbs = sorted(p.glob("*.pdb"))
    if not pdbs:
        raise ValueError(f"No .pdb files found in directory: {p} (glob '*.pdb')")
    if n_conf is not None:
        if not (1 <= n_conf <= len(pdbs)):
            raise ValueError(
                f"--conf {n_conf} out of range: {len(pdbs)} PDB file(s) found in {p}"
            )
        if seed is not None:
            random.Random(seed).shuffle(pdbs)
        pdbs = pdbs[:n_conf]
    return pdbs


# ---------------------------------------------------------------------------
# Filename helpers
# ---------------------------------------------------------------------------

def get_pdb_filename(chrom: Chromosome, conf_idx: Optional[int] = None) -> str:
    """Build the PDB filename of a chromosome.

    Format: ``gen{g:02d}_Chr{i:03d}_d2o{d:03d}_deut{n:02d}[_conf{c:04d}].pdb``

    The ``_confNNNN`` suffix is added only in ensemble mode, so single-PDB
    filenames are unchanged and the ``_d2o(\\d+)`` regexes of
    ``parallel_process_pdb.sh`` keep matching conf-suffixed names.
    """
    num_deut = sum(chrom.deuteration)
    stem = (f"gen{chrom.generation:02d}_Chr{chrom.index:03d}"
            f"_d2o{chrom.d2o:03d}_deut{num_deut:02d}")
    if conf_idx is not None:
        stem += f"_conf{conf_idx:04d}"
    return stem + ".pdb"


def chrom_key_from_filename(name: str) -> str:
    """Chromosome-level stem: filename without extension and without ``_confNNNN``.

    Used to match ``.dat`` files back to their chromosome and to clean up
    non-tier-1 files in ensemble mode.
    """
    return CONF_SUFFIX_RE.sub("", Path(name).stem)


# ---------------------------------------------------------------------------
# Templates and composition
# ---------------------------------------------------------------------------

def build_templates(pdb_paths: Sequence[Path], policy) -> list:
    """One reusable Template per input PDB (single mode = list of one)."""
    return [build_template(str(p), policy) for p in pdb_paths]


def aggregate_composition(templates) -> Dict[str, Any]:
    """Sum the composition stats over all conformations.

    ``missing_aa`` is the union over conformations, so a gene is disabled
    only when its amino-acid type is absent from *every* conformation.
    """
    merged: Dict[str, Any] = {
        "aa_count": {},
        "aa_hydrogen_count": {},
        "n_residues": 0,
        "missing_aa": set(),
    }
    for t in templates:
        comp = composition_stats(t)
        for code, n in comp["aa_count"].items():
            merged["aa_count"][code] = merged["aa_count"].get(code, 0) + n
        for code, n in comp["aa_hydrogen_count"].items():
            merged["aa_hydrogen_count"][code] = \
                merged["aa_hydrogen_count"].get(code, 0) + n
        merged["n_residues"] += comp["n_residues"]
        merged["missing_aa"].update(comp["missing_aa"])
    merged["missing_aa"] = sorted(merged["missing_aa"])
    return merged


# ---------------------------------------------------------------------------
# PDB generation (ensemble)
# ---------------------------------------------------------------------------

def generate_pdbs_for_chromosomes_multi(templates, chromosomes,
                                        output_dir: Path, run_seed: int,
                                        collect_trace: bool = False) -> List[str]:
    """Write one deuterated PDB per (conformation, chromosome).

    Single mode (one template) keeps the legacy filenames (no ``_conf``
    suffix).  In ensemble mode each conformation gets its own deterministic
    seed (``conf_idx`` is part of the seed payload in
    :func:`src.deuteration.engine.seed_for_spec`), so conformations sample
    independently but reproducibly.

    H/D stats are summed over conformations into the chromosome, so the %D
    computed from the sums stays a weighted average over the ensemble.
    """
    multi = len(templates) > 1
    generated: List[str] = []
    for chrom in chromosomes:
        h = d = nld = 0
        for conf_idx, template in enumerate(templates):
            suffix = conf_idx if multi else None
            filename = get_pdb_filename(chrom, conf_idx=suffix)
            try:
                result = deuterate_chromosome(
                    template, chrom.deuteration, chrom.d2o, run_seed,
                    conf_idx=suffix, collect_trace=collect_trace,
                )
                stats = result.stats
                h += stats["hydrogen_atoms"]
                d += stats["deuterium_atoms"]
                nld += stats["non_labile_D"]
                out_path = output_dir / filename
                result.write_pdb(str(out_path))
                generated.append(str(out_path))
                logger.debug(f"  Written: {filename}")
            except Exception as exc:
                logger.error(f"  FAILED {filename}: {exc}")
        chrom.H, chrom.D, chrom.non_labile_D = h, d, nld
    logger.info(
        f"  {len(generated)} PDB file(s) written "
        f"({len(chromosomes)} chromosome(s) x {len(templates)} conformation(s))."
    )
    return generated
