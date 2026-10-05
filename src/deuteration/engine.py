"""Thin helpers around the vendored DeuterPDB engine, expressed in OptiSANS terms.

Everything here is about *translation*:

- 18-gene chromosome -> :class:`Selection` (linked pairs stay linked),
- run seed -> stable per-chromosome seed (same PDBs from one run to the next),
- :class:`Template` -> the legacy composition statistics used by the plots.

The engine itself (labile classification, exact-count sampling, PDB writing)
lives in ``src/deuterpdb`` and is not modified here.
"""

from __future__ import annotations

import hashlib
import logging
from typing import Any, Dict, List, Optional, Sequence

from src.aa import AMINO_ACIDS, expand_linked_codes
from src.deuterpdb import (
    DeuterationResult,
    DeuterationSpec,
    Policy,
    Selection,
    Template,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Template
# ---------------------------------------------------------------------------

def build_template(pdb_path: str, policy: Policy) -> Template:
    """Load the PDB once (protonation and strip per policy) and reuse it.

    The template holds the parsed structure, the lability map and the seeds'
    hash; it is immutable, so a single instance serves every chromosome of a
    whole GA run.
    """
    return Template.from_pdb(pdb_path, policy)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def selection_from_codes(codes: Sequence[str], warn: bool = True) -> Selection:
    """Selection from 3-letter codes; linked pairs are expanded (ASN => ASP)."""
    expanded, pairs = expand_linked_codes(codes)
    if warn:
        for given, partner in pairs:
            logger.info(
                f"  {given} est lié à {partner} (paire ASN+ASP / GLU+GLN) : "
                f"les deux seront deutérés"
            )
    return Selection.by_residue_names(set(expanded))


def selection_from_chromosome(bits18: Sequence[bool]) -> Selection:
    """Selection from an 18-gene chromosome (linked pairs share one gene)."""
    from src.aa import codes_from_chromosome

    return Selection.by_residue_names(set(codes_from_chromosome(bits18)))


# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------

def seed_for_spec(run_seed: int, bits18: Sequence[bool], d2o: float) -> int:
    """Deterministic per-spec seed derived from the run seed and the spec.

    Two runs with the same seed, pattern and D2O produce byte-identical PDBs;
    two different chromosomes never share the same sampling sequence.
    """
    bits = "".join("1" if bit else "0" for bit in bits18)
    payload = f"{int(run_seed)}|{bits}|{int(d2o)}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16)


def seed_for_codes(run_seed: int, codes: Sequence[str], d2o: float) -> int:
    """Same derivation as :func:`seed_for_spec`, for a list of residue codes."""
    payload = f"{int(run_seed)}|{','.join(sorted(code.upper() for code in codes))}|{int(d2o)}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16)


# ---------------------------------------------------------------------------
# Deuteration
# ---------------------------------------------------------------------------

def deuterate_chromosome(
    template: Template,
    bits18: Sequence[bool],
    d2o: float,
    run_seed: int,
    collect_trace: bool = False,
) -> DeuterationResult:
    """Deuterate one GA chromosome (simple mode, exact-count, seeded)."""
    spec = DeuterationSpec(
        selection=selection_from_chromosome(bits18),
        mode="simple",
        d2o_percent=float(d2o),
        seed=seed_for_spec(run_seed, bits18, d2o),
        collect_trace=collect_trace,
    )
    return template.deuterate(spec)


def deuterate_codes(
    template: Template,
    codes: Optional[Sequence[str]] = None,
    all_residues: bool = False,
    d2o: Optional[float] = None,
    seed: Optional[int] = None,
    mode: Optional[str] = None,
    collect_trace: bool = True,
) -> DeuterationResult:
    """Deuterate an explicit selection (``optisans deuterate`` and references).

    ``mode=None`` lets the policy decide (simple by default).  An empty
    selection with ``all_residues=False`` gives labile-only deuteration, which
    is exactly what the internal 0 % / 100 % D2O references need.
    """
    if all_residues:
        selection = Selection.all()
    else:
        selection = selection_from_codes(codes or [])
    spec = DeuterationSpec(
        selection=selection,
        mode=mode,
        d2o_percent=None if d2o is None else float(d2o),
        seed=seed,
        collect_trace=collect_trace,
    )
    return template.deuterate(spec)


# ---------------------------------------------------------------------------
# Composition (legacy-compatible statistics for the plots)
# ---------------------------------------------------------------------------

def composition_stats(template: Template) -> Dict[str, Any]:
    """Residue counts, H counts and absent types over the 20 standard AAs.

    Mirrors the attributes the OptiSANS plotting code used to read from the
    old ``PdbDeuteration`` object.
    """
    comp = template.composition()
    counts: Dict[str, int] = comp["residue_counts"]
    h_counts: Dict[str, int] = comp["hydrogen_counts"]
    aa_count = {aa.code_3: counts.get(aa.code_3, 0) for aa in AMINO_ACIDS}
    aa_hydrogen_count = {aa.code_3: h_counts.get(aa.code_3, 0) for aa in AMINO_ACIDS}
    missing_aa = [aa.code_3 for aa in AMINO_ACIDS if aa_count[aa.code_3] == 0]
    return {
        "aa_count": aa_count,
        "aa_hydrogen_count": aa_hydrogen_count,
        "missing_aa": missing_aa,
        "n_residues": comp["n_residues"],
    }


def apply_stats_to_chromosome(chrom, result: DeuterationResult) -> None:
    """Copy the engine counters into the legacy chromosome fields (same keys)."""
    stats = result.stats
    chrom.H = stats["hydrogen_atoms"]
    chrom.D = stats["deuterium_atoms"]
    chrom.non_labile_D = stats["non_labile_D"]


__all__: List[str] = [
    "build_template",
    "selection_from_codes",
    "selection_from_chromosome",
    "seed_for_spec",
    "seed_for_codes",
    "deuterate_chromosome",
    "deuterate_codes",
    "composition_stats",
    "apply_stats_to_chromosome",
]
