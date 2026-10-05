"""Deuteration adapter between OptiSANS and the vendored DeuterPDB engine.

Responsibilities:

- build a :class:`Policy` from the OptiSANS INI files (no policy YAML needed),
- translate an 18-gene chromosome (or a list of residue codes) into a
  :class:`Selection`, keeping ASN+ASP and GLU+GLN linked,
- derive a stable per-chromosome seed from a run seed,
- build the reusable :class:`Template` once per run and deuterate chromosomes.
"""

from .engine import (
    apply_stats_to_chromosome,
    build_template,
    composition_stats,
    deuterate_chromosome,
    deuterate_codes,
    seed_for_codes,
    seed_for_spec,
    selection_from_chromosome,
    selection_from_codes,
)
from .policy_ini import (
    DeuterPDBConfigError,
    load_policy,
    load_policy_and_options,
    read_deuteration_ini,
)

__all__ = [
    "apply_stats_to_chromosome",
    "build_template",
    "composition_stats",
    "deuterate_chromosome",
    "deuterate_codes",
    "seed_for_codes",
    "seed_for_spec",
    "selection_from_chromosome",
    "selection_from_codes",
    "DeuterPDBConfigError",
    "load_policy",
    "load_policy_and_options",
    "read_deuteration_ini",
]
