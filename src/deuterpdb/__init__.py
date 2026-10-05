"""DeuterPDB: PDB parsing, protonation, lability and kinetic deuteration.

Public API:
    Policy              immutable configuration loaded from YAML
    Selection           residue/atom selection for non-labile deuteration
    Template            protonated structure with lability and accessibility maps
    DeuterationSpec     one deuteration run (selection, mode, parameters, seed)
    DeuterationResult   deuterated structure plus stats and per-atom trace
"""

from .deuterate import DeuterationResult, DeuterationSpec, Template
from .policy import Policy
from .select import Selection

__version__ = "0.1.0"

__all__ = [
    "Policy",
    "Selection",
    "Template",
    "DeuterationSpec",
    "DeuterationResult",
    "__version__",
]
