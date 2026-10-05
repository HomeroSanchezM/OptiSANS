"""Selection of residues and atoms for non-labile deuteration.

The package has no knowledge of the OptiSANS 18-group chromosome; callers
translate their pattern into residue names or an explicit atom mask.
An empty selection is allowed: it produces labile-only deuteration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import AbstractSet, FrozenSet, Tuple

AtomKey = Tuple[str, int, str]  # (chain id, residue number, atom name)


@dataclass(frozen=True)
class Selection:
    residue_names: FrozenSet[str] = frozenset()
    atom_mask: FrozenSet[AtomKey] = frozenset()
    all_residues: bool = False

    @classmethod
    def by_residue_names(cls, names: AbstractSet[str]) -> "Selection":
        return cls(residue_names=frozenset(n.strip().upper() for n in names if n.strip()))

    @classmethod
    def by_atom_mask(cls, mask: AbstractSet[AtomKey]) -> "Selection":
        return cls(atom_mask=frozenset(mask))

    @classmethod
    def all(cls) -> "Selection":
        return cls(all_residues=True)

    def is_empty(self) -> bool:
        return not self.residue_names and not self.atom_mask and not self.all_residues

    def selects_residue(self, residue_name: str) -> bool:
        return self.all_residues or residue_name.strip().upper() in self.residue_names

    def selects_atom(self, chain_id: str, residue_number: int, atom_name: str) -> bool:
        return (chain_id, residue_number, atom_name.strip().upper()) in self.atom_mask
