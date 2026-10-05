"""Labile H classification from the bonded heavy atom element.

An exchangeable hydrogen (or deuterium) is bonded to O, N or S. The bonded
partner comes from the coordinate-based topology, so the classification does
not depend on the PDB atom ordering.
"""

from __future__ import annotations

from typing import Dict, Optional

import gemmi

from .topology import AtomKey, Topology

HETERO_ATOMS = {"O", "N", "S"}


def is_labile(atom: gemmi.Atom, partner: Optional[gemmi.Atom]) -> bool:
    """True when atom is an exchangeable H/D, i.e. bonded to O, N or S."""
    if partner is None:
        return False
    return partner.element.name in HETERO_ATOMS


def classify(model: gemmi.Model, topo: Topology) -> Dict[AtomKey, bool]:
    """Map AtomKey -> labile flag for every H/D atom of the model."""
    labile_map: Dict[AtomKey, bool] = {}
    for ci, chain in enumerate(model):
        for ri, residue in enumerate(chain):
            for atom in residue:
                if atom.is_hydrogen():
                    key = (ci, ri, atom.name.strip())
                    partner = topo.partners.get(key)
                    labile_map[key] = is_labile(
                        atom, partner.heavy_atom if partner else None
                    )
    return labile_map
