"""Bonded-partner detection with gemmi NeighborSearch on a heavy-atom-only copy.

Each hydrogen (or deuterium) is assigned the nearest heavy atom within a
cutoff; the element of that partner (O/N/S vs C) drives lability. The result
is independent of the PDB atom ordering, unlike the ordering-based heuristic
of the legacy OptiSANS code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import gemmi

from .policy import Policy

AtomKey = Tuple[int, int, str]  # (chain index, residue index, atom name)


@dataclass
class Partner:
    heavy_atom: Optional[gemmi.Atom] = None
    partner_key: Optional[AtomKey] = None
    distance: float = float("nan")
    ambiguous: bool = False
    cross_residue: bool = False


@dataclass
class Topology:
    partners: Dict[AtomKey, Partner] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def check_quarantine(self, quarantine: str) -> None:
        if quarantine == "error" and self.warnings:
            raise ValueError("topology quarantine: " + "; ".join(self.warnings[:5]))


def build(structure: gemmi.Structure, policy: Policy) -> Topology:
    """Build the bonded-partner map for every H/D atom of structure[0]."""
    cutoff = policy.get("lability", "cutoff_angstrom")
    tolerance = policy.get("lability", "ambiguity_tolerance")
    model = structure[0]

    heavy = structure.clone()
    for m in heavy:
        for chain in m:
            for residue in chain:
                residue.remove_hydrogens()
    hmodel = heavy[0]
    ns = gemmi.NeighborSearch(heavy, 5.0).populate()

    topo = Topology()
    for ci, chain in enumerate(model):
        for ri, residue in enumerate(chain):
            for atom in residue:
                if not atom.is_hydrogen():
                    continue
                key = (ci, ri, atom.name.strip())
                mk = ns.find_nearest_atom(atom.pos, cutoff)
                if mk is None or mk.atom_idx is None:
                    topo.partners[key] = Partner()
                    topo.warnings.append(
                        f"chain {chain.name} {residue.name}{residue.seqid.num} {atom.name.strip()}: "
                        f"no heavy partner within {cutoff} A"
                    )
                    continue
                heavy_atom = hmodel[mk.chain_idx][mk.residue_idx][mk.atom_idx]
                dist = atom.pos.dist(heavy_atom.pos)
                cross = (mk.chain_idx, mk.residue_idx) != (ci, ri)
                amb = False
                if not cross:
                    second = None
                    for mk2 in ns.find_atoms(atom.pos, radius=cutoff):
                        if (mk2.chain_idx, mk2.residue_idx, mk2.atom_idx) == (
                            mk.chain_idx,
                            mk.residue_idx,
                            mk.atom_idx,
                        ):
                            continue
                        d2 = atom.pos.dist(mk2.pos)
                        if second is None or d2 < second:
                            second = d2
                    amb = second is not None and (second - dist) < tolerance
                if cross:
                    topo.warnings.append(
                        f"chain {chain.name} {residue.name}{residue.seqid.num} {atom.name.strip()}: "
                        f"nearest heavy {heavy_atom.name.strip()} is in another residue "
                        f"(d={dist:.2f} A)"
                    )
                if amb:
                    topo.warnings.append(
                        f"chain {chain.name} {residue.name}{residue.seqid.num} {atom.name.strip()}: "
                        f"ambiguous partner (two heavy atoms within {tolerance} A)"
                    )
                topo.partners[key] = Partner(
                    heavy_atom=heavy_atom,
                    partner_key=(mk.chain_idx, mk.residue_idx, heavy_atom.name.strip()),
                    distance=dist,
                    ambiguous=amb,
                    cross_residue=cross,
                )
    return topo
