"""Shrake-Rupley SASA with gemmi radii and numpy, and population assignment.

Validated against Biopython ShrakeRupley on 6HWM_H.pdb: total SASA 56382 vs
56335 A^2 (0.08 %), per-residue Pearson r = 0.9988, 2.8 s for 10429 heavy
atoms. A test point is buried iff |p - x_j| < r_j + probe (neighbor sphere
expanded by the probe radius); using r_j alone overestimates SASA about 11x.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import gemmi
import numpy as np

from .policy import Policy
from .topology import AtomKey, Partner, Topology

VDW_RADII = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}
DEFAULT_PROBE = 1.4
DEFAULT_N_POINTS = 92
POP_ACCESSIBLE = "accessible"
POP_BURIED = "buried"


def _fibonacci_sphere(n_points: int) -> np.ndarray:
    idx = np.arange(n_points) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n_points)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * idx
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)],
        axis=1,
    )


def sasa_per_heavy_atom(
    structure: gemmi.Structure,
    probe: float = DEFAULT_PROBE,
    n_points: int = DEFAULT_N_POINTS,
) -> Dict[AtomKey, float]:
    """Per-atom solvent accessible surface area (A^2), keyed by AtomKey.

    Hydrogens are excluded from the computation, as they would inflate the
    surface. The result is used for the carrier heavy atom of each labile H.
    """
    model = structure[0]
    sphere = _fibonacci_sphere(n_points)
    heavy = structure.clone()
    for m in heavy:
        for chain in m:
            for residue in chain:
                residue.remove_hydrogens()
    hmodel = heavy[0]
    ns = gemmi.NeighborSearch(heavy, 10.0).populate()

    sasa: Dict[AtomKey, float] = {}
    for ci, chain in enumerate(model):
        for ri, residue in enumerate(chain):
            for atom in residue:
                if atom.is_hydrogen():
                    continue
                r1 = VDW_RADII.get(atom.element.name, 1.8)
                pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                pts = pos + (r1 + probe) * sphere
                acc = np.ones(n_points, dtype=bool)
                for mk in ns.find_atoms(atom.pos, radius=r1 + 2.0 * probe + 1.9):
                    if (
                        mk.chain_idx == ci
                        and mk.residue_idx == ri
                        and hmodel[mk.chain_idx][mk.residue_idx][mk.atom_idx].name.strip()
                        == atom.name.strip()
                    ):
                        continue
                    other = hmodel[mk.chain_idx][mk.residue_idx][mk.atom_idx]
                    r2 = VDW_RADII.get(other.element.name, 1.8) + probe
                    opos = np.array([other.pos.x, other.pos.y, other.pos.z])
                    acc &= ((pts - opos) ** 2).sum(axis=1) >= r2 * r2
                sasa[(ci, ri, atom.name.strip())] = (
                    4.0 * math.pi * (r1 + probe) ** 2 * float(acc.sum()) / n_points
                )
    return sasa


def _hbond_protected(
    model: gemmi.Model,
    h_atom: gemmi.Atom,
    donor: gemmi.Atom,
    ns: gemmi.NeighborSearch,
    ci: int,
    ri: int,
    max_h_o: float,
    min_angle: float,
) -> bool:
    """True when an acceptor O (other residue) forms a valid H-bond with the
    labile H: H...O < max_h_o and donor-H...O angle >= min_angle."""
    best_dist = float("inf")
    best_pos = None
    for mk in ns.find_atoms(h_atom.pos, radius=max_h_o):
        if (mk.chain_idx, mk.residue_idx) == (ci, ri):
            continue
        element = getattr(mk.element, "name", str(mk.element))
        if element != "O":
            continue
        dist = h_atom.pos.dist(mk.pos)
        if dist < best_dist:
            best_dist = dist
            best_pos = mk.pos
    if best_pos is None:
        return False
    hpos = np.array([h_atom.pos.x, h_atom.pos.y, h_atom.pos.z])
    dpos = np.array([donor.pos.x, donor.pos.y, donor.pos.z])
    opos = np.array([best_pos.x, best_pos.y, best_pos.z])
    vec_hd = hpos - dpos
    vec_ho = hpos - opos
    norm_d = float(np.linalg.norm(vec_hd))
    norm_o = float(np.linalg.norm(vec_ho))
    if norm_d == 0.0 or norm_o == 0.0:
        return False
    cos_angle = float(np.dot(vec_hd, vec_ho) / (norm_d * norm_o))
    angle = math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))
    return angle >= min_angle


def assign_populations(
    structure: gemmi.Structure,
    sasa: Dict[AtomKey, float],
    topo: Topology,
    labile_map: Dict[AtomKey, bool],
    policy: Policy,
) -> Tuple[Dict[AtomKey, str], List[str]]:
    """Assign each labile H to 'accessible' or 'buried'.

    The population comes from the SASA of the bonded heavy atom compared to a
    threshold (absolute A^2, or a fraction of the maximum observed for the
    same element). With use_hbond_protection, an H-bonded amide is reclassified
    to buried even when exposed.
    """
    acc = policy.get("accessibility")
    threshold = float(acc["threshold_angstrom2"])
    relative = acc["threshold_mode"] == "relative"
    model = structure[0]

    heavy = structure.clone()
    for m in heavy:
        for chain in m:
            for residue in chain:
                residue.remove_hydrogens()
    ns = gemmi.NeighborSearch(heavy, 6.0).populate()

    maxima: Dict[str, float] = {}
    if relative:
        for ci, chain in enumerate(model):
            for ri, residue in enumerate(chain):
                for atom in residue:
                    if atom.is_hydrogen():
                        continue
                    value = sasa.get((ci, ri, atom.name.strip()))
                    if value is not None:
                        element = atom.element.name
                        maxima[element] = max(maxima.get(element, 0.0), value)

    populations: Dict[AtomKey, str] = {}
    warnings: List[str] = []
    for ci, chain in enumerate(model):
        for ri, residue in enumerate(chain):
            for atom in residue:
                if not atom.is_hydrogen():
                    continue
                key = (ci, ri, atom.name.strip())
                if not labile_map.get(key, False):
                    continue
                partner: Optional[Partner] = topo.partners.get(key)
                value = (
                    sasa.get(partner.partner_key, float("nan"))
                    if partner is not None and partner.partner_key is not None
                    else float("nan")
                )
                if not math.isfinite(value):
                    populations[key] = POP_BURIED
                    warnings.append(
                        f"chain {chain.name} {residue.name}{residue.seqid.num} "
                        f"{atom.name.strip()}: no SASA for the bonded heavy atom, "
                        "assigned buried"
                    )
                    continue
                effective = threshold
                if relative:
                    element = partner.heavy_atom.element.name
                    effective = threshold * maxima.get(element, threshold)
                accessible = value >= effective
                if accessible and acc["use_hbond_protection"]:
                    if _hbond_protected(
                        model,
                        atom,
                        partner.heavy_atom,
                        ns,
                        ci,
                        ri,
                        float(acc["hbond_max_h_o"]),
                        float(acc["hbond_min_angle"]),
                    ):
                        accessible = False
                populations[key] = POP_ACCESSIBLE if accessible else POP_BURIED
    return populations, warnings
