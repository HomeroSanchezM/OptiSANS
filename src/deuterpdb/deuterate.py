"""Deuteration engine: exact-count sampling, H->D conversion, stats and trace."""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import gemmi
import numpy as np

from . import kinetics
from . import lability as lability_mod
from .accessibility import (
    POP_ACCESSIBLE,
    POP_BURIED,
    apply_granularity,
    populations_per_atom,
    residue_populations,
    sasa_per_heavy_atom,
)
from .io import StructureError, load_structure, structure_hash
from .policy import Policy
from .protonate import protonate_structure, strip_water_and_ligand
from .report import write_aggregates, write_atom_trace, write_provenance
from .select import Selection
from .topology import AtomKey, Topology, build as build_topology

STATS_KEYS = [
    "total_atoms",
    "hydrogen_atoms",
    "deuterium_atoms",
    "labile_H",
    "labile_D",
    "non_labile_H",
    "non_labile_D",
]


def _exact_count(n: int, p: float, rounding: str) -> int:
    """Exact-count sampling size, matrix-compatible (floor by default)."""
    if n <= 0 or p <= 0.0:
        return 0
    if p >= 1.0:
        return n
    value = n * p
    if rounding == "floor":
        k = int(math.floor(value))
    elif rounding == "ceil":
        k = int(math.ceil(value))
    else:
        k = int(round(value))
    return max(0, min(n, k))


def _convert_h_to_d(atom: gemmi.Atom) -> None:
    atom.name = atom.name.replace("H", "D", 1)
    atom.element = gemmi.Element("D")


@dataclass(frozen=True)
class DeuterationSpec:
    """One deuteration run on a template. Non-None fields override the policy.

    ``collect_trace`` is not part of :meth:`hash`: turning the per-atom trace
    on or off must not change the derived seed nor the produced PDB.
    """

    selection: Selection = field(default_factory=Selection)
    mode: Optional[str] = None
    d2o_percent: Optional[float] = None
    seed: Optional[int] = None
    collect_trace: bool = True

    def hash(self) -> str:
        text = (
            f"{sorted(self.selection.residue_names)}|{sorted(self.selection.atom_mask)}|"
            f"{self.selection.all_residues}|{self.mode}|{self.d2o_percent}|{self.seed}"
        )
        return hashlib.sha256(text.encode()).hexdigest()[:16]


def _initial_counts(structure: gemmi.Structure, labile_map: Dict[AtomKey, bool]) -> Dict[str, int]:
    counts: Dict[str, int] = {
        "total_atoms": 0,
        "h_total": 0,
        "d_total": 0,
        "h_labile": 0,
        "h_nonlabile": 0,
        "d_labile": 0,
        "d_nonlabile": 0,
    }
    for ci, chain in enumerate(structure[0]):
        for ri, residue in enumerate(chain):
            for atom in residue:
                counts["total_atoms"] += 1
                labile = labile_map.get((ci, ri, atom.name.strip()), False)
                element = atom.element.name
                if element == "H":
                    counts["h_total"] += 1
                    if labile:
                        counts["h_labile"] += 1
                    else:
                        counts["h_nonlabile"] += 1
                elif element == "D":
                    counts["d_total"] += 1
                    if labile:
                        counts["d_labile"] += 1
                    else:
                        counts["d_nonlabile"] += 1
    return counts


def _needs_protonation(structure: gemmi.Structure, policy: Policy) -> bool:
    n_h = sum(
        1
        for model in structure
        for chain in model
        for residue in chain
        for atom in residue
        if atom.element.name == "H"
    )
    if n_h == 0:
        return True
    topo = build_topology(structure, policy)
    orphans = sum(1 for partner in topo.partners.values() if partner.heavy_atom is None)
    return orphans > 0


class Template:
    """Protein-only structure with precomputed lability and population maps.

    Build once per input structure; reuse for any number of DeuterationSpecs.
    `population_map` follows `accessibility.granularity` (residue by default),
    `atom_population_map` keeps the raw per-atom classification and
    `residue_population_map` holds one population per residue.
    """

    def __init__(
        self,
        structure: gemmi.Structure,
        policy: Policy,
        topo: Topology,
        labile_map: Dict[AtomKey, bool],
        population_map: Dict[AtomKey, str],
        dropped: Dict[str, Dict[str, int]],
        warnings: List[str],
        atom_population_map: Optional[Dict[AtomKey, str]] = None,
        residue_population_map: Optional[Dict[Tuple[int, int], str]] = None,
        granularity: str = "residue",
    ):
        self.structure = structure
        self.policy = policy
        self.topo = topo
        self.labile_map = labile_map
        self.population_map = population_map
        self.atom_population_map = atom_population_map or {}
        self.residue_population_map = residue_population_map or {}
        self.granularity = granularity
        self.dropped = dropped
        self.warnings = warnings
        self.hash = structure_hash(structure)
        self._initial = _initial_counts(structure, labile_map)

    @classmethod
    def from_pdb(cls, path: str, policy: Policy) -> "Template":
        structure = load_structure(path, policy)
        warnings_list: List[str] = []
        dropped: Dict[str, Dict[str, int]] = {}
        if policy.get("template", "strip_water_ligand"):
            structure, dropped = strip_water_and_ligand(structure)
        protonation = policy.get("template", "protonation")
        n_h_after_strip = sum(
            1
            for model in structure
            for chain in model
            for residue in chain
            for atom in residue
            if atom.element.name == "H"
        )
        if protonation == "never" and n_h_after_strip == 0:
            raise StructureError(
                f"{path}: no hydrogen atoms and template.protonation is 'never'; "
                "set template.protonation to 'auto' or 'always'"
            )
        if protonation == "always" or (
            protonation == "auto" and _needs_protonation(structure, policy)
        ):
            structure = protonate_structure(structure, policy)
        topo = build_topology(structure, policy)
        topo.check_quarantine(policy.get("lability", "quarantine"))
        warnings_list.extend(topo.warnings)
        labile_map = lability_mod.classify(structure[0], topo)
        population_map: Dict[AtomKey, str] = {}
        atom_population_map: Dict[AtomKey, str] = {}
        residue_population_map: Dict[Tuple[int, int], str] = {}
        granularity = policy.get("accessibility", "granularity")
        if (
            policy.get("labile_exchange", "mode") in ("kinetic", "dual")
            and policy.get("accessibility", "enabled")
        ):
            sasa = sasa_per_heavy_atom(structure)
            atom_population_map, pop_warnings = populations_per_atom(
                structure, sasa, topo, labile_map, policy
            )
            residue_population_map = residue_populations(atom_population_map)
            population_map = apply_granularity(atom_population_map, granularity)
            warnings_list.extend(pop_warnings)
        return cls(
            structure,
            policy,
            topo,
            labile_map,
            population_map,
            dropped,
            warnings_list,
            atom_population_map=atom_population_map,
            residue_population_map=residue_population_map,
            granularity=granularity,
        )

    @property
    def n_labile(self) -> int:
        return sum(1 for value in self.labile_map.values() if value)

    @property
    def n_nonlabile(self) -> int:
        return sum(1 for value in self.labile_map.values() if not value)

    def composition(self) -> Dict[str, Any]:
        """Residue and original-hydrogen counts per residue type.

        Returns a dict::

            {"residue_counts": {name: n}, "hydrogen_counts": {name: n_h},
             "n_residues": total}

        Counts are per *unique* residue (first model, keyed by chain, sequence
        number and name), matching the legacy OptiSANS composition analysis.
        The package stays amino-acid agnostic: callers that know the standard
        20 codes fill in the missing entries themselves.
        """
        residue_counts: Dict[str, int] = {}
        hydrogen_counts: Dict[str, int] = {}
        seen: Set[Tuple[str, str, str]] = set()
        model = self.structure[0]
        for chain in model:
            for residue in chain:
                name = residue.name.strip().upper()
                key = (chain.name, str(residue.seqid), name)
                if key in seen:
                    continue
                seen.add(key)
                residue_counts[name] = residue_counts.get(name, 0) + 1
                n_h = sum(1 for atom in residue if atom.element.name == "H")
                hydrogen_counts[name] = hydrogen_counts.get(name, 0) + n_h
        return {
            "residue_counts": residue_counts,
            "hydrogen_counts": hydrogen_counts,
            "n_residues": sum(residue_counts.values()),
        }

    def deuterate(self, spec: Optional[DeuterationSpec] = None) -> "DeuterationResult":
        spec = spec or DeuterationSpec()
        policy = self.policy
        lx = policy.get("labile_exchange")
        mode = spec.mode or lx["mode"]
        if mode not in ("simple", "kinetic", "dual"):
            raise ValueError(f"unknown labile_exchange.mode {mode!r}")
        if mode == "dual":
            if lx["probability_burried"] is None:
                raise ValueError("labile_exchange.probability_burried is required for dual mode")
            if not policy.get("accessibility", "enabled"):
                raise ValueError("labile_exchange.mode dual requires accessibility.enabled: true")
        d2o = lx["d2o_percent"] if spec.d2o_percent is None else float(spec.d2o_percent)
        if not 0.0 <= d2o <= 100.0:
            raise ValueError("d2o_percent must be in [0, 100]")
        sampling = lx["sampling"]
        rounding = lx["rounding"]
        kin = lx["kinetic"]
        seed = lx["seed"] if spec.seed is None else spec.seed
        if seed is None:
            seed = int(spec.hash() + self.hash, 16) % (2**32)
        rng = np.random.default_rng(seed)

        if mode == "simple":
            p_all = kinetics.simple_probability(d2o)
        elif mode == "kinetic":
            p_acc = kinetics.kinetic_probability(d2o, kin["accessible_rate"], kin["time_minutes"])
            p_bur = kinetics.kinetic_probability(d2o, kin["buried_rate"], kin["time_minutes"])
        else:
            p_acc = kinetics.simple_probability(d2o)
            p_bur = float(lx["probability_burried"])

        model = self.structure[0]

        # Population assignment: template SASA map, or fixed fraction drawn here.
        pop_map: Dict[AtomKey, str] = (
            dict(self.population_map) if mode in ("kinetic", "dual") else {}
        )
        if mode == "kinetic" and not pop_map and kin["buried_fraction_source"] == "fixed":
            labile_keys = [key for key, labile in self.labile_map.items() if labile]
            n_buried = _exact_count(len(labile_keys), float(kin["buried_fraction"]), rounding)
            chosen = set(rng.choice(len(labile_keys), size=n_buried, replace=False).tolist())
            pop_map = {
                key: (POP_BURIED if i in chosen else POP_ACCESSIBLE)
                for i, key in enumerate(labile_keys)
            }

        # Pass 1: pool the H atoms to convert.
        nonlabile_to_convert: List[Tuple[int, int, gemmi.Atom]] = []
        pools: Dict[str, List[Tuple[int, int, gemmi.Atom]]] = {}
        mask_labile_warnings: List[str] = []
        for ci, chain in enumerate(model):
            for ri, residue in enumerate(chain):
                resname = residue.name.strip().upper()
                selected = spec.selection.selects_residue(resname)
                for atom in residue:
                    if not atom.is_hydrogen() or atom.element.name != "H":
                        continue
                    key = (ci, ri, atom.name.strip())
                    labile = self.labile_map.get(key, False)
                    mask_hit = spec.selection.selects_atom(
                        chain.name, residue.seqid.num, atom.name
                    )
                    if labile:
                        if mask_hit:
                            mask_labile_warnings.append(
                                f"chain {chain.name} {residue.name}{residue.seqid.num} "
                                f"{atom.name.strip()}: atom mask targets a labile H; "
                                "labile atoms are handled by the exchange, mask ignored"
                            )
                        population = (
                            pop_map.get(key, POP_BURIED)
                            if mode in ("kinetic", "dual")
                            else "labile"
                        )
                        pools.setdefault(population, []).append((ci, ri, atom))
                    elif selected or mask_hit:
                        nonlabile_to_convert.append((ci, ri, atom))

        # Pass 2: exact-count selection per population.
        p_of: Dict[str, float] = {"labile": p_all} if mode == "simple" else {
            POP_ACCESSIBLE: p_acc,
            POP_BURIED: p_bur,
        }
        fallback_p = p_all if mode == "simple" else p_bur
        converted: Set[Tuple[int, int, str]] = {
            (ci, ri, atom.name.strip()) for ci, ri, atom in nonlabile_to_convert
        }
        for population, atoms in pools.items():
            p = p_of.get(population, fallback_p)
            if sampling == "exact":
                k = _exact_count(len(atoms), p, rounding)
                indices = set(rng.choice(len(atoms), size=k, replace=False).tolist())
                for i, (ci, ri, atom) in enumerate(atoms):
                    if i in indices:
                        converted.add((ci, ri, atom.name.strip()))
            else:
                for ci, ri, atom in atoms:
                    if rng.random() < p:
                        converted.add((ci, ri, atom.name.strip()))

        # Pass 3: clone, convert, trace and count.
        structure = self.structure.clone()
        model2 = structure[0]
        trace: List[Dict[str, Any]] = []
        converted_total = converted_labile = converted_nonlabile = 0
        for ci, chain in enumerate(model2):
            for ri, residue in enumerate(chain):
                for atom in residue:
                    if not atom.is_hydrogen():
                        continue
                    key = (ci, ri, atom.name.strip())
                    labile = self.labile_map.get(key, False)
                    element_before = atom.element.name
                    will_convert = key in converted and element_before == "H"
                    population = ""
                    population_atom = ""
                    probability = 0.0
                    reason = ""
                    if labile:
                        population = (
                            pop_map.get(key, POP_BURIED)
                            if mode in ("kinetic", "dual")
                            else "labile"
                        )
                        probability = p_of.get(population, fallback_p)
                        population_atom = self.atom_population_map.get(key, population)
                    elif element_before == "D":
                        reason = "already_deuterated"
                    if will_convert:
                        _convert_h_to_d(atom)
                        if labile:
                            converted_labile += 1
                        else:
                            converted_nonlabile += 1
                        converted_total += 1
                    if spec.collect_trace:
                        trace.append(
                            {
                                "chain": chain.name,
                                "residue_number": residue.seqid.num,
                                "residue_name": residue.name.strip(),
                                "atom_name": atom.name.strip(),
                                "element_after": atom.element.name,
                                "labile": int(labile),
                                "population": population,
                                "population_atom": population_atom,
                                "probability": probability,
                                "converted": int(will_convert),
                                "reason": reason,
                            }
                        )

        init = self._initial
        stats: Dict[str, Any] = {
            "total_atoms": init["total_atoms"],
            "hydrogen_atoms": init["h_total"] - converted_total,
            "deuterium_atoms": init["d_total"] + converted_total,
            "labile_H": init["h_labile"] - converted_labile,
            "labile_D": init["d_labile"] + converted_labile,
            "non_labile_H": init["h_nonlabile"] - converted_nonlabile,
            "non_labile_D": init["d_nonlabile"] + converted_nonlabile,
        }
        denominator = init["h_total"] + init["d_total"]
        stats["deuterium_pct"] = (
            100.0 * stats["deuterium_atoms"] / denominator if denominator else 0.0
        )
        nonlabile_denominator = init["h_nonlabile"] + init["d_nonlabile"]
        stats["non_labile_deut_pct"] = (
            100.0 * stats["non_labile_D"] / nonlabile_denominator
            if nonlabile_denominator
            else 0.0
        )
        labile_denominator = init["h_labile"] + init["d_labile"]
        stats["labile_deut_pct"] = (
            100.0 * stats["labile_D"] / labile_denominator if labile_denominator else 0.0
        )
        stats["mode"] = mode
        if mode in ("kinetic", "dual"):
            stats["accessible_probability"] = p_acc
            stats["buried_probability"] = p_bur
            if pop_map:
                applied = {
                    key: pop_map.get(key, POP_BURIED)
                    for key, labile in self.labile_map.items()
                    if labile
                }
                stats["accessibility_granularity"] = self.granularity
                stats["labile_H_accessible"] = sum(
                    1 for value in applied.values() if value == POP_ACCESSIBLE
                )
                stats["labile_H_buried"] = sum(
                    1 for value in applied.values() if value == POP_BURIED
                )
                by_residue = residue_populations(applied)
                stats["residues_accessible"] = sum(
                    1 for value in by_residue.values() if value == POP_ACCESSIBLE
                )
                stats["residues_buried"] = sum(
                    1 for value in by_residue.values() if value == POP_BURIED
                )
                # Mixing is a property of the atom-level classification, so it
                # is reported against the pre-aggregation map.
                reference = self.atom_population_map or applied
                per_residue: Dict[Tuple[int, int], Set[str]] = {}
                for (ci, ri, _name), value in reference.items():
                    per_residue.setdefault((ci, ri), set()).add(value)
                stats["residues_mixed"] = sum(
                    1 for values in per_residue.values() if len(values) > 1
                )
        else:
            stats["labile_probability"] = p_all
        stats["deuterium_pct_protein"] = stats["deuterium_pct"]

        run_warnings = list(self.warnings) + mask_labile_warnings
        return DeuterationResult(
            structure=structure,
            stats=stats,
            atom_trace=trace,
            warnings=run_warnings,
            policy_hash=policy.hash,
            template_hash=self.hash,
            seed=seed,
            dropped=self.dropped,
            mode=mode,
        )


@dataclass
class DeuterationResult:
    structure: gemmi.Structure
    stats: Dict[str, Any]
    atom_trace: List[Dict[str, Any]]
    warnings: List[str]
    policy_hash: str
    template_hash: str
    seed: int
    dropped: Dict[str, Dict[str, int]]
    mode: str

    def write_pdb(self, path: str) -> None:
        self.structure.write_pdb(path)

    def write_report(self, path_prefix: str) -> None:
        write_atom_trace(path_prefix + ".atom_trace.tsv", self.atom_trace)
        write_aggregates(path_prefix + ".aggregates.tsv", self.atom_trace)
        write_provenance(
            path_prefix + ".provenance.json",
            {
                "policy_hash": self.policy_hash,
                "template_hash": self.template_hash,
                "seed": self.seed,
                "mode": self.mode,
                "dropped": self.dropped,
                "warnings": self.warnings,
                "stats": self.stats,
            },
        )
