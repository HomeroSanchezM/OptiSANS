"""Immutable policy for a DeuterPDB run, loaded and validated from YAML.

All thresholds and modes are declared in a policy file. Any unknown or
mis-typed key is an error (DeuteRater-style YAML settings silently ignore
malformed keys; this package does not). The policy hash is embedded in every
output for full traceability.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

import yaml

_DEFAULTS: Dict[str, Any] = {
    "template": {
        "protonation": "never",
        "gmx_path": None,
        "force_field": "amber99sb-ildn",
        "strip_water_ligand": True,
        "model": None,
    },
    "lability": {
        "cutoff_angstrom": 1.4,
        "ambiguity_tolerance": 0.25,
        "quarantine": "warn",
    },
    "accessibility": {
        "enabled": False,
        "granularity": "residue",
        "sasa_probe": 1.4,
        "sasa_n_points": 92,
        "threshold_angstrom2": 5.0,
        "threshold_mode": "absolute",
        "use_hbond_protection": False,
        "hbond_max_h_o": 2.5,
        "hbond_min_angle": 120.0,
    },
    "labile_exchange": {
        "mode": "simple",
        "d2o_percent": 0.0,
        "probability_burried": None,
        "kinetic": {
            "accessible_rate": 1.0,
            "buried_rate": 0.01,
            "time_minutes": 60.0,
            "buried_fraction_source": "sasa",
            "buried_fraction": None,
        },
        "sampling": "exact",
        "rounding": "floor",
        "seed": None,
    },
    "selection": {"residue_names": [], "all": False},
    "output": {"write_pdb": True, "write_report": True},
}

_SECTIONS = set(_DEFAULTS)

_TYPES: Dict[str, Any] = {
    "template.protonation": str,
    "template.gmx_path": (str, type(None)),
    "template.force_field": str,
    "template.strip_water_ligand": bool,
    "template.model": (int, type(None)),
    "lability.cutoff_angstrom": float,
    "lability.ambiguity_tolerance": float,
    "lability.quarantine": str,
    "accessibility.enabled": bool,
    "accessibility.granularity": str,
    "accessibility.sasa_probe": float,
    "accessibility.sasa_n_points": int,
    "accessibility.threshold_angstrom2": float,
    "accessibility.threshold_mode": str,
    "accessibility.use_hbond_protection": bool,
    "accessibility.hbond_max_h_o": float,
    "accessibility.hbond_min_angle": float,
    "labile_exchange.mode": str,
    "labile_exchange.d2o_percent": float,
    "labile_exchange.probability_burried": (float, type(None)),
    "labile_exchange.kinetic.accessible_rate": float,
    "labile_exchange.kinetic.buried_rate": float,
    "labile_exchange.kinetic.time_minutes": float,
    "labile_exchange.kinetic.buried_fraction_source": str,
    "labile_exchange.kinetic.buried_fraction": (float, type(None)),
    "labile_exchange.sampling": str,
    "labile_exchange.rounding": str,
    "labile_exchange.seed": (int, type(None)),
    "selection.residue_names": list,
    "selection.all": bool,
    "output.write_pdb": bool,
    "output.write_report": bool,
}

_ENUMS: Dict[str, set] = {
    "template.protonation": {"never", "always", "auto"},
    "lability.quarantine": {"warn", "error"},
    "accessibility.threshold_mode": {"absolute", "relative"},
    "accessibility.granularity": {"atom", "residue"},
    "labile_exchange.mode": {"simple", "kinetic", "dual"},
    "labile_exchange.sampling": {"exact", "bernoulli"},
    "labile_exchange.rounding": {"floor", "round", "ceil"},
    "labile_exchange.kinetic.buried_fraction_source": {"sasa", "fixed"},
}


def _merge(base: Dict[str, Any], override: Dict[str, Any], section: str) -> Dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if key not in merged:
            raise ValueError(f"policy: unknown key '{section}.{key}'")
        if isinstance(value, dict):
            merged[key] = _merge(merged[key], value, f"{section}.{key}")
        else:
            merged[key] = value
    return merged


def _lookup(data: Dict[str, Any], path: str) -> Any:
    node = data
    for part in path.split("."):
        node = node[part]
    return node


def _type_ok(value: Any, expected: Any) -> bool:
    if expected is float:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected is int:
        return isinstance(value, int) and not isinstance(value, bool)
    return isinstance(value, expected)


def _validate(data: Dict[str, Any]) -> None:
    for path, expected in _TYPES.items():
        value = _lookup(data, path)
        if not _type_ok(value, expected):
            raise ValueError(
                f"policy: key '{path}' must be {expected}, got {type(value).__name__}"
            )
    for path, values in _ENUMS.items():
        value = _lookup(data, path)
        if value not in values:
            raise ValueError(f"policy: key '{path}' must be one of {sorted(values)}, got {value!r}")

    lx = data["labile_exchange"]
    if not 0.0 <= lx["d2o_percent"] <= 100.0:
        raise ValueError("policy: labile_exchange.d2o_percent must be in [0, 100]")

    probability_burried = lx["probability_burried"]
    if lx["mode"] == "dual":
        if probability_burried is None:
            raise ValueError(
                "policy: labile_exchange.probability_burried is required when mode is dual"
            )
        if not 0.0 <= probability_burried <= 1.0:
            raise ValueError(
                "policy: labile_exchange.probability_burried must be in [0, 1]"
            )
        if not data["accessibility"]["enabled"]:
            raise ValueError("policy: mode dual requires accessibility.enabled: true")
    elif probability_burried is not None:
        raise ValueError(
            "policy: labile_exchange.probability_burried is only valid with mode dual"
        )

    kin = lx["kinetic"]
    if kin["accessible_rate"] <= 0 or kin["buried_rate"] <= 0:
        raise ValueError("policy: kinetic rates must be strictly positive")
    if kin["time_minutes"] < 0:
        raise ValueError("policy: kinetic.time_minutes must be >= 0")
    if kin["buried_fraction_source"] == "fixed":
        frac = kin["buried_fraction"]
        if frac is None or not 0.0 <= frac <= 1.0:
            raise ValueError(
                "policy: kinetic.buried_fraction must be in [0, 1] when source is fixed"
            )
    elif kin["buried_fraction"] is not None:
        raise ValueError(
            "policy: kinetic.buried_fraction is only valid with buried_fraction_source: fixed"
        )

    if lx["mode"] == "kinetic":
        acc = data["accessibility"]
        if not acc["enabled"] and kin["buried_fraction_source"] != "fixed":
            raise ValueError(
                "policy: mode kinetic requires accessibility.enabled: true or "
                "labile_exchange.kinetic.buried_fraction_source: fixed"
            )

    if data["accessibility"]["enabled"]:
        if data["accessibility"]["sasa_n_points"] < 20:
            raise ValueError("policy: accessibility.sasa_n_points must be >= 20")
        if data["accessibility"]["threshold_angstrom2"] <= 0:
            raise ValueError("policy: accessibility.threshold_angstrom2 must be > 0")


def _canonical(data: Dict[str, Any]) -> Dict[str, Any]:
    def clean(node: Any) -> Any:
        if isinstance(node, dict):
            return {key: clean(value) for key, value in sorted(node.items())}
        if isinstance(node, float):
            return round(node, 8)
        return node

    return clean(data)


@dataclass(frozen=True)
class Policy:
    """Validated, immutable run configuration. `hash` identifies the exact content."""

    data: Dict[str, Any] = field(default_factory=dict, compare=False)
    hash: str = ""

    @classmethod
    def _from_raw(cls, raw: Dict[str, Any]) -> "Policy":
        """Merge *raw* over the defaults, validate and hash (shared by all factories)."""
        merged: Dict[str, Any] = {}
        for section, defaults in _DEFAULTS.items():
            merged[section] = _merge(defaults, raw.get(section, {}), section)
        for section in raw:
            if section not in _SECTIONS:
                raise ValueError(f"policy: unknown section '{section}'")
        _validate(merged)
        canonical = yaml.safe_dump(_canonical(merged), sort_keys=True)
        return cls(data=merged, hash=hashlib.sha256(canonical.encode()).hexdigest()[:16])

    @classmethod
    def from_yaml(cls, path: str) -> "Policy":
        raw = yaml.safe_load(Path(path).read_text())
        if raw is None:
            raw = {}
        if not isinstance(raw, dict):
            raise ValueError(f"policy file {path} must contain a YAML mapping")
        return cls._from_raw(raw)

    @classmethod
    def from_dict(cls, overrides: Dict[str, Any] | None = None) -> "Policy":
        """Build a policy from an in-memory mapping of partial overrides.

        Used by host tools (OptiSANS) that keep their configuration in an INI
        file: no policy YAML is needed on disk.  Unknown sections or keys are
        errors, exactly as in ``from_yaml``.
        """
        raw = overrides or {}
        if not isinstance(raw, dict):
            raise ValueError("policy overrides must be a mapping")
        return cls._from_raw(raw)

    @classmethod
    def defaults(cls) -> "Policy":
        """Built-in defaults (labile_exchange.mode = simple, no accessibility)."""
        return cls._from_raw({})

    def get(self, *path: str) -> Any:
        node = self.data
        for key in path:
            node = node[key]
        return node
