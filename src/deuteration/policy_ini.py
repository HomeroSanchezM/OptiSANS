"""Build a DeuterPDB ``Policy`` from OptiSANS INI files.

OptiSANS keeps its configuration in INI files: ``config.ini`` for the genetic
algorithm (with an optional ``[DEUTERPDB.*]`` block) and
``config_deuterate.ini`` for ``optisans deuterate``.  No policy YAML is
required: the INI values are merged over OptiSANS defaults and handed to
``Policy.from_dict()``.

Escape hatch: ``[DEUTERPDB] policy_file = <path>`` loads an explicit YAML
instead (all the other ``[DEUTERPDB.*]`` keys are then ignored).
"""

from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from src.deuterpdb import Policy

# ---------------------------------------------------------------------------
# OptiSANS defaults for the DeuterPDB engine.
#
# These are applied BEFORE the INI values (so an absent config file still
# gives a working, reproducible setup) and differ from the standalone
# DeuterPDB defaults on purpose:
#   * protonation = auto       (add hydrogens when the input is not protonated)
#   * seed = 42                (same PDBs from one run to the next)
# ---------------------------------------------------------------------------
OPTISANS_DEFAULTS: Dict[str, Any] = {
    "template": {
        "protonation": "auto",
        "force_field": "amber99sb-ildn",
        "strip_water_ligand": True,
    },
    "lability": {},
    "accessibility": {
        "enabled": False,
    },
    "labile_exchange": {
        "mode": "simple",
        "sampling": "exact",
        "rounding": "floor",
        "seed": 42,
        "kinetic": {},
    },
    "output": {
        "write_pdb": True,
        "write_report": False,
    },
}


class DeuterPDBConfigError(ValueError):
    """Raised when a [DEUTERPDB.*] section holds an unknown key or bad value."""


# section name -> (policy section, {ini key: type})
_SECTION_MAP: Dict[str, Tuple[str, Dict[str, type]]] = {
    "DEUTERPDB.TEMPLATE": (
        "template",
        {
            "protonation": str,
            "gmx_path": str,
            "force_field": str,
            "strip_water_ligand": bool,
            "model": int,
        },
    ),
    "DEUTERPDB.LABILITY": (
        "lability",
        {
            "cutoff_angstrom": float,
            "ambiguity_tolerance": float,
            "quarantine": str,
        },
    ),
    "DEUTERPDB.ACCESSIBILITY": (
        "accessibility",
        {
            "enabled": bool,
            "granularity": str,
            "sasa_probe": float,
            "sasa_n_points": int,
            "threshold_angstrom2": float,
            "threshold_mode": str,
            "use_hbond_protection": bool,
            "hbond_max_h_o": float,
            "hbond_min_angle": float,
        },
    ),
    "DEUTERPDB.LABILE_EXCHANGE": (
        "labile_exchange",
        {
            "mode": str,
            "sampling": str,
            "rounding": str,
            "seed": int,
            "probability_burried": float,
        },
    ),
    "DEUTERPDB.KINETIC": (
        "labile_exchange.kinetic",
        {
            "accessible_rate": float,
            "buried_rate": float,
            "time_minutes": float,
            "buried_fraction_source": str,
            "buried_fraction": float,
        },
    ),
}


def _read_value(parser: ConfigParser, section: str, key: str, kind: type) -> Any:
    """Read one INI value with a type; empty string means 'not set' (None)."""
    try:
        if kind is bool:
            return parser.getboolean(section, key)
        if kind is int:
            raw = parser.get(section, key)
            return None if raw.strip() == "" else int(raw)
        if kind is float:
            raw = parser.get(section, key)
            return None if raw.strip() == "" else float(raw)
        raw = parser.get(section, key)
        return None if raw.strip() == "" else raw.strip()
    except ValueError as exc:  # bad literal for the requested type
        raise DeuterPDBConfigError(
            f"[{section}] {key}: invalid value ({exc})"
        ) from exc


def _sections_to_overrides(parser: ConfigParser) -> Dict[str, Any]:
    """Convert every known [DEUTERPDB.*] section into a policy override dict."""
    overrides: Dict[str, Any] = {}
    for section in parser.sections():
        if not section.upper().startswith("DEUTERPDB"):
            continue
        if section == "DEUTERPDB":
            continue  # handled by the caller (policy_file / write_report)
        spec = _SECTION_MAP.get(section.upper())
        if spec is None:
            raise DeuterPDBConfigError(
                f"unknown section [{section}] "
                f"(expected one of: {', '.join(sorted(_SECTION_MAP))})"
            )
        policy_path, keys = spec
        unknown = set(parser.options(section)) - set(keys)
        if unknown:
            raise DeuterPDBConfigError(
                f"[{section}] unknown key(s): {', '.join(sorted(unknown))}"
            )
        values: Dict[str, Any] = {}
        for key, kind in keys.items():
            if parser.has_option(section, key):
                value = _read_value(parser, section, key, kind)
                if value is not None:
                    values[key] = value
        if values:
            node = overrides
            for part in policy_path.split("."):
                node = node.setdefault(part, {})
            node.update(values)
    return overrides


def _deep_merge(base: Dict[str, Any], extra: Dict[str, Any]) -> Dict[str, Any]:
    out = {key: (dict(value) if isinstance(value, dict) else value) for key, value in base.items()}
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_policy_and_options(
    ini_path: Optional[str] = None,
    overrides: Optional[Dict[str, Any]] = None,
    default_write_report: bool = False,
) -> Tuple[Policy, Dict[str, Any]]:
    """Build the policy and the extra OptiSANS-side options.

    Args:
        ini_path:  optional INI file (``config.ini`` or ``config_deuterate.ini``).
        overrides: policy overrides applied last (CLI wins over the INI).
        default_write_report: value used when neither the INI nor *overrides*
                              sets ``write_report`` (False for ``run``, True for
                              ``optisans deuterate``).

    Returns:
        ``(policy, options)`` where *options* carries ``write_report`` (bool)
        and ``write_pdb`` (bool) as read from the file.

    Raises:
        FileNotFoundError: if *ini_path* does not exist.
        DeuterPDBConfigError: for unknown sections/keys or bad values.
    """
    # Inline comments ("key = value   # comment") are allowed in the INI files.
    parser = ConfigParser(inline_comment_prefixes=("#", ";"))
    if ini_path is not None:
        if not Path(ini_path).exists():
            raise FileNotFoundError(f"Config file not found: {ini_path}")
        parser.read(ini_path)

    options: Dict[str, Any] = {
        "write_report": default_write_report,
        "write_pdb": OPTISANS_DEFAULTS["output"]["write_pdb"],
    }

    # --- escape hatch: an explicit policy YAML ---
    if parser.has_section("DEUTERPDB") and parser.has_option("DEUTERPDB", "policy_file"):
        raw = parser.get("DEUTERPDB", "policy_file").strip()
        if raw:
            policy = Policy.from_yaml(raw)
            return policy, options

    # --- [DEUTERPDB] scalars handled here (not part of the policy schema) ---
    if parser.has_section("DEUTERPDB"):
        allowed = {"policy_file", "write_report", "write_pdb"}
        unknown = set(parser.options("DEUTERPDB")) - allowed
        if unknown:
            raise DeuterPDBConfigError(
                f"[DEUTERPDB] unknown key(s): {', '.join(sorted(unknown))} "
                f"(expected: {', '.join(sorted(allowed))})"
            )
        if parser.has_option("DEUTERPDB", "write_report"):
            options["write_report"] = parser.getboolean("DEUTERPDB", "write_report")
        if parser.has_option("DEUTERPDB", "write_pdb"):
            options["write_pdb"] = parser.getboolean("DEUTERPDB", "write_pdb")

    sections = _sections_to_overrides(parser)
    merged = _deep_merge(OPTISANS_DEFAULTS, sections)
    merged["output"] = {
        "write_pdb": options["write_pdb"],
        "write_report": options["write_report"],
    }
    if overrides:
        merged = _deep_merge(merged, overrides)

    try:
        policy = Policy.from_dict(merged)
    except (ValueError, TypeError) as exc:
        raise DeuterPDBConfigError(f"invalid DEUTERPDB configuration: {exc}") from exc
    return policy, options


def load_policy(
    ini_path: Optional[str] = None,
    overrides: Optional[Dict[str, Any]] = None,
) -> Policy:
    """Convenience wrapper returning only the policy."""
    policy, _ = load_policy_and_options(ini_path, overrides)
    return policy


def read_deuteration_ini(path: Optional[str]) -> Dict[str, Any]:
    """Read the legacy-style ``[DEUTERATION]`` and ``[AMINO_ACIDS]`` sections.

    Returns a dict with any of: ``input_pdb``, ``output_pdb``, ``d2o_percent``,
    ``amino_acids`` (list of 3-letter codes set to true).  Missing file or
    sections give an empty dict; unknown keys in those historical sections are
    tolerated.
    """
    out: Dict[str, Any] = {}
    if not path:
        return out
    parser = ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(path)
    if parser.has_section("DEUTERATION"):
        for key in ("input_pdb", "output_pdb"):
            if parser.has_option("DEUTERATION", key):
                value = parser.get("DEUTERATION", key).strip()
                if value:
                    out[key] = value
        if parser.has_option("DEUTERATION", "d2o_percent"):
            out["d2o_percent"] = parser.getfloat("DEUTERATION", "d2o_percent")
    if parser.has_section("AMINO_ACIDS"):
        out["amino_acids"] = [
            code.upper()
            for code in parser.options("AMINO_ACIDS")
            if parser.getboolean("AMINO_ACIDS", code)
        ]
    return out
