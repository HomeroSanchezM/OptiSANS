"""PDB loading, model selection, altloc handling and writing with gemmi."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

import gemmi

from .policy import Policy


class StructureError(ValueError):
    """Raised when the input PDB cannot be used as a deuteration template."""


def load_structure(path: str, policy: Policy) -> gemmi.Structure:
    """Read a PDB, keep the requested model (error if several and none chosen),
    and collapse alternate conformations to the highest occupancy."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")
    st = gemmi.read_structure(str(p))
    if len(st) == 0:
        raise StructureError(f"{path}: no model in structure")
    if len(st) > 1:
        model_idx: Optional[int] = policy.get("template", "model")
        if model_idx is None:
            raise StructureError(
                f"{path}: {len(st)} models present; set template.model in the policy"
            )
        if not 0 <= model_idx < len(st):
            raise StructureError(f"{path}: template.model {model_idx} out of range")
        st = _single_model(st, model_idx)
    for model in st:
        model.remove_alternative_conformations()
    return st


def _single_model(st: gemmi.Structure, idx: int) -> gemmi.Structure:
    out = st.clone()
    del out[1:]
    return out


def write_pdb(structure: gemmi.Structure, path: str) -> None:
    """Write a structure to a PDB file."""
    structure.write_pdb(path)


def structure_hash(structure: gemmi.Structure) -> str:
    """Content hash of a single-model structure: chain, residue, atom name,
    element and coordinates."""
    h = hashlib.sha256()
    model = structure[0]
    for chain in model:
        for residue in chain:
            for atom in residue:
                line = (
                    f"{chain.name}|{residue.name}|{residue.seqid.num}|"
                    f"{atom.name.strip()}|{atom.element.name}|"
                    f"{atom.pos.x:.3f}|{atom.pos.y:.3f}|{atom.pos.z:.3f}"
                )
                h.update(line.encode())
    return h.hexdigest()[:16]
