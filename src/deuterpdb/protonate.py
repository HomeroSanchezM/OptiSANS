"""GROMACS protonation with pdbfixer fallback. Water and ligands are stripped first."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import gemmi

from .policy import Policy


class ProtonationError(RuntimeError):
    """Raised when gromacs or pdbfixer cannot produce a protonated structure."""


def strip_water_and_ligand(
    structure: gemmi.Structure,
) -> Tuple[gemmi.Structure, Dict[str, Dict[str, int]]]:
    """Remove HETATM residues (waters and ligands).

    Returns (structure, dropped) where dropped maps residue name to
    {"count": number of residues, "h": number of hydrogens removed}.
    """
    st = structure.clone()
    dropped: Dict[str, Dict[str, int]] = {}
    for model in st:
        for chain in model:
            for i in range(len(chain) - 1, -1, -1):
                residue = chain[i]
                if residue.het_flag == "H":
                    name = residue.name.strip() or "HET"
                    entry = dropped.setdefault(name, {"count": 0, "h": 0})
                    entry["count"] += 1
                    entry["h"] += sum(1 for atom in residue if atom.is_hydrogen())
                    del chain[i]
    return st, dropped


def _rename_his_to_hid(structure: gemmi.Structure) -> None:
    """pdb2gmx can prompt interactively for HIS protonation; renaming HIS to
    HID (ND1 protonated) makes the run deterministic."""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name.strip().upper() == "HIS":
                    residue.name = "HID"


def run_gromacs(
    structure: gemmi.Structure,
    gmx_path: str | None = None,
    force_field: str = "amber99sb-ildn",
) -> gemmi.Structure:
    """Protonate with gromacs pdb2gmx. The gmx path must be explicit when the
    environment wrapper is broken (known for the OptiSANS pixi env)."""
    gmx = gmx_path or shutil.which("gmx")
    if not gmx:
        raise ProtonationError("gmx executable not found; set template.gmx_path in the policy")
    with tempfile.TemporaryDirectory(prefix="deuterpdb_") as td:
        td = Path(td)
        in_pdb = td / "in.pdb"
        out_pdb = td / "out.pdb"
        structure.write_pdb(str(in_pdb))
        cmd = [
            gmx,
            "pdb2gmx",
            "-f",
            str(in_pdb),
            "-o",
            str(out_pdb),
            "-ff",
            force_field,
            "-water",
            "none",
            "-ignh",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if proc.returncode != 0:
            raise ProtonationError(
                f"pdb2gmx failed (rc={proc.returncode}): {proc.stderr[-2000:]}"
            )
        return gemmi.read_structure(str(out_pdb))


def run_pdbfixer(structure: gemmi.Structure) -> gemmi.Structure:
    """Fallback protonation with pdbfixer + openmm (adds hydrogens at pH 7)."""
    try:
        import pdbfixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise ProtonationError("pdbfixer/openmm not installed") from exc
    with tempfile.TemporaryDirectory(prefix="deuterpdb_") as td:
        td = Path(td)
        in_pdb = td / "in.pdb"
        out_pdb = td / "out.pdb"
        structure.write_pdb(str(in_pdb))
        fixer = pdbfixer.PDBFixer(filename=str(in_pdb))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        with open(out_pdb, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
        return gemmi.read_structure(str(out_pdb))


def protonate_structure(structure: gemmi.Structure, policy: Policy) -> gemmi.Structure:
    """Strip water/ligands, then protonate with gromacs (pdbfixer fallback)."""
    structure, _ = strip_water_and_ligand(structure)
    _rename_his_to_hid(structure)
    try:
        return run_gromacs(
            structure,
            policy.get("template", "gmx_path"),
            policy.get("template", "force_field"),
        )
    except ProtonationError as gmx_error:
        try:
            return run_pdbfixer(structure)
        except ProtonationError:
            raise gmx_error


_RESOURCE_DIR = Path(__file__).resolve().parent.parent.parent / "resources"


def load_residue_h_counts(path: str | None = None) -> Dict[str, Dict[str, int]]:
    """Expected total/labile/non-labile H per residue type (internal residues)."""
    p = Path(path) if path else _RESOURCE_DIR / "residue_h_counts.tsv"
    out: Dict[str, Dict[str, int]] = {}
    for line in p.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] == "residue":
            continue
        out[parts[0]] = {
            "h_total": int(parts[1]),
            "h_labile": int(parts[2]),
            "h_nonlabile": int(parts[3]),
        }
    return out


def validate_protonation(
    structure: gemmi.Structure,
    tolerance: int = 1,
    expected: Dict[str, Dict[str, int]] | None = None,
) -> List[str]:
    """Check that no HETATM remains and that per-residue H counts match the
    reference table within tolerance. Returns a list of issues (empty is OK)."""
    expected = expected if expected is not None else load_residue_h_counts()
    issues: List[str] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.het_flag == "H":
                    issues.append(
                        f"HETATM residue {residue.name.strip()} still present after stripping"
                    )
                    continue
                name = residue.name.strip().upper()
                if name in ("HID", "HIE", "HIP"):
                    name = "HIS"
                if name not in expected:
                    continue
                n_h = sum(1 for atom in residue if atom.element.name == "H")
                diff = n_h - expected[name]["h_total"]
                if abs(diff) > tolerance:
                    issues.append(
                        f"chain {chain.name} {name}{residue.seqid.num}: {n_h} H, "
                        f"expected {expected[name]['h_total']} (+-{tolerance})"
                    )
    return issues
