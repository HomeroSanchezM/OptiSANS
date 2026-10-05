"""Command line interface: run, inspect, protonate, batch."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .deuterate import DeuterationSpec, Template
from .io import load_structure
from .policy import Policy
from .protonate import (
    ProtonationError,
    _rename_his_to_hid,
    run_gromacs,
    run_pdbfixer,
    strip_water_and_ligand,
    validate_protonation,
)
from .select import Selection


def _parse_selection(text: str) -> Selection:
    names = [name.strip().upper() for name in (text or "").split(",") if name.strip()]
    return Selection.by_residue_names(names)


def cmd_run(args: argparse.Namespace) -> None:
    policy = Policy.from_yaml(args.policy)
    template = Template.from_pdb(args.input, policy)
    spec = DeuterationSpec(
        selection=_parse_selection(args.selection),
        mode=args.mode,
        d2o_percent=args.d2o,
        seed=args.seed,
    )
    result = template.deuterate(spec)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    if policy.get("output", "write_pdb"):
        result.write_pdb(str(out / "deuterated.pdb"))
    if policy.get("output", "write_report"):
        result.write_report(str(out / "report"))
    print(json.dumps(result.stats, indent=2))


def cmd_inspect(args: argparse.Namespace) -> None:
    policy = Policy.from_yaml(args.policy)
    template = Template.from_pdb(args.input, policy)
    print(f"template_hash: {template.hash}")
    print(f"dropped: {json.dumps(template.dropped)}")
    print(f"n_labile: {template.n_labile}  n_nonlabile: {template.n_nonlabile}")
    for warning in template.warnings:
        print(f"warning: {warning}")
    if template.population_map:
        from collections import Counter

        print(f"populations: {dict(Counter(template.population_map.values()))}")


def cmd_protonate(args: argparse.Namespace) -> None:
    policy = Policy.from_yaml(args.policy)
    structure = load_structure(args.input, policy)
    structure, dropped = strip_water_and_ligand(structure)
    _rename_his_to_hid(structure)
    try:
        protonated = run_gromacs(
            structure,
            policy.get("template", "gmx_path"),
            policy.get("template", "force_field"),
        )
    except ProtonationError as gmx_error:
        print(f"gromacs failed ({gmx_error}); falling back to pdbfixer", file=sys.stderr)
        protonated = run_pdbfixer(structure)
    issues = validate_protonation(protonated)
    protonated.write_pdb(args.output)
    print(json.dumps({"dropped": dropped, "validation_issues": issues}, indent=2))


def cmd_batch(args: argparse.Namespace) -> None:
    policy = Policy.from_yaml(args.policy)
    template = Template.from_pdb(args.input, policy)
    specs = json.loads(Path(args.specs).read_text())
    if not isinstance(specs, list):
        raise ValueError("specs file must contain a JSON list of spec objects")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for i, item in enumerate(specs):
        spec = DeuterationSpec(
            selection=Selection.by_residue_names(item.get("selection", [])),
            mode=item.get("mode"),
            d2o_percent=item.get("d2o_percent"),
            seed=item.get("seed"),
        )
        result = template.deuterate(spec)
        prefix = f"spec_{i:04d}"
        if policy.get("output", "write_pdb"):
            result.write_pdb(str(out / f"{prefix}.pdb"))
        if policy.get("output", "write_report"):
            result.write_report(str(out / prefix))
        results.append({"index": i, "stats": result.stats})
    print(json.dumps(results, indent=2))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="deuterpdb",
        description="PDB parsing, protonation, lability and simple, kinetic or dual deuteration",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="deuterate one structure")
    p_run.add_argument("input", help="input PDB")
    p_run.add_argument("--policy", required=True, help="policy YAML")
    p_run.add_argument("-o", "--output", default="out", help="output directory")
    p_run.add_argument("--selection", default="", help="comma-separated residue names")
    p_run.add_argument("--mode", choices=["simple", "kinetic", "dual"], default=None)
    p_run.add_argument("--d2o", type=float, default=None)
    p_run.add_argument("--seed", type=int, default=None)
    p_run.set_defaults(func=cmd_run)

    p_inspect = sub.add_parser("inspect", help="report lability and populations without writing")
    p_inspect.add_argument("input")
    p_inspect.add_argument("--policy", required=True)
    p_inspect.set_defaults(func=cmd_inspect)

    p_protonate = sub.add_parser(
        "protonate", help="strip water/ligands and protonate with gromacs"
    )
    p_protonate.add_argument("input")
    p_protonate.add_argument("--policy", required=True)
    p_protonate.add_argument("-o", "--output", required=True)
    p_protonate.set_defaults(func=cmd_protonate)

    p_batch = sub.add_parser("batch", help="run many specs on one template")
    p_batch.add_argument("input")
    p_batch.add_argument("--policy", required=True)
    p_batch.add_argument("--specs", required=True, help="JSON list of spec objects")
    p_batch.add_argument("-o", "--output", default="out")
    p_batch.set_defaults(func=cmd_batch)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
