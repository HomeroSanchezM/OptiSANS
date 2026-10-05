"""Report writing: per-atom trace, aggregates and provenance."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from typing import Any, Dict, List

TRACE_COLUMNS = [
    "chain",
    "residue_number",
    "residue_name",
    "atom_name",
    "element_after",
    "labile",
    "population",
    "probability",
    "converted",
    "reason",
]


def write_atom_trace(path: str, trace: List[Dict[str, Any]]) -> None:
    """Write the per-atom decision trace. Every row carries the inputs of the
    decision (population, probability), so figures are reproducible from the
    table alone."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TRACE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in trace:
            writer.writerow({key: row.get(key, "") for key in TRACE_COLUMNS})


def write_aggregates(path: str, trace: List[Dict[str, Any]]) -> None:
    """Aggregated counts per residue and per amino-acid type."""
    by_residue: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"n": 0, "n_labile": 0, "n_d": 0}
    )
    by_aatype: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"n": 0, "n_labile": 0, "n_d": 0}
    )
    for row in trace:
        res_key = f"{row['chain']}|{row['residue_number']}|{row['residue_name']}"
        for table, key in ((by_residue, res_key), (by_aatype, row["residue_name"])):
            table[key]["n"] += 1
            if row["labile"]:
                table[key]["n_labile"] += 1
            if row["element_after"] == "D":
                table[key]["n_d"] += 1
    with open(path, "w", newline="") as fh:
        fh.write("level\tkey\tn\tn_labile\tn_D\tpct_D\n")
        for level, table in (("residue", by_residue), ("aa_type", by_aatype)):
            for key in sorted(table):
                entry = table[key]
                pct = 100.0 * entry["n_d"] / entry["n"] if entry["n"] else 0.0
                fh.write(
                    f"{level}\t{key}\t{entry['n']}\t{entry['n_labile']}\t"
                    f"{entry['n_d']}\t{pct:.2f}\n"
                )


def write_provenance(path: str, meta: Dict[str, Any]) -> None:
    """Write policy hash, template hash, seed, dropped residues and stats."""
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
