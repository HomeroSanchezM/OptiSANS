#!/usr/bin/env python3
"""
Amino-acid tables for OptiSANS: single source of truth.
=======================================================

Every module that needs the canonical amino-acid list, the 18 linked-group
genes (ASN+ASP and GLU+GLN share one gene each) or the chromosome/restriction
vector helpers imports them from here.

Linked pairs are a hard rule of the tool: selecting (or restricting) one
member of a pair always applies to both.  ``linked_partners`` and the
``*_linked_*`` helpers below implement that rule, so callers can reuse the
same logic and print consistent warnings.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass
class AminoAcid:
    """
    Represents an amino acid (or a linked pair of amino acids) with its codes.

    Attributes:
        name (str): Full name, e.g. "Alanine" or "Asparagine/Aspartic acid"
        code_3 (str): 3-letter code, e.g. "ALA" or "ASN+ASP" for linked pairs
        code_1 (str): 1-letter code, e.g. "A" or "ND" for linked pairs
    """

    name: str
    code_3: str
    code_1: str


# Full canonical list of 20 standard amino acids (used for PDB operations).
# Indices are referenced throughout the codebase.
AMINO_ACIDS = [
    AminoAcid("Alanine", "ALA", "A"),  # index 0
    AminoAcid("Arginine", "ARG", "R"),  # index 1
    AminoAcid("Asparagine", "ASN", "N"),  # index 2
    AminoAcid("Aspartic acid", "ASP", "D"),  # index 3
    AminoAcid("Cysteine", "CYS", "C"),  # index 4
    AminoAcid("Glutamic acid", "GLU", "E"),  # index 5
    AminoAcid("Glutamine", "GLN", "Q"),  # index 6
    AminoAcid("Glycine", "GLY", "G"),  # index 7
    AminoAcid("Histidine", "HIS", "H"),  # index 8
    AminoAcid("Isoleucine", "ILE", "I"),  # index 9
    AminoAcid("Leucine", "LEU", "L"),  # index 10
    AminoAcid("Lysine", "LYS", "K"),  # index 11
    AminoAcid("Methionine", "MET", "M"),  # index 12
    AminoAcid("Phenylalanine", "PHE", "F"),  # index 13
    AminoAcid("Proline", "PRO", "P"),  # index 14
    AminoAcid("Serine", "SER", "S"),  # index 15
    AminoAcid("Threonine", "THR", "T"),  # index 16
    AminoAcid("Tryptophan", "TRP", "W"),  # index 17
    AminoAcid("Tyrosine", "TYR", "Y"),  # index 18
    AminoAcid("Valine", "VAL", "V"),  # index 19
]

# 3-letter code -> position in AMINO_ACIDS (0-19).  Both names exist for
# backwards compatibility (GA.py used AA_DICT, pdb_deuteration.py AA_INDEX).
AA_DICT: Dict[str, int] = {aa.code_3: i for i, aa in enumerate(AMINO_ACIDS)}
AA_INDEX: Dict[str, int] = AA_DICT

# LINKED AMINO ACID GROUPS
# Each entry is a list of 3-letter codes that always deuterate together.
# Linked pairs:
#   ASN + ASP
#   GLU + GLN
LINKED_AA_GROUPS: List[List[str]] = [
    ["ALA"],  # group  0  -> AMINO_ACIDS index 0
    ["ARG"],  # group  1  -> AMINO_ACIDS index 1
    ["ASN", "ASP"],  # group  2  -> AMINO_ACIDS indices 2, 3
    ["CYS"],  # group  3  -> AMINO_ACIDS index 4
    ["GLU", "GLN"],  # group  4  -> AMINO_ACIDS indices 5, 6
    ["GLY"],  # group  5  -> AMINO_ACIDS index 7
    ["HIS"],  # group  6  -> AMINO_ACIDS index 8
    ["ILE"],  # group  7  -> AMINO_ACIDS index 9
    ["LEU"],  # group  8  -> AMINO_ACIDS index 10
    ["LYS"],  # group  9  -> AMINO_ACIDS index 11
    ["MET"],  # group 10  -> AMINO_ACIDS index 12
    ["PHE"],  # group 11  -> AMINO_ACIDS index 13
    ["PRO"],  # group 12  -> AMINO_ACIDS index 14
    ["SER"],  # group 13  -> AMINO_ACIDS index 15
    ["THR"],  # group 14  -> AMINO_ACIDS index 16
    ["TRP"],  # group 15  -> AMINO_ACIDS index 17
    ["TYR"],  # group 16  -> AMINO_ACIDS index 18
    ["VAL"],  # group 17  -> AMINO_ACIDS index 19
]

def _build_group_index() -> Dict[str, int]:
    """3-letter code -> index of its group in LINKED_AA_GROUPS."""
    index: Dict[str, int] = {}
    for group_idx, group in enumerate(LINKED_AA_GROUPS):
        for code in group:
            index[code] = group_idx
    return index


def _build_effective_amino_acids() -> List[AminoAcid]:
    """One entry per linked group (18 total); linked pairs get a synthetic entry."""
    effective: List[AminoAcid] = []
    for group in LINKED_AA_GROUPS:
        if len(group) == 1:
            # Standalone AA — use the canonical entry from AMINO_ACIDS
            effective.append(next(aa for aa in AMINO_ACIDS if aa.code_3 == group[0]))
        else:
            # Linked pair — build a synthetic AminoAcid
            members = [next(aa for aa in AMINO_ACIDS if aa.code_3 == c) for c in group]
            effective.append(
                AminoAcid(
                    name="/".join(aa.name for aa in members),
                    code_3="+".join(aa.code_3 for aa in members),
                    code_1="".join(aa.code_1 for aa in members),
                )
            )
    return effective


# Dictionary 3-letter code -> group index in LINKED_AA_GROUPS
AA_GROUP_INDEX: Dict[str, int] = _build_group_index()


# EFFECTIVE_AMINO_ACIDS
# 18-element list, one entry per linked group.
# This is the list the Chromosome and PopulationGenerator work with.
EFFECTIVE_AMINO_ACIDS: List[AminoAcid] = _build_effective_amino_acids()


# Number of effective genes in a chromosome (18)
N_EFFECTIVE_AA: int = len(EFFECTIVE_AMINO_ACIDS)  # 18

# Default: all effective genes are modifiable
restrictions: List[bool] = [True] * N_EFFECTIVE_AA


# ---------------------------------------------------------------------------
# Vector helpers
# ---------------------------------------------------------------------------

def expand_deuteration_vector(deut_18: Sequence[bool]) -> List[bool]:
    """
    Expand an 18-element chromosome deuteration vector to the full 20-element
    vector aligned with AMINO_ACIDS ordering.

    For linked groups (ASN+ASP, GLU+GLN), both amino acids in the pair receive
    the same deuteration state as their shared gene.

    Args:
        deut_18: 18-element boolean list (one entry per effective gene /
                 linked group in LINKED_AA_GROUPS).

    Returns:
        20-element boolean list aligned with AMINO_ACIDS ordering.

    Raises:
        ValueError: If the input length is neither 18 nor 20.

    Note:
        If a 20-element vector is passed it is returned unchanged.
    """
    if len(deut_18) == 20:
        return list(deut_18)  # already full-length, pass through

    if len(deut_18) != N_EFFECTIVE_AA:
        raise ValueError(
            f"expand_deuteration_vector expects {N_EFFECTIVE_AA} or 20 elements, "
            f"got {len(deut_18)}"
        )

    deut_20 = [False] * len(AMINO_ACIDS)  # 20 elements, default False
    for gene_idx, group in enumerate(LINKED_AA_GROUPS):
        for aa_code in group:
            aa_idx = AA_DICT[aa_code]  # position in AMINO_ACIDS (0-19)
            deut_20[aa_idx] = deut_18[gene_idx]
    return deut_20


def merge_restrictions_to_18(restrictions_20: Sequence[bool]) -> List[bool]:
    """
    Convert a 20-element restriction vector (one per canonical AA) to an
    18-element vector (one per effective gene / linked group).

    For linked pairs the group is considered modifiable if ANY member AA is
    modifiable in the input (OR logic).

    Args:
        restrictions_20: 20-element boolean list aligned with AMINO_ACIDS.

    Returns:
        18-element boolean list aligned with EFFECTIVE_AMINO_ACIDS.

    Note:
        If an 18-element list is passed it is returned unchanged.
    """
    if len(restrictions_20) == N_EFFECTIVE_AA:
        return list(restrictions_20)  # already 18 elements

    if len(restrictions_20) != len(AMINO_ACIDS):
        raise ValueError(
            f"merge_restrictions_to_18 expects {len(AMINO_ACIDS)} or "
            f"{N_EFFECTIVE_AA} elements, got {len(restrictions_20)}"
        )

    result: List[bool] = []
    for group in LINKED_AA_GROUPS:
        group_modifiable = any(restrictions_20[AA_DICT[code]] for code in group)
        result.append(group_modifiable)
    return result


# ---------------------------------------------------------------------------
# Linked-pair helpers (selection and restrictions)
# ---------------------------------------------------------------------------

# group label -> residue names, e.g. "ASN+ASP" -> ("ASN", "ASP")
RESIDUE_NAMES_BY_GROUP: Dict[str, Tuple[str, ...]] = {
    "+".join(group): tuple(group) for group in LINKED_AA_GROUPS
}

# residue name -> group label, e.g. "ASP" -> "ASN+ASP"
GROUP_LABEL_BY_RESIDUE: Dict[str, str] = {
    code: "+".join(group) for group in LINKED_AA_GROUPS for code in group
}


def linked_partners(code: str) -> Tuple[str, ...]:
    """Return the other member(s) of the linked group of *code* (may be empty)."""
    code = code.strip().upper()
    label = GROUP_LABEL_BY_RESIDUE.get(code)
    if label is None:
        return ()
    return tuple(c for c in RESIDUE_NAMES_BY_GROUP[label] if c != code)


def expand_linked_codes(codes: Iterable[str]) -> Tuple[List[str], List[Tuple[str, str]]]:
    """
    Expand a list of 3-letter codes so that linked pairs stay linked.

    Returns:
        (expanded, pairs) where *expanded* preserves the input order (deduplicated)
        and *pairs* lists the ``(given_code, added_partner)`` additions made by
        the linked rule, for user-facing warnings.
    """
    expanded: List[str] = []
    pairs: List[Tuple[str, str]] = []
    for raw in codes:
        code = raw.strip().upper()
        if not code:
            continue
        if code not in expanded:
            expanded.append(code)
        for partner in linked_partners(code):
            if partner not in expanded:
                expanded.append(partner)
                pairs.append((code, partner))
    return expanded, pairs


def codes_from_chromosome(bits18: Sequence[bool]) -> List[str]:
    """Residue names deuterated by an 18-element chromosome (linked pairs included)."""
    if len(bits18) != N_EFFECTIVE_AA:
        raise ValueError(
            f"codes_from_chromosome expects {N_EFFECTIVE_AA} elements, got {len(bits18)}"
        )
    names: List[str] = []
    for bit, group in zip(bits18, LINKED_AA_GROUPS):
        if bit:
            names.extend(group)
    return names


def chromosome_from_codes(codes: Iterable[str]) -> List[bool]:
    """
    Build an 18-element chromosome vector from 3-letter codes.

    Linked pairs are handled by the group table: setting any member sets the
    whole group (so ASN implies ASP and GLU implies GLN).
    """
    bits = [False] * N_EFFECTIVE_AA
    expanded, _ = expand_linked_codes(codes)
    for code in expanded:
        group_idx = AA_GROUP_INDEX.get(code)
        if group_idx is None:
            raise ValueError(
                f"Unknown amino acid code: {code!r} "
                f"(valid codes: {', '.join(sorted(AA_DICT))})"
            )
        bits[group_idx] = True
    return bits


def restrict_linked_codes(
    restrictions_18: Sequence[bool],
    codes: Iterable[str],
) -> Tuple[List[bool], List[Tuple[str, str]]]:
    """
    Disable the effective gene of every code in *codes* (linked rule applied).

    Args:
        restrictions_18: current 18-element restriction list (True = modifiable).
        codes: 3-letter codes to restrict (forbid).

    Returns:
        (updated_restrictions, pairs) with *pairs* the ``(given, partner)``
        additions caused by the linked rule.
    """
    updated = list(restrictions_18)
    expanded, pairs = expand_linked_codes(codes)
    for code in expanded:
        group_idx = AA_GROUP_INDEX[code]
        updated[group_idx] = False
    return updated, pairs
