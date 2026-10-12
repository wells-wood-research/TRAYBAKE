"""drMD-compatible custom selections for choosing TRAYBAKE targets."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from traybake.pdb import Atom, PDBStructure, ResidueKey

_CUSTOM_FIELDS = ("CHAIN_ID", "RES_NAME", "RES_ID", "ATOM_NAME")


def load_selections(path: str | Path) -> list[dict[str, Any]]:
    """Load a list of drMD-style custom-selection clauses from YAML."""
    selection_path = Path(path)
    try:
        with selection_path.open(encoding="utf-8") as handle:
            document = yaml.safe_load(handle)
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise ValueError(f"Could not read selections file {selection_path}: {error}") from error

    if not isinstance(document, list):
        raise ValueError("Selections YAML must contain a list of custom-selection mappings")
    for index, clause in enumerate(document, start=1):
        _validate_clause(clause, index)
    return document


def selected_residues(
    structure: PDBStructure, selections: list[dict[str, Any]]
) -> set[ResidueKey]:
    """Return residues containing atoms matched by any custom-selection clause."""
    selected: set[ResidueKey] = set()
    for atom in structure.atoms:
        if any(_matches(atom, clause) for clause in selections):
            selected.add(atom.key)
    return selected


def _validate_clause(clause: Any, index: int) -> None:
    if not isinstance(clause, dict):
        raise ValueError(f"Selection {index} must be a mapping")
    missing = set(_CUSTOM_FIELDS) - set(clause)
    if missing:
        raise ValueError(f"Selection {index} is missing: {', '.join(sorted(missing))}")
    for field in _CUSTOM_FIELDS:
        value = clause[field]
        allowed = (int,) if field == "RES_ID" else (str,)
        if value == "_":
            continue
        values = value if isinstance(value, list) else [value]
        if not values or any(
            isinstance(item, bool) or not isinstance(item, allowed) for item in values
        ):
            expected = "integer" if field == "RES_ID" else "string"
            raise ValueError(
                f"Selection {index} {field} must be '_', an {expected}, "
                f"or a list of {expected}s"
            )


def _matches(atom: Atom, clause: dict[str, Any]) -> bool:
    values = {
        "CHAIN_ID": atom.key.chain,
        "RES_NAME": atom.key.name,
        "RES_ID": atom.key.number,
        "ATOM_NAME": atom.name,
    }
    return all(
        clause[field] == "_"
        or values[field] in (clause[field] if isinstance(clause[field], list) else [clause[field]])
        for field in _CUSTOM_FIELDS
    )
