"""Titratable-residue selection and local neighborhood geometry."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import dist

from traybake.pdb import Atom, PDBStructure, ResidueKey

TITRATABLE_RESIDUES = frozenset({"ASP", "GLU", "LYS", "ARG", "HIS"})
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})


@dataclass(frozen=True)
class Bubble:
    center: ResidueKey
    center_of_mass: tuple[float, float, float]
    residues: tuple[ResidueKey, ...]


def group_residues(structure: PDBStructure) -> dict[ResidueKey, tuple[Atom, ...]]:
    grouped: dict[ResidueKey, list[Atom]] = defaultdict(list)
    for atom in structure.atoms:
        grouped[atom.key].append(atom)
    return {key: tuple(atoms) for key, atoms in grouped.items()}


def sidechain_center(atoms: tuple[Atom, ...]) -> tuple[float, float, float]:
    """Return the unweighted centroid of a residue's observed heavy sidechain atoms."""
    sidechain = [
        atom for atom in atoms
        if atom.name not in BACKBONE_ATOMS and atom.element not in {"H", "D"}
    ]
    if not sidechain:
        key = atoms[0].key if atoms else "unknown residue"
        raise ValueError(f"No observed heavy sidechain atoms for {key}")
    count = len(sidechain)
    return tuple(sum(getattr(atom, axis) for atom in sidechain) / count for axis in ("x", "y", "z"))


def build_bubbles(
    structure: PDBStructure,
    radius: float = 5.0,
    target_residues: set[ResidueKey] | None = None,
) -> tuple[Bubble, ...]:
    """Build a bubble for each titratable residue.

    A residue is included when any observed atom is within or exactly on the
    radius from the target's unweighted heavy-sidechain-atom centroid. All
    parsed residues, including HETATM records, are eligible for inclusion.
    """
    if radius <= 0:
        raise ValueError("radius must be greater than zero")
    grouped = group_residues(structure)
    bubbles = []
    for key in sorted(grouped):
        if key.name not in TITRATABLE_RESIDUES or (
            target_residues is not None and key not in target_residues
        ):
            continue
        center = sidechain_center(grouped[key])
        members = tuple(
            candidate
            for candidate in sorted(grouped)
            if any(dist(center, atom.coordinates) <= radius for atom in grouped[candidate])
        )
        bubbles.append(Bubble(key, center, members))
    return tuple(bubbles)
