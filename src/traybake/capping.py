"""Build bubble fragments with proton caps at exposed covalent cuts."""

from __future__ import annotations

from dataclasses import dataclass
from math import dist, sqrt

from traybake.geometry import Bubble, group_residues
from traybake.pdb import Atom, PDBStructure, ResidueKey

BACKBONE_COMPLETE = frozenset({"N", "CA", "C", "O"})
BACKBONE_ONLY = frozenset({
    "N", "CA", "C", "O", "OXT", "H", "HN", "H1", "H2", "H3",
    "HT", "HT1", "HT2", "HT3", "HA", "HA2", "HA3", "1H", "2H", "3H",
})
PEPTIDE_BOND_MAX = 1.8
CAP_BOND_LENGTH = 1.01
CARBON_CAP_BOND_LENGTH = 1.09
SIDECHAIN_CAP_BOND_LENGTH = 1.0


@dataclass(frozen=True)
class CapRecord:
    residue: ResidueKey
    atom_name: str
    attached_to: str
    coordinates: tuple[float, float, float]
    kind: str


@dataclass(frozen=True)
class CappedBubble:
    atoms: tuple[Atom, ...]
    residue_modes: tuple[tuple[ResidueKey, str], ...]
    caps: tuple[CapRecord, ...]


def build_capped_bubble(structure: PDBStructure, bubble: Bubble, radius: float) -> CappedBubble:
    """Extract a bubble, trimming incomplete backbones and proton-capping cuts.

    A residue keeps all its atoms only when N, CA, C, and O are present and
    every one lies inside the bubble. Otherwise, only CB and downstream
    sidechain atoms are retained, and CB receives a proton in place of CA.
    Exposed peptide N and C atoms receive one proton each.
    """
    if radius <= 0:
        raise ValueError("radius must be greater than zero")

    grouped = group_residues(structure)
    selected: dict[ResidueKey, tuple[Atom, ...]] = {}
    modes: list[tuple[ResidueKey, str]] = []
    for key in bubble.residues:
        atoms = grouped[key]
        by_name = {atom.name: atom for atom in atoms}
        complete_backbone = BACKBONE_COMPLETE <= by_name.keys() and all(
            dist(bubble.center_of_mass, by_name[name].coordinates) <= radius
            for name in BACKBONE_COMPLETE
        )
        if complete_backbone:
            selected[key] = atoms
            modes.append((key, "full"))
            continue

        sidechain = tuple(atom for atom in atoms if atom.name not in BACKBONE_ONLY)
        if "CB" not in {atom.name for atom in sidechain}:
            selected[key] = ()
            modes.append((key, "omitted_no_cb"))
            continue
        selected[key] = sidechain
        modes.append((key, "sidechain_only"))

    atoms = [atom for key in bubble.residues for atom in selected[key]]
    caps: list[CapRecord] = []
    serial = max((_serial_number(atom) for atom in structure.atoms), default=0)

    # Replace the omitted CA-CB bond with one C-H cap on each sidechain-only CB.
    for key, mode in modes:
        if mode != "sidechain_only":
            continue
        by_name = {atom.name: atom for atom in grouped[key]}
        if "CA" not in by_name:
            raise ValueError(f"Cannot place CB cap for {key}: CA is missing")
        cb = by_name["CB"]
        ca = by_name["CA"]
        # Replace the omitted CA with H, then slide H toward CB until the
        # replacement bond is exactly 1.0 Å long.
        coord = _move_to_bond_length(ca.coordinates, cb.coordinates, SIDECHAIN_CAP_BOND_LENGTH)
        cap = _make_cap(key, "HCB", "CB", coord, "sidechain", serial + 1)
        serial += 1
        atoms.append(cap)
        caps.append(CapRecord(key, cap.name, "CB", coord, "sidechain"))

    # Identify covalently connected neighboring residues in the source chain.
    links = _peptide_links(grouped)
    full_keys = {key for key, mode in modes if mode == "full"}
    for key in sorted(full_keys):
        by_name = {atom.name: atom for atom in selected[key]}
        previous_links = {left for left, right in links if right == key}
        next_links = {right for left, right in links if left == key}

        if "N" in by_name and not (previous_links & full_keys):
            coord = _n_cap_position(grouped[key])
            name = _available_name("HNT", selected[key])
            cap = _make_cap(key, name, "N", coord, "N_terminus", serial + 1)
            serial += 1
            atoms.append(cap)
            caps.append(CapRecord(key, cap.name, "N", coord, "N_terminus"))

        if "C" in by_name and not (next_links & full_keys):
            coord = _carbon_cap_position(grouped[key])
            name = _available_name("HCT", selected[key])
            cap = _make_cap(key, name, "C", coord, "C_terminus", serial + 1)
            serial += 1
            atoms.append(cap)
            caps.append(CapRecord(key, cap.name, "C", coord, "C_terminus"))

    return CappedBubble(tuple(atoms), tuple(modes), tuple(caps))


def _peptide_links(grouped: dict[ResidueKey, tuple[Atom, ...]]) -> set[tuple[ResidueKey, ResidueKey]]:
    chains: dict[tuple[int, str], list[ResidueKey]] = {}
    for key in grouped:
        chains.setdefault((key.model, key.chain), []).append(key)
    links: set[tuple[ResidueKey, ResidueKey]] = set()
    for keys in chains.values():
        ordered = sorted(keys, key=lambda key: (key.number, key.insertion_code, key.name))
        for left, right in zip(ordered, ordered[1:]):
            c_atom = next((atom for atom in grouped[left] if atom.name == "C"), None)
            n_atom = next((atom for atom in grouped[right] if atom.name == "N"), None)
            if c_atom and n_atom and dist(c_atom.coordinates, n_atom.coordinates) <= PEPTIDE_BOND_MAX:
                links.add((left, right))
    return links


def _n_cap_position(atoms: tuple[Atom, ...]) -> tuple[float, float, float]:
    by_name = {atom.name: atom for atom in atoms}
    if "N" not in by_name or "CA" not in by_name:
        raise ValueError(f"Cannot place N-terminal cap for {atoms[0].key}: N or CA is missing")
    n_coord = by_name["N"].coordinates
    ca_coord = by_name["CA"].coordinates
    toward_ca = _unit(_subtract(ca_coord, n_coord), "N-CA vector")
    if "CB" in by_name:
        ca_to_cb = _subtract(by_name["CB"].coordinates, ca_coord)
    else:
        # Glycine has no CB; choose a stable perpendicular for its terminal cap.
        ca_to_cb = (1.0, 0.0, 0.0) if abs(toward_ca[0]) < 0.9 else (0.0, 1.0, 0.0)
    projection = _subtract(ca_to_cb, _scale(toward_ca, _dot(ca_to_cb, toward_ca)))
    perpendicular = _unit(projection, "N-terminal cap plane")
    direction = _unit(
        _add(_scale(toward_ca, -0.5), _scale(perpendicular, sqrt(3.0) / 2.0)),
        "N-terminal cap direction",
    )
    return _add(n_coord, _scale(direction, CAP_BOND_LENGTH))


def _carbon_cap_position(atoms: tuple[Atom, ...]) -> tuple[float, float, float]:
    by_name = {atom.name: atom for atom in atoms}
    if not {"C", "CA", "O"} <= by_name.keys():
        raise ValueError(f"Cannot place C-terminal cap for {atoms[0].key}: C, CA, or O is missing")
    c_coord = by_name["C"].coordinates
    toward_ca = _unit(_subtract(by_name["CA"].coordinates, c_coord), "C-CA vector")
    toward_o = _unit(_subtract(by_name["O"].coordinates, c_coord), "C-O vector")
    direction = _unit(_scale(_add(toward_ca, toward_o), -1.0), "C-terminal cap direction")
    return _add(c_coord, _scale(direction, CARBON_CAP_BOND_LENGTH))


def _make_cap(
    key: ResidueKey,
    name: str,
    attached_to: str,
    coord: tuple[float, float, float],
    kind: str,
    serial: int,
) -> Atom:
    atom_name = name[:4]
    line = (
        f"ATOM  {serial % 100000:5d} {atom_name:>4s} {key.name:>3s} "
        f"{(key.chain or ' '):1s}{key.number:4d}{(key.insertion_code or ' '):1s}   "
        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{0.0:6.2f}"
        f"          {'H':>2s}"
    )
    return Atom(key, atom_name, "H", *coord, 1.0, "", line)


def _available_name(preferred: str, atoms: tuple[Atom, ...]) -> str:
    existing = {atom.name for atom in atoms}
    if preferred not in existing:
        return preferred
    for suffix in range(1, 10):
        candidate = f"{preferred[:3]}{suffix}"
        if candidate not in existing:
            return candidate
    raise ValueError(f"No available hydrogen name for {atoms[0].key}")


def _serial_number(atom: Atom) -> int:
    try:
        return int(atom.line[6:11])
    except ValueError:
        return 0


def _subtract(left, right):
    return tuple(a - b for a, b in zip(left, right))


def _add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def _scale(vector, factor):
    return tuple(value * factor for value in vector)


def _dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def _unit(vector, label: str):
    length = sum(value * value for value in vector) ** 0.5
    if length < 1e-10:
        raise ValueError(f"Cannot determine cap geometry: degenerate {label}")
    return tuple(value / length for value in vector)


def _move_to_bond_length(
    seed: tuple[float, float, float],
    anchor: tuple[float, float, float],
    bond_length: float,
) -> tuple[float, float, float]:
    """Move a seed atom along seed→anchor to a fixed distance from anchor."""
    direction = _unit(_subtract(anchor, seed), "replacement-atom direction")
    return _subtract(anchor, _scale(direction, bond_length))
