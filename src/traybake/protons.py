"""Add residue-local hydrogen coordinates to fixed-coordinate PDB bubbles."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from math import acos, atan2, cos, sin, sqrt
from typing import Iterable

from traybake.geometry import group_residues
from traybake.pdb import Atom, PDBStructure, ResidueKey

Vector = tuple[float, float, float]


def add_residue_protons(atoms: tuple[Atom, ...]) -> tuple[Atom, ...]:
    """Add implicit hydrogens using separately embedded residue templates.

    Each residue is embedded and locally optimized as an independent RDKit
    molecule. For every hydrogen required by the complete bubble molecule, its
    bond length, polar angle, and azimuth in the embedded residue's local frame
    are transferred onto the original PDB heavy-atom coordinates. Existing
    heavy atoms and explicit cap hydrogens are left untouched.
    """
    if not atoms:
        return atoms

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as error:
        raise RuntimeError("RDKit is required to add residue-local protons") from error

    groups = group_residues(PDBStructure(atoms))
    whole_mol = _molecule_from_atoms(Chem, atoms)
    if whole_mol is None:
        raise ValueError("RDKit could not read the capped bubble PDB")
    hydrogenated = Chem.AddHs(whole_mol, addCoords=True)
    conformer = hydrogenated.GetConformer()
    original_atom_count = whole_mol.GetNumAtoms()

    key_by_identity = {
        (key.chain, key.number, key.insertion_code, key.name): key for key in groups
    }
    required: dict[tuple[ResidueKey, str], list[Vector]] = defaultdict(list)
    for hydrogen in hydrogenated.GetAtoms():
        if hydrogen.GetIdx() < original_atom_count or hydrogen.GetSymbol() != "H":
            continue
        parent, = hydrogen.GetNeighbors()
        info = parent.GetPDBResidueInfo()
        if info is None:
            continue
        identity = (
            info.GetChainId().strip(),
            info.GetResidueNumber(),
            info.GetInsertionCode().strip(),
            info.GetResidueName().strip().upper(),
        )
        key = key_by_identity.get(identity)
        if key is None:
            continue
        required[(key, info.GetName().strip())].append(
            _position(conformer, hydrogen.GetIdx())
        )

    new_atoms: list[Atom] = []
    serial = max((_atom_serial(atom) for atom in atoms), default=0)
    for key, residue_atoms in groups.items():
        residue_hydrogens = {
            atom_name: positions
            for (residue_key, atom_name), positions in required.items()
            if residue_key == key
        }
        if not residue_hydrogens:
            continue

        template = _molecule_from_atoms(Chem, residue_atoms)
        if template is None:
            raise ValueError(f"RDKit could not build residue template for {key}")
        original_template_count = template.GetNumAtoms()
        template = Chem.AddHs(template, addCoords=True)
        seed = _residue_seed(key)
        embedded = AllChem.EmbedMolecule(
            template,
            useRandomCoords=True,
            ignoreSmoothingFailures=True,
            randomSeed=seed,
        )
        if embedded != 0:
            raise ValueError(f"RDKit could not embed residue {key}")
        if AllChem.MMFFHasAllMoleculeParams(template):
            AllChem.MMFFOptimizeMolecule(template, maxIters=200)

        template_conf = template.GetConformer()
        template_names = {
            info.GetName().strip(): atom.GetIdx()
            for atom in template.GetAtoms()
            if atom.GetAtomicNum() > 1
            and (info := atom.GetPDBResidueInfo()) is not None
        }
        source_by_name = {atom.name: atom for atom in residue_atoms}
        used_names = {atom.name for atom in residue_atoms}

        template_hydrogens: dict[str, list[int]] = defaultdict(list)
        for hydrogen in template.GetAtoms():
            if hydrogen.GetIdx() < original_template_count or hydrogen.GetSymbol() != "H":
                continue
            parent, = hydrogen.GetNeighbors()
            info = parent.GetPDBResidueInfo()
            if info is not None:
                template_hydrogens[info.GetName().strip()].append(hydrogen.GetIdx())

        for parent_name, fallback_positions in residue_hydrogens.items():
            parent_idx = template_names.get(parent_name)
            parent_source = source_by_name.get(parent_name)
            if parent_idx is None or parent_source is None:
                raise ValueError(f"No template atom {parent_name!r} for {key}")

            anchors = _anchor_names(
                template,
                template_conf,
                parent_idx,
                template_names,
                source_by_name,
            )
            candidates: list[Vector] = []
            for hydrogen_idx in template_hydrogens.get(parent_name, ()):
                if anchors is None:
                    break
                internal = _internal_coordinates(
                    template_conf,
                    parent_idx,
                    hydrogen_idx,
                    template_names[anchors[0]],
                    template_names[anchors[1]],
                )
                target_anchors = (
                    source_by_name[anchors[0]].coordinates,
                    source_by_name[anchors[1]].coordinates,
                )
                candidates.append(
                    _place_internal_coordinates(
                        parent_source.coordinates, target_anchors, internal
                    )
                )

            count = len(fallback_positions)
            candidates = _select_hydrogen_positions(
                candidates, count, parent_source.coordinates, residue_atoms
            )
            fallback = fallback_positions[:count]
            if len(candidates) >= count and len(fallback) == count:
                if _hydrogen_clash_score(candidates, parent_source.coordinates, residue_atoms) > _hydrogen_clash_score(
                    fallback, parent_source.coordinates, residue_atoms
                ):
                    candidates = fallback
            if len(candidates) < count:
                candidates.extend(fallback[len(candidates) : count])

            names = _hydrogen_names(parent_name, count, used_names)
            for name, coordinates in zip(names, candidates):
                serial += 1
                new_atoms.append(
                    _make_hydrogen(parent_source, name, coordinates, serial)
                )

    return (*atoms, *new_atoms)


def _molecule_from_atoms(Chem, atoms: Iterable[Atom]):
    block = "\n".join(atom.line for atom in atoms) + "\nEND\n"
    return Chem.MolFromPDBBlock(
        block, sanitize=True, removeHs=False, proximityBonding=True
    )


def _position(conformer, atom_index: int) -> Vector:
    point = conformer.GetAtomPosition(atom_index)
    return point.x, point.y, point.z


def _residue_seed(key: ResidueKey) -> int:
    seed = key.model * 1009 + key.number * 917 + sum(map(ord, key.chain + key.name))
    return seed % 2_000_000_000 + 1


def _anchor_names(
    molecule,
    conformer,
    parent_idx: int,
    atom_names: dict[str, int],
    source_by_name: dict[str, Atom],
) -> tuple[str, str] | None:
    parent_name = next(name for name, idx in atom_names.items() if idx == parent_idx)
    parent_template = _position(conformer, parent_idx)
    parent_source = source_by_name[parent_name].coordinates
    heavy_names = list(atom_names)
    direct_names = [
        neighbor.GetPDBResidueInfo().GetName().strip()
        for neighbor in molecule.GetAtomWithIdx(parent_idx).GetNeighbors()
        if neighbor.GetAtomicNum() > 1 and neighbor.GetPDBResidueInfo() is not None
    ]
    if len(direct_names) >= 2:
        candidates = list(combinations(direct_names, 2))
    elif direct_names:
        candidates = [
            (direct_names[0], name)
            for name in heavy_names
            if name not in {parent_name, direct_names[0]}
        ]
    else:
        candidates = [
            pair for pair in combinations(heavy_names, 2) if parent_name not in pair
        ]

    best_pair = None
    best_score = 0.0
    for first, second in candidates:
        template_score = _frame_area(
            parent_template,
            _position(conformer, atom_names[first]),
            _position(conformer, atom_names[second]),
        )
        source_score = _frame_area(
            parent_source,
            source_by_name[first].coordinates,
            source_by_name[second].coordinates,
        )
        score = min(template_score, source_score)
        if score > best_score:
            best_pair, best_score = (first, second), score
    return best_pair


def _frame_area(origin: Vector, first: Vector, second: Vector) -> float:
    first_vector = _subtract(first, origin)
    second_vector = _subtract(second, origin)
    denominator = _length(first_vector) * _length(second_vector)
    if denominator < 1e-10:
        return 0.0
    return _length(_cross(first_vector, second_vector)) / denominator


def _internal_coordinates(
    conformer, parent_idx: int, hydrogen_idx: int, first_idx: int, second_idx: int
) -> tuple[float, float, float]:
    parent = _position(conformer, parent_idx)
    first = _position(conformer, first_idx)
    second = _position(conformer, second_idx)
    basis = _local_basis(parent, first, second)
    hydrogen_vector = _subtract(_position(conformer, hydrogen_idx), parent)
    bond_length = _length(hydrogen_vector)
    if bond_length < 1e-10:
        raise ValueError("Embedded residue has a zero-length hydrogen bond")
    axial = max(-1.0, min(1.0, _dot(hydrogen_vector, basis[0]) / bond_length))
    polar_angle = acos(axial)
    azimuth = atan2(
        _dot(hydrogen_vector, basis[2]), _dot(hydrogen_vector, basis[1])
    )
    return bond_length, polar_angle, azimuth


def _place_internal_coordinates(
    parent: Vector, anchors: tuple[Vector, Vector], internal: tuple[float, float, float]
) -> Vector:
    basis = _local_basis(parent, anchors[0], anchors[1])
    bond_length, polar_angle, azimuth = internal
    direction = _add(
        _scale(basis[0], cos(polar_angle)),
        _scale(
            _add(
                _scale(basis[1], cos(azimuth)),
                _scale(basis[2], sin(azimuth)),
            ),
            sin(polar_angle),
        ),
    )
    return _add(parent, _scale(direction, bond_length))


def _local_basis(origin: Vector, first: Vector, second: Vector) -> tuple[Vector, Vector, Vector]:
    e1 = _unit(_subtract(first, origin))
    second_vector = _subtract(second, origin)
    in_plane = _subtract(second_vector, _scale(e1, _dot(second_vector, e1)))
    e2 = _unit(in_plane)
    e3 = _unit(_cross(e1, e2))
    return e1, e2, e3


def _select_hydrogen_positions(
    candidates: list[Vector],
    count: int,
    parent: Vector,
    residue_atoms: tuple[Atom, ...],
) -> list[Vector]:
    if len(candidates) <= count:
        return candidates
    heavy_coordinates = [
        atom.coordinates for atom in residue_atoms if atom.element != "H"
    ]
    best = min(
        combinations(candidates, count),
        key=lambda group: sum(
            max(0.0, 1.7 - _length(_subtract(hydrogen, other))) ** 2
            for hydrogen in group
            for other in heavy_coordinates
            if _length(_subtract(hydrogen, other)) > 1.15
        )
        + sum(
            max(0.0, 1.5 - _length(_subtract(left, right))) ** 2
            for left, right in combinations(group, 2)
        ),
    )
    return list(best)


def _hydrogen_clash_score(
    positions: list[Vector], parent: Vector, residue_atoms: tuple[Atom, ...]
) -> float:
    heavy = [
        atom.coordinates
        for atom in residue_atoms
        if atom.element != "H" and _length(_subtract(atom.coordinates, parent)) > 1e-6
    ]
    return sum(
        max(0.0, 1.7 - _length(_subtract(hydrogen, other))) ** 2
        for hydrogen in positions
        for other in heavy
    )


def _hydrogen_names(parent_name: str, count: int, used: set[str]) -> list[str]:
    if parent_name == "N":
        base = "H"
        suffixes = [""] if count == 1 else [str(index) for index in range(1, count + 1)]
    else:
        base = "H" + parent_name[1:]
        start = 1 if parent_name[0] in {"N", "O", "S"} else 2
        suffixes = [""] if count == 1 else [str(start + index) for index in range(count)]

    names = []
    for suffix in suffixes:
        candidate = (base + suffix)[:4]
        serial = 1
        while candidate in used:
            candidate = (base[:3] + str(serial))[:4]
            serial += 1
        names.append(candidate)
        used.add(candidate)
    return names


def _make_hydrogen(
    parent: Atom, name: str, coordinates: Vector, serial: int
) -> Atom:
    x, y, z = coordinates
    atom_name = f"{name:>4s}"
    key = parent.key
    line = (
        f"ATOM  {serial % 100000:5d} {atom_name} {key.name:>3s} "
        f"{(key.chain or ' '):1s}{key.number:4d}{(key.insertion_code or ' '):1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{parent.occupancy:6.2f}{0.0:6.2f}"
        f"          {'H':>2s}"
    )
    return Atom(key, name, "H", x, y, z, parent.occupancy, "", line)


def _atom_serial(atom: Atom) -> int:
    try:
        return int(atom.line[6:11])
    except ValueError:
        return 0


def _subtract(left: Vector, right: Vector) -> Vector:
    return tuple(a - b for a, b in zip(left, right))


def _add(left: Vector, right: Vector) -> Vector:
    return tuple(a + b for a, b in zip(left, right))


def _scale(vector: Vector, factor: float) -> Vector:
    return tuple(value * factor for value in vector)


def _dot(left: Vector, right: Vector) -> float:
    return sum(a * b for a, b in zip(left, right))


def _cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _length(vector: Vector) -> float:
    return sqrt(_dot(vector, vector))


def _unit(vector: Vector) -> Vector:
    length = _length(vector)
    if length < 1e-10:
        raise ValueError("Cannot define residue-local hydrogen coordinates")
    return _scale(vector, 1.0 / length)
