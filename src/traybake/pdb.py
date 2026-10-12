"""Minimal, dependency-free parsing helpers for fixed-column PDB files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, order=True)
class ResidueKey:
    model: int
    chain: str
    number: int
    insertion_code: str
    name: str

    def label(self) -> str:
        chain = self.chain or "_"
        insertion = self.insertion_code.strip()
        suffix = f"{self.number}{insertion}"
        return f"model{self.model}_{chain}_{self.name}_{suffix}"


@dataclass(frozen=True)
class Atom:
    key: ResidueKey
    name: str
    element: str
    x: float
    y: float
    z: float
    occupancy: float
    altloc: str
    line: str

    @property
    def coordinates(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


@dataclass(frozen=True)
class PDBStructure:
    atoms: tuple[Atom, ...]


def read_pdb(path: str | Path) -> PDBStructure:
    """Read ATOM/HETATM records, retaining coordinates and source lines.

    For alternate locations of the same atom, blank locations take precedence;
    otherwise the highest-occupancy location is kept. MODEL records are
    represented in residue keys, so identical chain/residue IDs in separate
    models do not collide.
    """
    atoms: list[Atom] = []
    model = 1
    seen_model = False
    with Path(path).open(encoding="ascii", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            record = line[:6].strip()
            if record == "MODEL":
                try:
                    model = int(line[10:14].strip())
                except ValueError:
                    model = 1
                seen_model = True
                continue
            if record == "ENDMDL":
                continue
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                key = ResidueKey(
                    model=model if seen_model else 1,
                    chain=line[21:22].strip(),
                    number=int(line[22:26]),
                    insertion_code=line[26:27].strip(),
                    name=line[17:20].strip().upper(),
                )
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper()
                if not element:
                    element = "".join(char for char in atom_name if char.isalpha())[:1].upper()
                occupancy_text = line[54:60].strip()
                occupancy = float(occupancy_text) if occupancy_text else 1.0
                atom = Atom(
                    key=key,
                    name=atom_name,
                    element=element,
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=occupancy,
                    altloc=line[16:17].strip(),
                    line=line,
                )
            except (ValueError, IndexError) as error:
                raise ValueError(f"Invalid atom record at {path}:{line_number}: {line!r}") from error
            atoms.append(atom)

    # Prefer the shared (blank altloc) atom; for nonblank alternatives keep
    # only the highest-occupancy record for each model/residue/atom name.
    selected: dict[tuple[ResidueKey, str], Atom] = {}
    for atom in atoms:
        identity = atom.key, atom.name
        current = selected.get(identity)
        if current is None or (current.altloc != "" and (atom.altloc == "" or atom.occupancy > current.occupancy)):
            selected[identity] = atom
    return PDBStructure(tuple(selected.values()))


def write_residues(path: str | Path, residues: set[ResidueKey], structure: PDBStructure) -> None:
    """Write selected source atom records as a standalone PDB fragment."""
    with Path(path).open("w", encoding="ascii") as handle:
        for atom in structure.atoms:
            if atom.key in residues:
                handle.write(atom.line + "\n")
        handle.write("END\n")


def write_atoms(path: str | Path, atoms: tuple[Atom, ...]) -> None:
    """Write selected source and generated atom records as a standalone PDB."""
    with Path(path).open("w", encoding="ascii") as handle:
        for atom in atoms:
            handle.write(atom.line + "\n")
        handle.write("END\n")
