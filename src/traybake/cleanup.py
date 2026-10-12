"""Remove non-protein residues before building local neighborhoods."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from traybake.geometry import group_residues
from traybake.pdb import PDBStructure, ResidueKey

# Standard amino-acid residue names recognized as protein for this initial
# cleanup stage. Ligands, cofactors, waters, nucleic acids, and other residue
# types are excluded; future chemistry-aware handling can retain selected ones.
PROTEIN_RESIDUES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})


@dataclass(frozen=True)
class CleanupResult:
    structure: PDBStructure
    removed_residues: tuple[ResidueKey, ...]


def remove_nonprotein(structure: PDBStructure) -> CleanupResult:
    """Keep canonical amino-acid residues and remove every other residue."""
    residue_groups = group_residues(structure)
    kept_keys = {key for key in residue_groups if key.name in PROTEIN_RESIDUES}
    cleaned = PDBStructure(tuple(atom for atom in structure.atoms if atom.key in kept_keys))
    removed = tuple(sorted(set(residue_groups) - kept_keys))
    return CleanupResult(cleaned, removed)


def write_cleaned_pdb(
    input_path: str | Path, output_path: str | Path, cleanup: CleanupResult
) -> Path:
    """Write the cleaned atom set while retaining useful PDB header/chain records."""
    source_path = Path(input_path)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    retained_lines = Counter(atom.line for atom in cleanup.structure.atoms)
    has_end = False
    with source_path.open(encoding="ascii", errors="replace", newline="") as source:
        lines = source.readlines()
    with destination.open("w", encoding="ascii", newline="") as output:
        for raw_line in lines:
            line = raw_line.rstrip("\r\n")
            record = line[:6].strip()
            if record in {"ATOM", "HETATM"}:
                if retained_lines[line] > 0:
                    retained_lines[line] -= 1
                    output.write(raw_line)
            elif record in {
                "ANISOU",
                "CONECT",
                "FORMUL",
                "HET",
                "HETNAM",
                "HETSYN",
                "MASTER",
            }:
                # These records refer to discarded atoms/residues or source atom counts.
                continue
            else:
                output.write(raw_line)
                has_end = has_end or record == "END"
        if not has_end:
            output.write("END\n")
    return destination
