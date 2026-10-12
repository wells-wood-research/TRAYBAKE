"""Rename protonation states and prepare a PDB with Amber LEaP."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any


# GLUP is the tally's protonated-glutamate label. Amber's ff19SB library uses
# the three-character residue name GLH, which fits the PDB residue-name field.
AMBER_PDB_RESNAMES = {
    "ASP": "ASP",
    "ASH": "ASH",
    "GLU": "GLU",
    "GLUP": "GLH",
    "GLH": "GLH",
    "LYS": "LYS",
    "LYN": "LYN",
    "ARG": "ARG",
    "AR0": "AR0",
    "HID": "HID",
    "HIE": "HIE",
    "HIP": "HIP",
}


def select_residue_assignments(
    counts: Mapping[str, Mapping[str, int]],
) -> list[dict[str, Any]]:
    """Choose each residue's unique highest-weight assignment.

    The count-file key format is ``CHAIN_ORIGINAL-NAME_RESID[ICODE]``. Since
    the current key format omits PDB model number, callers must reject input
    structures where a key would identify residues in multiple models.
    """
    assignments = []
    seen: set[tuple[str, str, int, str]] = set()
    for label, state_counts in sorted(counts.items()):
        chain, original_name, residue_number, insertion_code = _parse_label(label)
        identity = (chain, original_name, residue_number, insertion_code)
        if identity in seen:
            raise ValueError(f"Duplicate residue assignment for {label}")
        seen.add(identity)
        if not state_counts:
            raise ValueError(f"No state counts were supplied for {label}")
        highest_count = max(state_counts.values())
        winners = [state for state, count in state_counts.items() if count == highest_count]
        if len(winners) != 1:
            raise ValueError(
                f"Residue {label} has tied assignments at weight {highest_count}: "
                f"{', '.join(sorted(winners))}"
            )
        assignment = winners[0]
        try:
            amber_resname = AMBER_PDB_RESNAMES[assignment]
        except KeyError as error:
            raise ValueError(
                f"Unsupported Amber residue assignment {assignment!r} for {label}"
            ) from error
        assignments.append(
            {
                "label": label,
                "chain": chain,
                "original_resname": original_name,
                "residue_number": residue_number,
                "insertion_code": insertion_code,
                "weighted_counts": dict(state_counts),
                "assignment": assignment,
                "amber_resname": amber_resname,
            }
        )
    return assignments


def write_renamed_pdb(
    input_pdb: str | Path,
    output_pdb: str | Path,
    assignments: list[dict[str, Any]],
) -> dict[str, Any]:
    """Copy a PDB, changing only assigned residue-name columns.

    All source records, including waters, remarks, TER records, and coordinates,
    are retained. The input file is never modified.
    """
    input_path = Path(input_pdb).resolve()
    output_path = Path(output_pdb).resolve()
    if input_path == output_path:
        raise ValueError("The renamed PDB must not overwrite the input PDB")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    by_identity = {
        (
            item["chain"],
            item["original_resname"],
            item["residue_number"],
            item["insertion_code"],
        ): item
        for item in assignments
    }
    matches: dict[tuple[str, str, int, str], set[int]] = {
        identity: set() for identity in by_identity
    }
    model = 1
    seen_models = False
    with input_path.open(encoding="ascii", errors="replace", newline="") as source:
        lines = source.readlines()

    renamed_lines = []
    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip()
        if record == "MODEL":
            seen_models = True
            try:
                model = int(line[10:14].strip())
            except ValueError as error:
                raise ValueError(
                    f"Invalid MODEL record at {input_path}:{line_number}"
                ) from error
            renamed_lines.append(line)
            continue
        if record == "ENDMDL":
            renamed_lines.append(line)
            continue
        if record not in {"ATOM", "HETATM"} or len(line) < 27:
            renamed_lines.append(line)
            continue

        try:
            chain = line[21:22].strip()
            original_name = line[17:20].strip().upper()
            residue_number = int(line[22:26])
            insertion_code = line[26:27].strip()
        except ValueError as error:
            raise ValueError(
                f"Invalid residue identity at {input_path}:{line_number}: {line.rstrip()}"
            ) from error
        identity = (chain, original_name, residue_number, insertion_code)
        assignment = by_identity.get(identity)
        if assignment is None:
            renamed_lines.append(line)
            continue
        matches[identity].add(model if seen_models else 1)
        replacement = assignment["amber_resname"]
        renamed_lines.append(line[:17] + replacement.rjust(3) + line[20:])

    missing = [
        item["label"] for identity, item in by_identity.items() if not matches[identity]
    ]
    if missing:
        raise ValueError(
            "Assigned residues were not found in the input PDB: "
            + ", ".join(sorted(missing))
        )
    repeated_models = [
        by_identity[identity]["label"]
        for identity, model_ids in matches.items()
        if len(model_ids) > 1
    ]
    if repeated_models:
        raise ValueError(
            "Assignment keys do not distinguish PDB models for: "
            + ", ".join(sorted(repeated_models))
        )

    with output_path.open("w", encoding="ascii", newline="") as output:
        output.writelines(renamed_lines)
    return {
        "path": str(output_path),
        "renamed_residue_count": sum(
            item["original_resname"] != item["amber_resname"]
            for item in assignments
        ),
        "assignments": assignments,
    }


def run_tleap(
    input_pdb: str | Path,
    counts_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Run Amber LEaP on a residue-renamed copy of the input PDB."""
    if shutil.which("tleap") is None:
        raise FileNotFoundError("TLEaP executable 'tleap' was not found on PATH")

    counts_file = Path(counts_path)
    try:
        counts = json.loads(counts_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Could not read assignment counts from {counts_file}: {error}") from error
    if not isinstance(counts, dict):
        raise ValueError(f"Assignment counts in {counts_file} must be a JSON object")
    assignments = select_residue_assignments(counts)

    work_dir = Path(output_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    renamed_pdb = work_dir / "amber-renamed.pdb"
    final_pdb = work_dir / "amber-protonated.pdb"
    leap_input = work_dir / "tleap.in"
    stdout_log = work_dir / "tleap.stdout.log"
    stderr_log = work_dir / "tleap.stderr.log"
    rename_result = write_renamed_pdb(input_pdb, renamed_pdb, assignments)

    leap_input.write_text(
        "source leaprc.protein.ff19SB\n"
        "source leaprc.water.tip3p\n"
        "mol = loadpdb amber-renamed.pdb\n"
        "check mol\n"
        "savepdb mol amber-protonated.pdb\n"
        "quit\n",
        encoding="ascii",
    )
    final_pdb.unlink(missing_ok=True)
    try:
        completed = subprocess.run(
            ["tleap", "-f", leap_input.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise RuntimeError(f"Could not start TLEaP: {error}") from error
    stdout_log.write_text(completed.stdout, encoding="utf-8")
    stderr_log.write_text(completed.stderr, encoding="utf-8")

    error_summaries = re.findall(r"Errors\s*=\s*(\d+)", completed.stdout)
    warning_summaries = re.findall(r"Warnings\s*=\s*(\d+)", completed.stdout)
    leap_errors = int(error_summaries[-1]) if error_summaries else 0
    leap_warnings = int(warning_summaries[-1]) if warning_summaries else 0
    if completed.returncode or leap_errors or not final_pdb.is_file() or not final_pdb.stat().st_size:
        raise RuntimeError(
            f"TLEaP failed (exit status {completed.returncode}, reported errors "
            f"{leap_errors}); see {stdout_log} and {stderr_log}"
        )

    assignment_map = work_dir / "residue-assignments.json"
    assignment_map.write_text(
        json.dumps(rename_result["assignments"], indent=2) + "\n", encoding="utf-8"
    )
    return {
        "status": "completed",
        "renamed_input_pdb": str(renamed_pdb),
        "tleap_input": str(leap_input),
        "tleap_stdout": str(stdout_log),
        "tleap_stderr": str(stderr_log),
        "tleap_log": str(work_dir / "leap.log"),
        "assignment_map": str(assignment_map),
        "output_pdb": str(final_pdb),
        "renamed_residue_count": rename_result["renamed_residue_count"],
        "assigned_residue_count": len(assignments),
        "tleap_warnings": leap_warnings,
    }


def _parse_label(label: str) -> tuple[str, str, int, str]:
    try:
        chain, original_name, residue_label = label.rsplit("_", 2)
    except ValueError as error:
        raise ValueError(f"Invalid assignment residue key {label!r}") from error
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", residue_label)
    if not match or len(chain) > 1 or len(original_name) != 3:
        raise ValueError(f"Invalid assignment residue key {label!r}")
    return chain, original_name.upper(), int(match.group(1)), match.group(2)
