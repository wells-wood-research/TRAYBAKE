"""Command-line dry run for selecting titratable-residue neighborhoods."""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import re
import shutil
from pathlib import Path

from traybake import __version__
from traybake.capping import build_capped_bubble
from traybake.cleanup import remove_nonprotein, write_cleaned_pdb
from traybake.geometry import build_bubbles
from traybake.pdb import read_pdb, write_atoms
from traybake.protonation import VOTE_AUDIT_FIELDS, make_protonation_votes
from traybake.protons import add_residue_protons
from traybake.qupkake import _PROGRESS_COLOURS, prepare_job, run_jobs
from traybake.selections import load_selections, selected_residues
from traybake.splash import cmap as SPLASH_CMAP, splashText as SPLASH_TEXT
from traybake.tleap import run_tleap


def main() -> None:
    _print_splash()
    parser = argparse.ArgumentParser(description="Plan 5 Å local PDB neighborhoods for titratable residues.")
    parser.add_argument("pdb", type=Path, help="input PDB file")
    parser.add_argument("--radius", type=float, default=5.0, help="bubble radius in Å (default: 5.0)")
    parser.add_argument(
        "--selections",
        type=Path,
        help="YAML list of drMD-style custom selections for titratable target residues",
    )
    parser.add_argument(
        "--ph",
        "--pH",
        type=float,
        default=7.4,
        dest="ph",
        help="pH used to translate predicted pKa values to residue states (default: 7.4)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="directory for per-target fragments and automatic Amber preparation",
    )
    parser.add_argument(
        "--qupkake-executable",
        default="qupkake",
        help="QupKake CLI executable (default: qupkake)",
    )
    parser.add_argument(
        "--n-workers",
        "--n_workers",
        type=int,
        default=min(4, mp.cpu_count()),
        dest="n_workers",
        help="maximum QupKake bubbles run concurrently (default: up to 4)",
    )
    parser.add_argument(
        "--clean-up",
        type=int,
        choices=(0, 1),
        default=0,
        help="remove amber_prep/, bubbles/, and amber-prep/ after preparation (default: 0)",
    )
    args = parser.parse_args()
    if not args.output_dir:
        parser.error("--output-dir is required because Amber preparation runs automatically")
    if args.n_workers < 1:
        parser.error("--n-workers must be at least 1")
    if not 0 <= args.ph <= 14:
        parser.error("--ph must be between 0 and 14")

    original_structure = read_pdb(args.pdb)
    cleanup = remove_nonprotein(original_structure)
    structure = cleanup.structure
    selections = None
    matched_residues = None
    if args.selections:
        try:
            selections = load_selections(args.selections)
        except ValueError as error:
            parser.error(str(error))
        matched_residues = selected_residues(structure, selections)
    bubbles = build_bubbles(structure, args.radius, matched_residues)
    bubble_dir = args.output_dir / "bubbles" if args.output_dir else None
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        bubble_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "qupkake").mkdir(parents=True, exist_ok=True)

    result = {
        "traybake_version": __version__,
        "input": str(args.pdb),
        "cleanup": {
            "rule": "retain canonical amino-acid residues; remove all other residue types",
            "removed_residue_count": len(cleanup.removed_residues),
            "removed_residues": [
                {
                    "model": key.model,
                    "chain": key.chain,
                    "residue_number": key.number,
                    "insertion_code": key.insertion_code,
                    "residue_name": key.name,
                }
                for key in cleanup.removed_residues
            ],
        },
        "radius_angstrom": args.radius,
        "ph": args.ph,
        "selection": None if args.selections is None else {
            "file": str(args.selections),
            "matching_residue_count": len(matched_residues),
            "target_residue_count": len(bubbles),
            "target_residues": [
                {
                    "model": bubble.center.model,
                    "chain": bubble.center.chain,
                    "residue_number": bubble.center.number,
                    "insertion_code": bubble.center.insertion_code,
                    "residue_name": bubble.center.name,
                }
                for bubble in bubbles
            ],
            "applies_to": "titratable bubble targets only; nearby residues remain context",
        },
        "center_definition": "unweighted centroid of observed heavy sidechain atoms",
        "membership_rule": "include a residue if any observed atom is within or on the radius",
        "targets": [],
    }
    qupkake_jobs = []
    qupkake_job_targets = []
    for bubble in bubbles:
        capped = build_capped_bubble(structure, bubble, args.radius)
        fragment_atoms = add_residue_protons(capped.atoms)
        entry = {
            "target": {
                "model": bubble.center.model,
                "chain": bubble.center.chain,
                "residue_number": bubble.center.number,
                "insertion_code": bubble.center.insertion_code,
                "residue_name": bubble.center.name,
            },
            "center_angstrom": bubble.center_of_mass,
            "residue_count": len(bubble.residues),
            "fragment_residue_count": sum(mode != "omitted_no_cb" for _, mode in capped.residue_modes),
            "fragment_atom_count": len(fragment_atoms),
            "residue_modes": [
                {
                    "model": key.model,
                    "chain": key.chain,
                    "residue_number": key.number,
                    "insertion_code": key.insertion_code,
                    "residue_name": key.name,
                    "mode": mode,
                }
                for key, mode in capped.residue_modes
            ],
            "caps": [
                {
                    "model": cap.residue.model,
                    "chain": cap.residue.chain,
                    "residue_number": cap.residue.number,
                    "insertion_code": cap.residue.insertion_code,
                    "residue_name": cap.residue.name,
                    "atom_name": cap.atom_name,
                    "attached_to": cap.attached_to,
                    "kind": cap.kind,
                    "coordinates_angstrom": cap.coordinates,
                }
                for cap in capped.caps
            ],
            "residues": [
                {
                    "model": key.model,
                    "chain": key.chain,
                    "residue_number": key.number,
                    "insertion_code": key.insertion_code,
                    "residue_name": key.name,
                }
                for key in bubble.residues
            ],
        }
        if args.output_dir:
            fragment = bubble_dir / f"{bubble.center.label()}.pdb"
            write_atoms(fragment, fragment_atoms)
            entry["fragment"] = str(fragment)
            try:
                job = prepare_job(
                    fragment,
                    args.output_dir / "qupkake" / bubble.center.label(),
                    bubble.center,
                    args.qupkake_executable,
                )
            except (OSError, RuntimeError, ValueError) as error:
                entry["qupkake"] = {"status": "failed", "error": str(error)}
            else:
                qupkake_jobs.append(job)
                qupkake_job_targets.append(entry)
        result["targets"].append(entry)

    try:
        qupkake_results = run_jobs(qupkake_jobs, args.n_workers)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
    for entry, qupkake_result in zip(qupkake_job_targets, qupkake_results):
        entry["qupkake"] = qupkake_result
    result["qupkake"] = {
        "executable": args.qupkake_executable,
        "bubble_workers": min(args.n_workers, len(qupkake_jobs)) if qupkake_jobs else 0,
        "cpus_per_run": 1,
    }

    ph_label = format(args.ph, ".15g")
    try:
        audit_rows, assignment_counts = make_protonation_votes(result["targets"], args.ph)
    except ValueError as error:
        parser.error(str(error))
    votes_path = args.output_dir / "protonation-assignment-votes.csv"
    assignments_path = args.output_dir / f"protonation-assignment-counts-pH{ph_label}.json"
    with votes_path.open("w", newline="", encoding="utf-8") as votes_file:
        writer = csv.DictWriter(votes_file, fieldnames=VOTE_AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(audit_rows)
    assignments_path.write_text(
        json.dumps(assignment_counts, indent=2) + "\n", encoding="utf-8"
    )
    result["protonation_votes"] = {
        "ph": args.ph,
        "weighting": "centered bubble x3; other bubbles x1",
        "audit_csv": str(votes_path),
        "weighted_counts_json": str(assignments_path),
    }
    amber_dir = args.output_dir / "amber-prep"
    cleaned_pdb = write_cleaned_pdb(
        args.pdb, amber_dir / "cleaned.pdb", cleanup
    )
    try:
        result["amber"] = run_tleap(
            input_pdb=cleaned_pdb,
            counts_path=assignments_path,
            output_dir=amber_dir,
        )
        final_pdb = Path(result["amber"]["output_pdb"])
        published_pdb = args.output_dir / f"{args.pdb.stem}_protonated.pdb"
        if final_pdb.resolve() != published_pdb.resolve():
            shutil.copy2(final_pdb, published_pdb)
        result["amber"]["output_pdb"] = str(published_pdb.resolve())
        result["amber"]["assignment_counts"] = str(assignments_path.resolve())
        result["amber"]["ph"] = args.ph
        if args.clean_up:
            removed_directories = []
            for directory_name in ("amber_prep", "bubbles", "amber-prep"):
                directory = args.output_dir / directory_name
                if directory.is_dir() and not directory.is_symlink():
                    shutil.rmtree(directory)
                    removed_directories.append(directory_name)
            for key in (
                "renamed_input_pdb",
                "tleap_input",
                "tleap_stdout",
                "tleap_stderr",
                "tleap_log",
                "assignment_map",
                "cleaned_input_pdb",
            ):
                result["amber"].pop(key, None)
            result["amber"]["removed_intermediate_directories"] = removed_directories
        else:
            result["amber"]["cleaned_input_pdb"] = str(cleaned_pdb.resolve())
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
    report_path = args.output_dir / "traybake-report.json"
    report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def _print_splash() -> None:
    print(_color_splash(SPLASH_TEXT, SPLASH_CMAP), end="")


def _color_splash(splash: str, color_map: list[int]) -> str:
    output = []
    for line in splash.splitlines(keepends=True):
        markers = list(re.finditer(r"\d+", line))
        if not markers:
            output.append(line)
            continue

        output.append(line[:markers[0].start()])
        for index, marker in enumerate(markers):
            end = markers[index + 1].start() if index + 1 < len(markers) else len(line)
            color_index = color_map[int(marker.group())]
            output.append(_color_section(line[marker.end():end], color_index))
    return "".join(output)


def _color_section(section: str, palette_index: int) -> str:
    colour = _PROGRESS_COLOURS[palette_index % len(_PROGRESS_COLOURS)]
    red, green, blue = (int(colour[index:index + 2], 16) for index in (1, 3, 5))
    return "".join(
        character
        if character.isspace()
        else f"\033[38;2;{red};{green};{blue}m{character}\033[0m"
        for character in section
    )


if __name__ == "__main__":
    main()
