"""Run QupKake independently on capped PDB bubbles."""

from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tqdm import tqdm

from traybake.pdb import ResidueKey

_PROGRESS_COLOURS = (
    "#ffb3ba", "#ffdfba", "#ffffba", "#baffc9",
    "#baffff", "#bae1ff", "#bdb2ff", "#ffc6ff",
    "#f4acb7", "#f7d6bf", "#f9f7a5", "#c7f0bd",
    "#a8e6cf", "#a2d2ff", "#cdb4db", "#f1c0e8",
    "#f6bdc0", "#f8d5a3", "#e9edc9", "#d8f3dc",
    "#caffbf", "#9bf6ff", "#cde7f0", "#d7c0f1",
    "#fbc4ab", "#ffe5ec", "#fce1a8", "#d0f4de",
    "#bde0fe", "#cddafd", "#e4c1f9", "#f2c6de",
)


@dataclass(frozen=True)
class QupKakeJob:
    executable: str
    input_sdf: Path
    root: Path
    center: dict[str, Any]
    atom_residues: dict[int, dict[str, Any]]
    input_pdb: Path | None = None

    @property
    def output_name(self) -> str:
        return f"{self.input_sdf.stem}-predictions.sdf"


def prepare_job(
    fragment: str | Path,
    root: str | Path,
    center: ResidueKey,
    executable: str,
) -> QupKakeJob:
    """Describe a capped PDB job, deferring its prediction files until start."""
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError(
            "RDKit is required to convert PDB bubbles for QupKake; run "
            "./setup_environment.sh to create the project environment."
        ) from error

    fragment = Path(fragment).resolve()
    root = Path(root).resolve()
    molecule = Chem.MolFromPDBFile(
        str(fragment), sanitize=True, removeHs=False, proximityBonding=True
    )
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError(f"RDKit could not read a molecule from {fragment}")

    center_dict = _residue_dict(center)
    molecule.SetProp("_Name", center.label())
    atom_residues: dict[int, dict[str, Any]] = {}
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        key = ResidueKey(
            model=center.model,
            chain=info.GetChainId().strip(),
            number=info.GetResidueNumber(),
            insertion_code=info.GetInsertionCode().strip(),
            name=info.GetResidueName().strip().upper(),
        )
        atom_residues[atom.GetIdx()] = {
            **_residue_dict(key),
            "atom_name": info.GetName().strip(),
            "element": atom.GetSymbol(),
            "is_center": key == center,
        }

    input_sdf = root / "input.sdf"
    return QupKakeJob(
        executable=executable,
        input_sdf=input_sdf,
        root=root / "data",
        center=center_dict,
        atom_residues=atom_residues,
        input_pdb=fragment,
    )


def run_jobs(jobs: list[QupKakeJob], workers: int) -> list[dict[str, Any]]:
    """Run jobs in parallel, showing one protonation-option bar per bubble."""
    if not jobs:
        return []
    if workers < 1:
        raise ValueError("workers must be at least 1")

    executable = shutil.which(jobs[0].executable)
    if executable is None:
        raise FileNotFoundError(
            f"QupKake executable {jobs[0].executable!r} was not found on PATH"
        )
    normalized = [
        QupKakeJob(
            executable=executable,
            input_sdf=job.input_sdf,
            root=job.root,
            center=job.center,
            atom_residues=job.atom_residues,
            input_pdb=job.input_pdb,
        )
        for job in jobs
    ]
    results: list[dict[str, Any] | None] = [None] * len(normalized)
    worker_count = min(workers, len(normalized))
    bars = [
        tqdm(
            desc=_progress_label(job.center, index),
            position=index,
            unit="option",
            colour=_progress_colour(index),
            leave=True,
        )
        for index, job in enumerate(normalized)
    ]
    try:
        with mp.Pool(processes=worker_count) as pool:
            if worker_count == 1:
                # Start jobs one at a time so a failure does not create files for
                # later bubbles, while polling the active job's option artifacts.
                for index, job in enumerate(normalized):
                    pending = pool.apply_async(_run_job, (job,))
                    while not pending.ready():
                        _refresh_progress(job, bars[index])
                        time.sleep(0.1)
                    result = pending.get()
                    _finish_progress(job, bars[index], result)
                    if result["status"] == "failed":
                        pool.terminate()
                        _raise_job_failure(result)
                    results[index] = result
            else:
                pending = {
                    index: pool.apply_async(_run_indexed_job, ((index, job),))
                    for index, job in enumerate(normalized)
                }
                while pending:
                    for index, async_result in list(pending.items()):
                        job = normalized[index]
                        _refresh_progress(job, bars[index])
                        if not async_result.ready():
                            continue
                        result_index, result = async_result.get()
                        _finish_progress(job, bars[index], result)
                        if result["status"] == "failed":
                            pool.terminate()
                            _raise_job_failure(result)
                        results[result_index] = result
                        del pending[index]
                    if pending:
                        time.sleep(0.1)
    finally:
        for bar in bars:
            bar.close()
    return [result for result in results if result is not None]


def _refresh_progress(job: QupKakeJob, bar: Any) -> None:
    """Update a bubble bar from QupKake's completed protonation-pair files."""
    completed = _completed_options(job)
    if completed > bar.n:
        bar.update(completed - bar.n)

    # QupKake writes this SDF after predicting sites and before processing its
    # protonation/deprotonation pairs. Its records are the exact bar total.
    options_path = job.root / "raw" / job.output_name
    if options_path.is_file():
        total = max(bar.n, _count_sdf_records(options_path))
        if total != bar.total:
            bar.total = total
            bar.refresh()


def _finish_progress(job: QupKakeJob, bar: Any, result: dict[str, Any]) -> None:
    _refresh_progress(job, bar)
    if result.get("status") == "failed":
        bar.set_postfix_str("failed")
        bar.refresh()
        return
    # QupKake removes the temporary option SDF at the end of a successful run.
    # In that case, every option it attempted is represented by a processed
    # pair file or a recorded xTB failure. Cached runs can instead use outputs.
    total = max(bar.n, len(result.get("predictions", [])))
    if bar.total is None or bar.total < total:
        bar.total = total
    if bar.n < bar.total:
        bar.update(bar.total - bar.n)
    bar.set_postfix_str(result.get("status", "completed"))
    bar.refresh()


def _completed_options(job: QupKakeJob) -> int:
    processed_dir = job.root / "processed"
    completed = sum(1 for _ in processed_dir.glob("*.pt")) if processed_dir.is_dir() else 0
    error_log = job.root / "logs" / "error_log.txt"
    if error_log.is_file():
        try:
            log_text = error_log.read_text(encoding="utf-8", errors="replace")
            completed += log_text.count("Skipping protonation pair row ")
            completed += log_text.count("Error processing row ")
        except OSError:
            pass
    return completed


def _count_sdf_records(path: Path) -> int:
    try:
        return path.read_text(encoding="utf-8", errors="replace").count("$$$$")
    except OSError:
        return 0


def _progress_label(center: dict[str, Any], index: int) -> str:
    try:
        return _center_label(center)
    except (KeyError, TypeError):
        residue = center.get("residue_name", "residue")
        number = center.get("residue_number", index + 1)
        return f"bubble {index + 1}: {residue} {number}"


def _progress_colour(index: int) -> str:
    """Return a pastel colour, cycling through the 32-colour palette."""
    return _PROGRESS_COLOURS[index % len(_PROGRESS_COLOURS)]


def _run_indexed_job(
    indexed_job: tuple[int, QupKakeJob],
) -> tuple[int, dict[str, Any]]:
    index, job = indexed_job
    return index, _run_job(job)


def _raise_job_failure(result: dict[str, Any]) -> None:
    center = result.get("center", {})
    target = " ".join(
        str(value)
        for value in (
            center.get("residue_name"),
            center.get("chain"),
            center.get("residue_number"),
        )
        if value is not None
    )
    raise RuntimeError(
        f"QupKake failed for {target}: {result.get('error', 'unknown error')}"
    )


def _run_job(job: QupKakeJob) -> dict[str, Any]:
    job.root.mkdir(parents=True, exist_ok=True)
    if job.input_pdb is not None:
        try:
            from rdkit import Chem

            molecule = Chem.MolFromPDBFile(
                str(job.input_pdb), sanitize=True, removeHs=False, proximityBonding=True
            )
            if molecule is None or molecule.GetNumAtoms() == 0:
                raise ValueError(f"RDKit could not read a molecule from {job.input_pdb}")
            molecule.SetProp("_Name", _center_label(job.center))
            job.input_sdf.parent.mkdir(parents=True, exist_ok=True)
            writer = Chem.SDWriter(str(job.input_sdf))
            try:
                writer.write(molecule)
            finally:
                writer.close()
        except (ImportError, OSError, ValueError, RuntimeError) as error:
            return _failed(job, f"Could not prepare QupKake input: {error}")

    output_path = job.root / "output" / job.output_name
    checkpoint_path = job.root / "traybake.checkpoint.json"
    input_digest = hashlib.sha256(job.input_sdf.read_bytes()).hexdigest()
    cached_result = _load_checkpoint(
        job, output_path, checkpoint_path, input_digest
    )
    if cached_result is not None:
        return cached_result

    # A non-checkpointed/invalid output may be a partial file from an interrupted
    # run. Remove only this generated prediction artifact before recalculating.
    output_path.unlink(missing_ok=True)
    command = [
        job.executable,
        "file",
        str(job.input_sdf),
        "--root",
        str(job.root),
        "--multiprocessing",
        "1",
        "--output",
        job.output_name,
    ]
    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    environment["OMP_THREAD_LIMIT"] = "1"
    environment["QUPKAKE_XTB_LOG_DIR"] = str(job.root)

    try:
        completed = subprocess.run(
            command,
            cwd=job.root,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as error:
        return _failed(job, str(error))

    (job.root / "traybake.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (job.root / "traybake.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        return _failed(
            job,
            f"QupKake exited with status {completed.returncode}; see "
            f"{job.root / 'traybake.stderr.log'}",
        )

    if not output_path.is_file():
        result = {
            "status": "no_sites",
            "center": job.center,
            "predictions": [],
            "output": None,
            "stdout_log": str(job.root / "traybake.stdout.log"),
            "stderr_log": str(job.root / "traybake.stderr.log"),
        }
        _write_checkpoint(checkpoint_path, input_digest, result)
        return result
    try:
        predictions = _read_predictions(output_path, job.atom_residues)
    except (OSError, ValueError) as error:
        return _failed(job, f"Could not parse QupKake output: {error}")

    result = {
        "status": "completed",
        "center": job.center,
        "predictions": predictions,
        "output": str(output_path),
        "stdout_log": str(job.root / "traybake.stdout.log"),
        "stderr_log": str(job.root / "traybake.stderr.log"),
    }
    _write_checkpoint(checkpoint_path, input_digest, result)
    return result


def _load_checkpoint(
    job: QupKakeJob, output_path: Path, checkpoint_path: Path, input_digest: str
) -> dict[str, Any] | None:
    """Reuse a complete result when its generated input and output artifacts agree."""
    checkpoint_exists = checkpoint_path.exists()
    if checkpoint_exists:
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            status = checkpoint.get("status")
            if checkpoint.get("input_sha256") == input_digest:
                if status == "completed" and output_path.is_file():
                    predictions = _read_predictions(output_path, job.atom_residues)
                    if predictions:
                        return _cached_result(job, output_path, predictions)
                if (
                    status == "no_sites"
                    and not output_path.exists()
                    and (job.root / "traybake.stdout.log").is_file()
                    and (job.root / "traybake.stderr.log").is_file()
                ):
                    return {
                        "status": "no_sites",
                        "center": job.center,
                        "predictions": [],
                        "output": None,
                        "stdout_log": str(job.root / "traybake.stdout.log"),
                        "stderr_log": str(job.root / "traybake.stderr.log"),
                        "skipped": True,
                        "checkpoint": str(checkpoint_path),
                    }
        except (OSError, ValueError, TypeError, AttributeError):
            # A truncated checkpoint is not evidence of a completed calculation.
            pass
        # Do not let an old output SDF bypass a stale or invalid checkpoint.
        return None

    # Support outputs from earlier TRAYBAKE runs, before checkpoint files existed.
    # A prediction SDF is accepted only when it contains at least one valid record.
    if output_path.is_file():
        try:
            predictions = _read_predictions(output_path, job.atom_residues)
        except (OSError, ValueError):
            return None
        if predictions:
            result = _cached_result(job, output_path, predictions)
            _write_checkpoint(checkpoint_path, input_digest, result)
            return result
    return None


def _cached_result(
    job: QupKakeJob, output_path: Path, predictions: list[dict[str, Any]]
) -> dict[str, Any]:
    return {
        "status": "completed",
        "center": job.center,
        "predictions": predictions,
        "output": str(output_path),
        "stdout_log": str(job.root / "traybake.stdout.log"),
        "stderr_log": str(job.root / "traybake.stderr.log"),
        "skipped": True,
        "checkpoint": str(job.root / "traybake.checkpoint.json"),
    }


def _write_checkpoint(
    checkpoint_path: Path, input_digest: str, result: dict[str, Any]
) -> None:
    checkpoint = {
        "input_sha256": input_digest,
        "status": result["status"],
    }
    temporary_path = checkpoint_path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(checkpoint, sort_keys=True) + "\n", encoding="utf-8")
    temporary_path.replace(checkpoint_path)


def _read_predictions(
    path: Path, atom_residues: dict[int, dict[str, Any]]
) -> list[dict[str, Any]]:
    """Read QupKake SDF records and map reactive atom indexes to PDB residues."""
    predictions = []
    for block in path.read_text(encoding="utf-8", errors="replace").split("$$$$"):
        properties = _sdf_properties(block)
        if not properties:
            continue
        if not {"idx", "pka_type", "pka"}.issubset(properties):
            raise ValueError("prediction record is missing idx, pka_type, or pka")
        try:
            atom_index = int(properties["idx"])
            pka_text = properties["pka"]
            tensor_match = re.fullmatch(
                r"tensor\(\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
                r"(?:\s*,[^)]*)?\)",
                pka_text,
            )
            pka = float(tensor_match.group(1) if tensor_match else pka_text)
        except ValueError as error:
            raise ValueError(f"invalid prediction fields: {properties}") from error
        if not math.isfinite(pka):
            raise ValueError(f"non-finite pKa value: {properties['pka']}")
        predictions.append(
            {
                "atom_index": atom_index,
                "pka_type": properties["pka_type"],
                "pka": pka,
                "residue": atom_residues.get(atom_index),
            }
        )
    return predictions


_PROPERTY_LINE = re.compile(r"^>\s*<([^>]+)>")


def _sdf_properties(block: str) -> dict[str, str]:
    lines = block.splitlines()
    properties: dict[str, str] = {}
    index = 0
    while index < len(lines):
        match = _PROPERTY_LINE.match(lines[index])
        if match:
            properties[match.group(1)] = lines[index + 1].strip() if index + 1 < len(lines) else ""
            index += 2
        else:
            index += 1
    return properties


def _failed(job: QupKakeJob, error: str) -> dict[str, Any]:
    return {
        "status": "failed",
        "center": job.center,
        "predictions": [],
        "output": None,
        "error": error,
        "stdout_log": str(job.root / "traybake.stdout.log"),
        "stderr_log": str(job.root / "traybake.stderr.log"),
    }


def _residue_dict(key: ResidueKey) -> dict[str, Any]:
    return {
        "model": key.model,
        "chain": key.chain,
        "residue_number": key.number,
        "insertion_code": key.insertion_code,
        "residue_name": key.name,
    }


def _center_label(center: dict[str, Any]) -> str:
    chain = center.get("chain") or "_"
    insertion_code = center.get("insertion_code", "").strip()
    return (
        f"model{center['model']}_{chain}_{center['residue_name']}_"
        f"{center['residue_number']}{insertion_code}"
    )
