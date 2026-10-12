"""Translate each QupKake bubble's pKa predictions into weighted votes."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Any


TITRATABLE_SIDECHAIN_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    "HIS": {"ND1", "NE2"},
}
PKA_TYPES = {"ASP": "acidic", "GLU": "acidic", "LYS": "basic", "ARG": "basic", "HIS": "basic"}
VOTE_AUDIT_FIELDS = (
    "bubble", "residue", "residue_name", "pka_type_used", "site_pkas_used",
    "pka_mean", "assignment", "status", "centered_on_residue", "vote_weight",
)


def make_protonation_votes(
    bubble_results: list[dict[str, Any]], ph: float
) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]]]:
    """Assign a pH-specific state per residue/bubble and sum weighted votes.

    Predictions at equivalent atoms are averaged first. Each residue contributes
    one vote per bubble with a relevant sidechain pKa; its own centered bubble
    has weight 3 and all neighboring bubbles have weight 1.
    """
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError(f"pH must be a finite value between 0 and 14, got {ph!r}")

    audit_rows: list[dict[str, Any]] = []
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for bubble in bubble_results:
        target = bubble["target"]
        bubble_label = _residue_label(target)
        qresult = bubble.get("qupkake") or {}
        predictions = qresult.get("predictions", [])
        by_residue: dict[tuple[Any, ...], dict[str, Any]] = {}
        for prediction in predictions:
            residue = prediction.get("residue")
            if not residue:
                continue
            key = _residue_identity(residue)
            target_data = by_residue.setdefault(key, {"residue": residue, "sites": defaultdict(list)})
            target_data["sites"][residue["atom_name"].strip()].append(prediction)

        titratable_residues = [
            residue
            for residue in bubble.get("residues", [])
            if residue.get("residue_name", "").upper() in TITRATABLE_SIDECHAIN_ATOMS
        ]
        for residue_info in titratable_residues:
            residue_name = residue_info["residue_name"].upper()
            residue_label = _residue_label(residue_info)
            identity = _residue_identity(residue_info)
            prediction_group = by_residue.get(identity)
            site_pkas: dict[str, float] = {}
            if prediction_group:
                for atom_name, atom_predictions in prediction_group["sites"].items():
                    if atom_name not in TITRATABLE_SIDECHAIN_ATOMS[residue_name]:
                        continue
                    matching = [
                        float(item["pka"])
                        for item in atom_predictions
                        if item.get("pka_type", "").lower() == PKA_TYPES[residue_name]
                    ]
                    if matching:
                        if not all(math.isfinite(value) for value in matching):
                            raise ValueError(f"Non-finite pKa prediction for {residue_label}")
                        site_pkas[atom_name] = sum(matching) / len(matching)

            centered = identity == _residue_identity(target)
            weight = 3 if centered else 1
            pka_type = PKA_TYPES[residue_name]
            if not site_pkas:
                audit_rows.append(
                    _audit_row(
                        bubble_label, residue_label, residue_name, pka_type, {}, None,
                        None, "no_relevant_sidechain_pka", centered, weight,
                    )
                )
                continue

            pka_mean = sum(site_pkas.values()) / len(site_pkas)
            state = _state_at_ph(residue_name, site_pkas, pka_mean, ph)
            counts[residue_label][state] += weight
            audit_rows.append(
                _audit_row(
                    bubble_label, residue_label, residue_name, pka_type, site_pkas,
                    pka_mean, state, "assigned", centered, weight,
                )
            )

    return audit_rows, {
        residue: dict(sorted(states.items()))
        for residue, states in sorted(counts.items())
    }


def _state_at_ph(
    residue_name: str, site_pkas: dict[str, float], pka_mean: float, ph: float
) -> str:
    if residue_name in {"ASP", "GLU"}:
        if ph < pka_mean:
            return "ASH" if residue_name == "ASP" else "GLUP"
        return residue_name
    if residue_name == "LYS":
        return "LYS" if ph < pka_mean else "LYN"
    if residue_name == "ARG":
        return "ARG" if ph < pka_mean else "AR0"
    if ph < pka_mean:
        return "HIP"

    # For a neutral histidine, the predicted N-H site determines its tautomer.
    # When both nitrogens are predicted, the lower-pKa site is the N-H site.
    nh_site = min(site_pkas, key=site_pkas.get)
    return "HID" if nh_site == "ND1" else "HIE"


def _audit_row(
    bubble: str,
    residue: str,
    residue_name: str,
    pka_type: str,
    site_pkas: dict[str, float],
    pka_mean: float | None,
    assignment: str | None,
    status: str,
    centered: bool,
    weight: int,
) -> dict[str, Any]:
    return {
        "bubble": bubble,
        "residue": residue,
        "residue_name": residue_name,
        "pka_type_used": pka_type,
        "site_pkas_used": ";".join(
            f"{atom}:{value:.6f}" for atom, value in sorted(site_pkas.items())
        ),
        "pka_mean": "" if pka_mean is None else f"{pka_mean:.6f}",
        "assignment": assignment or "",
        "status": status,
        "centered_on_residue": str(centered).lower(),
        "vote_weight": weight,
    }


def _residue_identity(residue: dict[str, Any]) -> tuple[Any, ...]:
    return (
        residue.get("model", 1),
        residue.get("chain", ""),
        residue.get("residue_number"),
        residue.get("insertion_code", ""),
        residue.get("residue_name", "").upper(),
    )


def _residue_label(residue: dict[str, Any]) -> str:
    chain = residue.get("chain", "") or "_"
    insertion = residue.get("insertion_code", "")
    number = residue.get("residue_number")
    name = residue.get("residue_name", "").upper()
    return f"{chain}_{name}_{number}{insertion}"
