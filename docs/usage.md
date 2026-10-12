# Usage

## Basic run

```bash
./traybake.sh examples/UBQ1.pdb --output-dir test-results/ubq1
```

`--output-dir` is required. By default, all supported titratable residues are considered as bubble centers. To use the included selection (ASP/GLU in chain A):

```bash
./traybake.sh examples/UBQ1.pdb \
  --output-dir test-results/ubq1-selected \
  --selections examples/selections.yaml
```

Common options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--radius ANGSTROM` | `5.0` | Bubble radius around each target sidechain center. |
| `--selections FILE` | all titratable residues | YAML target filters; nearby residues still provide context. |
| `--ph VALUE` (or `--pH`) | `7.4` | pH used to translate predicted pKa values into protonation assignments; allowed range is 0–14. |
| `--n-workers COUNT` | up to `4` | Maximum number of QupKake bubbles run concurrently. Each QupKake process uses one CPU. |
| `--qupkake-executable NAME` | `qupkake` | QupKake CLI executable to invoke. |
| `--clean-up 1` | `0` | Remove intermediate `amber_prep/`, `bubbles/`, and `amber-prep/` directories after successful Amber preparation. |

## Output files

The output directory contains:

- `traybake-report.json`: input metadata, selected targets, fragment and QupKake results, pH, and Amber preparation status.
- `<input-stem>_protonated.pdb`: final PDB produced by TLEaP (also retained inside `amber-prep/` by default).
- `protonation-assignment-votes.csv`: per-bubble pKa-to-state audit.
- `protonation-assignment-counts-pH<VALUE>.json`: weighted state counts used for residue renaming.
- `bubbles/`: capped and hydrogenated PDB fragments.
- `qupkake/`: QupKake inputs, predictions, checkpoints, and per-run stdout/stderr logs.
- `amber-prep/`: cleaned/renamed inputs, TLEaP script and logs, and assignment map.

Successful QupKake results are checkpointed and reused when the corresponding input has not changed. The final file is written only after the AmberTools/TLEaP stage succeeds. TLEaP prepares the protein with ff19SB/TIP3P parameters; it does not solvate the structure or add ions. Review warnings and validate the output before downstream use.

## Example with a different radius and pH

```bash
./traybake.sh input.pdb \
  --output-dir results/run-01 \
  --radius 6.0 \
  --ph 6.5 \
  --n-workers 2
```

The original input PDB is not edited. Paths are interpreted relative to the current working directory.
