# Installation

TRAYBAKE is currently run from a source checkout; it is not published as a pip-installable package. Its Python and scientific dependencies are declared in [`environment.yml`](../environment.yml). Setup also obtains QupKake and installs a local compatibility patch.

## Requirements

- Linux or macOS with Conda (or a compatible `conda` command) and Git.
- Network access for the first setup, to obtain Conda packages and the QupKake source repository.
- Enough disk space for the scientific Python stack and QupKake/xTB. The environment can be several gigabytes.

AmberTools, including `tleap`, is installed as part of the declared Conda environment. TRAYBAKE runs TLEaP automatically at the end of a normal workflow.

## Create the environment

From the repository root:

```bash
./setup_environment.sh
```

The script creates `.conda-env/` from `environment.yml`, checks out QupKake at the pinned upstream commit `f7a294e91929e4e4f0a7100ad5cf5a64a34a7d4f` in the ignored `dependencies/QupKake/` directory, applies `patches/qupkake-traybake-compat.patch`, and installs QupKake into the project environment. If a QupKake checkout already exists at a different commit, setup stops without changing it. It also stops if `.conda-env/` already exists; setup will not replace an existing environment.

Run the wrapper without activating Conda:

```bash
./traybake.sh --help
```

Or activate the environment and run the Python module directly:

```bash
conda activate ./.conda-env
PYTHONPATH=src python -m traybake --help
```

The wrapper sets the source path and limits numerical-library threads for you. Use `./traybake.sh` for normal runs.

## Reproducibility

The environment file specifies Python 3.9 and the direct Conda/pip dependencies and version constraints. Conda may resolve newer compatible builds over time where versions are not pinned; the QupKake source revision is pinned, and TRAYBAKE's compatibility changes are maintained as a patch rather than distributing a modified QupKake checkout. The patch covers the Wiberg bond-order fallback, preserves existing 3D coordinates, disables QupKake's xTB geometry optimization for these prepared fragments, and improves xTB limits, diagnostics, and failure handling. For precise reproducibility, record the Conda platform and resolved package list alongside each scientific run.
