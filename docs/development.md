# Development and tests

Development-only material is excluded from a normal source distribution by `.gitignore`; generated run products belong in `test-results/` or another chosen output directory. In a development checkout that includes `tests/`, run the unit suite in the project environment:

## Repository layout

- `src/traybake/`: CLI orchestration plus focused modules for PDB parsing/cleanup, selections, geometry, capping/hydrogenation, QupKake, protonation voting, and TLEaP.
- `examples/`: redistributable example PDB and selection YAML.
- `environment.yml`: Python and scientific software environment declaration.
- `setup_environment.sh`: create the environment and install the patched QupKake checkout.
- `traybake.sh`: launch TRAYBAKE from the source tree with the local environment.
- `patches/`: compatibility patch applied to QupKake during setup.
- `docs/`: installation, input, usage, and workflow documentation.

## Run tests

```bash
conda run --prefix .conda-env python -m pytest
```

The tests exercise PDB parsing and cleanup, geometry, selection matching, capping, proton addition, QupKake result handling, pH voting, TLEaP integration through a fake executable, and CLI output. Most external processes are mocked; tests do not constitute an end-to-end scientific validation of QupKake or Amber results.

For a local environment that already has the declared test dependencies, this also works:

```bash
PYTHONPATH=src python -m pytest
```

The production launch path remains `./traybake.sh`; see the [installation guide](installation.md).
