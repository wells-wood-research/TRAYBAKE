# TRAYBAKE

TRAYBAKE prepares local protein fragments around titratable residues, runs QupKake to predict protonation sites, aggregates pH-dependent assignments, and uses AmberTools/TLEaP to write a protonated PDB.

## Quick start

```bash
./setup_environment.sh
./traybake.sh examples/UBQ1.pdb --output-dir test-results/ubq1 \
  --selections examples/selections.yaml
```

Setup creates the project-local Conda environment and installs QupKake. The example command writes a report, intermediate fragments, audit files, and a protonated PDB under the output directory. See the [installation guide](docs/installation.md) and [usage guide](docs/usage.md) before processing your own structures.

## Documentation

- [Install and set up dependencies](docs/installation.md)
- [Run TRAYBAKE and inspect outputs](docs/usage.md)
- [Input files and selection syntax](docs/inputs.md)
- [Workflow, assumptions, and limitations](docs/workflow.md)
- [Development and tests](docs/development.md)

The checked-in example inputs are [`examples/UBQ1.pdb`](examples/UBQ1.pdb) and [`examples/selections.yaml`](examples/selections.yaml). Generated results and development-only directories are not part of the distributed source tree.

**License:** TRAYBAKE is distributed under the [MIT License](LICENSE). Review the terms for QupKake and the scientific dependencies separately.
