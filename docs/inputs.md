# Inputs and selections

## PDB input

TRAYBAKE accepts a PDB file. The included [`examples/UBQ1.pdb`](../examples/UBQ1.pdb) is the ubiquitin structure 1UBQ and is intended as a small runnable example. Input interpretation follows standard fixed-column PDB records; separate models are treated independently.

Before bubble generation, TRAYBAKE keeps the 20 canonical amino-acid residue names and removes other residue types, including waters, ligands, cofactors, nucleic acids, and modified residues. Existing hydrogens are not used as bubble-center atoms; heavy sidechain atoms determine the target center. Alternate locations are reduced to one record per atom name, preferring a blank alternate-location field and otherwise the highest occupancy.

## Target selection YAML

The included [`examples/selections.yaml`](../examples/selections.yaml) selects ASP and GLU residues in chain A:

```yaml
- {CHAIN_ID: A, RES_NAME: [ASP, GLU], RES_ID: _, ATOM_NAME: _}
```

Each YAML document is a list of selection dictionaries. Every dictionary must include all four keys: `CHAIN_ID`, `RES_NAME`, `RES_ID`, and `ATOM_NAME`.

- `_` matches any value.
- A list matches any of its listed values.
- Fields in one dictionary are combined with AND.
- Separate dictionaries are combined with OR.
- Matching any atom selects its residue as a bubble target; neighboring residues remain in the bubble as structural context.

Without `--selections`, all ASP, GLU, LYS, ARG, and HIS residues are used as targets. Selection restricts bubble centers only; it does not remove nearby residues from the fragment.
