# Workflow and scientific assumptions

1. **Read and clean the PDB.** Retain canonical amino acids and remove all other residue types. Model and residue identifiers (model, chain, number, insertion code, and name) are carried through the workflow.
2. **Choose titration targets.** By default, use ASP, GLU, LYS, ARG, and HIS; an optional YAML selection further restricts targets.
3. **Build local bubbles.** The target center is the unweighted centroid of its observed heavy sidechain atoms. A residue is included when any observed atom is within or on the requested radius (5 Å by default).
4. **Cap and hydrogenate fragments.** Exposed peptide N/C atoms receive caps. A residue is retained whole only when N, CA, C, and O are present and in range; otherwise CB and sidechain atoms are retained and a hydrogen cap is placed at the sidechain cut. Hydrogens are added from residue templates while keeping observed heavy-atom coordinates fixed.
5. **Run QupKake.** Each fragment is processed in its own directory. Predictions are mapped back to original residue identities. Successful runs can be resumed from checkpoints; failures and tool logs are retained.
6. **Aggregate predicted pKa values.** Values are classified at the requested pH. A target's own centered bubble contributes weight 3 and other bubbles weight 1. Histidine's neutral tautomer is chosen from the predicted N-H site. The CSV and JSON assignment files preserve this tally.
7. **Prepare the output with TLEaP.** The cleaned structure is renamed according to the winning weighted assignments and passed to TLEaP with Amber ff19SB/TIP3P. GLUP is translated to Amber's GLH residue name. The original input is untouched.

## Important limitations

- Noncanonical residues, waters, ligands, and cofactors are removed. Such species can influence protonation, so this cleanup may not be appropriate for every structure.
- Fragment capping and generated hydrogen geometry should be chemically validated for the input structure.
- The weighted center-bubble consensus and pKa-to-state rules are implementation choices, not calibrated probabilities. Review the per-bubble audit and validate assignments for the scientific system.
- TLEaP warnings should be reviewed. The workflow does not add solvent or ions, and the final PDB should be validated before simulation.
- QupKake's applicability to capped peptide fragments and the resulting protonation predictions should be independently evaluated.
