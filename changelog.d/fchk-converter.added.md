Add a Gaussian formatted checkpoint (fchk) converter: `import_fchk` and
`export_fchk` in `qvf.convert`, and `qvf convert --from fchk` / `--to fchk`.
Import reads the geometry, charge and multiplicity, the basis (s to l = 7,
pure and Cartesian, SP shells split), restricted, restricted open and
unrestricted orbitals with occupations derived from the electron counts
(checked against the SCF density), the SCF densities, the total energy and
Mulliken charges, maps Gaussian's AO order and per component Cartesian
normalization onto Appendix A, and writes the AO overlap computed from the
basis so the validator checks it. Export writes these fields back with a
loss report. Synthetic fixtures with their generation script, cclib and
IOData readings of PySCF calculations, and PySCF orbital values check the
conventions (WP-48, partial).
