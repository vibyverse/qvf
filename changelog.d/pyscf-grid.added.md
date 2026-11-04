WP-62: `qvf.adapters.pyscf.add_density_grid` and `add_orbital_grid` write
the SCF electron density (`volume.density`), the spin density of open-shell
references (`volume.spin`) and selected orbitals (`volume.orbital`) on a
regular grid. Values are evaluated with PySCF's `numint` in blocks of at
most 64 MiB of AO values. `molecular_grid` builds the default box (4 bohr
margin, 0.2 bohr steps) or, for a given shape, the `tools.cubegen` grid.
Volumes name their structure at section and grid level and carry the
wavefunction id as `group`. Tests compare with the Appendix A evaluator,
`cubegen` and the electron count, with measured tolerances that depend on the grid,
as given in the PySCF adapter guide. Status: on main. Remaining in WP-62 are
triplet NTOs, unrestricted response and periodic Γ-point calculations.
Existing kinds only. No schema, registry or validator change.
