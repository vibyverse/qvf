Accepted by the maintainer on 2026-11-05 (#20, GOVERNANCE.md record of
2026-11-04 on the data unit of a volume grid): a new error `E_UNIT_CONTRADICTED` in spec Appendix B. The Python, C++,
JavaScript and Rust validators report a grid `units.data` that is a unit name
other than the one spec § 4.3 fixes for the section's kind: `e_per_bohr3`
for `volume.density` and `volume.transition_density`, `hartree_per_e` for
`volume.potential`, `dimensionless` for `volume.rdg` and `e_per_bohr5` for
`volume.laplacian`. The check is semantic, so it runs with or without the
JSON Schema check. New negative archive `unit_grid_data_contradicts_kind.qvf`.
No `qvf_version`, schema or package version change.
