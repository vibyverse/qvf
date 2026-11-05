Close the two open items of the Quantum ESPRESSO converters (WP-48).
`espresso-filplot` imports the pp.x `filplot` intermediate file as a
structure and one grid periodic without the endpoint plane, with the kind
and unit from `plot_num`: densities in e/bohr³ as written, electrostatic
potentials (plot_num 11, 12) as −V/2 in hartree/e and the other potentials
(plot_num 1, 2) as V/2 in hartree. `espresso-in` now expands `space_group`
input with `ATOMIC_POSITIONS crystal_sg` sites (`uniqueb`, `origin_choice`,
`rhombohedral`) to the conventional unit cell through the Hall symbols of
all 230 space groups. Wyckoff labels in `crystal_sg` lines remain open and
are refused. The tests compare with the CIF import, `ase.spacegroup.crystal`,
spglib and the espresso-cube import of the same data. No format, schema or
validator change.
