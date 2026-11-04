WP-62: `qvf.adapters.pyscf.add_td` / `write_td` take
`transition_densities=True` (or a dict of grid options) to write each
excited state's spin-summed transition density as a
`volume.transition_density` on a shared grid, linked through
`excited_states_ref`, `state_index` and the state's `density_ref`. Its first
moment is PySCF's transition dipole exactly. Restricted singlets and
UHF/UKS results are supported. Triplets are refused because their
spin-summed transition density vanishes. Uses existing kinds; no schema,
registry or validator change.
