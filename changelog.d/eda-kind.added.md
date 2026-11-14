Proposed for maintainer decision (GOVERNANCE.md, 2026-11-14 record): the
canonical kind `eda` (provisional, promoted from `x_vibeqc.peda`) for energy
decomposition analyses of any scheme — Ziegler–Rauk/ETS and periodic EDA,
LMO-EDA, GKS-EDA, Kitaura–Morokuma, ALMO-EDA, SAPT, local energy
decomposition and many-body expansions. Its strict `decomposition` payload
holds the fragments (charge, multiplicity, spin flip, basis, sub-fragments),
the terms (whole-system, per-fragment, pair and n-body scope, spin, a parent
hierarchy), declared sum rules that validators check (new finding code
`E_SUM_RULE`), an optional NOCV/COVP pair table with per-k values, and k-point
sampling. `volume.difference` gains `eda_ref`, `eda_role`, `pair_index` and
`pair_spin` for deformation densities. The Python writer gains `add_eda` and
the eda fields of `add_volume`, the JavaScript writer `addEda`; the C++, C,
Fortran and Rust typed kinds are generated. The corpus gains six positive and
ten negative archives. No `qvf_version` or package version change.
