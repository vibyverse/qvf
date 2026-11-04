Proposed format change, awaiting maintainer acceptance (GOVERNANCE.md,
2026-11-04): `volume.density`, `volume.spin` and `volume.orbital` sections
may name the `wavefunction.gto` section they were computed from
(`wavefunction_ref`), and a `volume.orbital` its 0-based coefficient row
(`mo_index`) and, for an unrestricted wavefunction, its channel (`spin`),
spec § 5.2, *Volumes of a wavefunction*. The Python, C++, JavaScript and Rust
validators check the reference, the row and the presence of `spin` and
report violations as `E_REF_UNRESOLVED`. Python `add_volume` gains
`wavefunction_ref=`, `mo_index=` and `spin=`, C++ gains
`QvfWriter::link_volume_wavefunction`, and the PySCF grid exporter writes
the fields beside `group`. Corpus: `volume_wavefunction_ref.qvf` and five
negative archives. Optional fields only, no `qvf_version` bump. Status:
draft PR, not on main.
