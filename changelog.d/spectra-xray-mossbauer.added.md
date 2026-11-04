Proposed, awaiting maintainer acceptance (GOVERNANCE.md, 2026-11-04 record):
three canonical kinds for format review P13. `spectra.xas` and `spectra.xes`
are one-dimensional spectra on the shared `SpectrumPayload` (x in eV), which
gains a strict optional `xray` block (`edge`, `element`, `atom_indices`,
`core_orbital`, `line`, `intensity`, `energy_shift_ev`) and the peak
`assignment.atom_index`. `spectra.mossbauer` carries per-site isomer shift,
quadrupole splitting, asymmetry parameter, EFG tensor, contact density,
hyperfine field, linewidth and population in a strict payload, 57Fe by
default. The Python writer gains `add_xray_spectrum` and `add_mossbauer`,
the C++ library the generated typed kinds `SpectraXAS`, `SpectraXES` and
`SpectraMossbauer`. The corpus gains three positive and nine negative
archives. No `qvf_version` or package version change.
