The `energy_decomposition` term vocabulary gains `kinetic`, `coulomb` and
`residual` (the kinetic and Coulomb parts of Pauli and orbital terms in ADF and
BAND, and ADF's printed remainder). The specification says what `interaction`
means for a local energy decomposition without monomer calculations. New
informative page `docs/energy_decomposition_producers.md`: how ADF, BAND, ORCA
(EDA-NOCV, LED), Psi4 (SAPT0, F-SAPT), CP2K (ALMO-EDA) and vibe-qc output maps
onto the kind, each mapping checked on a real output.
