Add explicit Molden producer dialects to `import_molden` (`dialect=`) and
`qvf convert --molden-dialect`: `standard`, `renormalize`, `orca`,
`turbomole`, `cfour`, `psi4-1.3` and `xtb`. ORCA contraction scaling and
f/g phases, Turbomole Cartesian contraction scaling, and CFOUR and old Psi4
Cartesian MO scaling map to QVF Appendix A as documented by IOData 1.0.1.
Import also accepts Molpro-style atom labels such as `O1`, integral-float
nuclear charges and `[MO]` keys in any order. Negative or excess
occupations, ECP effective charges and `[core]`/`[Pseudo]` sections are
refused with a reason. PySCF-written files and hand-built dialect variants
give PySCF's orbital values after import (WP-48, partial).
