Displacement parameters and disorder groups of crystal-structure atom sites
(spec § 5.1, *Crystallographic atom sites*). A `structure` atom gains the
optional `adp` (`$defs/AtomDisplacement`: `Uiso`, `Uani`, `Biso` or `Bani`
after CIF `_atom_site.ADP_type`; the six anisotropic components in the CIF
order, in the CIF `crystal` frame or `cartesian`; Å²) and the optional
`disorder_assembly` and `disorder_group` codes of CIF
`_atom_site.disorder_assembly` / `_atom_site.disorder_group`. The Python,
C++, JavaScript and Rust validators check that a tensor is positive definite
and that the crystal frame has a three-periodic structure
(`E_MEMBER_MALFORMED`). The Python and JavaScript writers pass the fields
through, and `qvf convert` CIF import and export carry them, each symmetry
copy with the tensor of its own orientation. Additive; existing archives are
unchanged.
