Proposed, awaiting maintainer acceptance (GOVERNANCE.md, 2026-11-04 record):
twenty compound unit names appended to the spec § 4.3 vocabulary and
`$defs/UnitName`, from `angstrom-1` and `hartree_per_bohr` to
`e2_bohr2_per_hartree` and `states_per_ev`, each with its SI factor from the
CODATA 2018 constants the spec already uses. Fields whose fixed unit had no
name can now restate it (trajectory `velocities` and `gradient`, Fermi
`reciprocal_lattice`, the electric multipoles and polarizabilities, the
elastic compliance, the piezoelectric tensor and the electric transition
dipole), and open fields such as QTAIM `rho` and `laplacian` or a volume
grid's `data` can declare it. The 25 earlier names and every existing archive
keep their meaning. New corpus archives `units_compound_molecular.qvf` and
`units_compound_periodic.qvf`, negative archives `unit_compound_unknown.qvf`
and `unit_compound_contradicts_fixed.qvf`. No `qvf_version` or package
version change.
