WP-62: `qvf.adapters.pyscf.add_td` / `write_td` now write paired NTOs for
closed-shell triplet TDA, TDHF and TDDFT results. The `M_S = 0` spin-orbital
transitions use the existing unrestricted NTO layout (two sections per state,
α channel `T/sqrt(2)`, β channel with reversed particle sign) beside the
restricted reference; weight sums equal the singlet-path convention. Checked
against PySCF's independent UHF response amplitudes. Uses the existing NTO
contract; no schema, registry or validator change.
