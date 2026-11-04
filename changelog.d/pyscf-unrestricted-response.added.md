WP-62: `qvf.adapters.pyscf.add_td` / `write_td` now export UHF/UKS TDA,
TDHF/TDDFT and Casida response of molecular references. Contributions carry
`spin` and index the alpha or beta reference rows, normalized over both
spins. States omit `multiplicity` because unrestricted response states are
not spin eigenfunctions. NTOs use the existing unrestricted layout, one
transition per spin channel (`X_s` for TDA, `X_s + Y_s` unscaled for full
response), and the pairs reproduce PySCF's transition dipoles. Restricted-open
and spin-flip response remain refused. Uses the existing contract; no schema,
registry or validator change.
