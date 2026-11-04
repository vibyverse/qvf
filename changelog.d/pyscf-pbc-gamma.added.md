The PySCF adapter writes periodic Γ-point `pbc.scf.RHF`/`UHF` and
`pbc.dft.RKS`/`UKS` results of a `pbc.gto.Cell` (WP-62): a periodic
`structure` with `pbc` and `lattice_vectors` in Å, and a `wavefunction.gto`
with `k_point` `[0, 0, 0]`, the Appendix A AO map, the spec's AO `overlap`
(without images) and the densities. Populations use the lattice-summed Γ
overlap. GTH pseudopotentials are recorded by name and core count only, and
the Ewald nuclear energy is `energies.components.nuclear_repulsion`.
k-point sampling, a nonzero `kpt`, complex orbitals and periodic response
are refused. Molecular output is unchanged.
