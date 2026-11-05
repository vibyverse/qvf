The PySCF adapter writes k-point sampled periodic `pbc.scf.KRHF`/`KUHF` and
`pbc.dft.KRKS`/`KUKS` results (WP-62) as a periodic `structure` and a
`wavefunction.bloch` section: the k-points in fractional reciprocal-lattice
coordinates with their weights, `mesh` and `shift` for a uniform or shifted
mesh, the complex orbitals per k through the Appendix A AO map, and their
energies and occupations. PySCF's Bloch phase convention was checked to be
the spec's: S(k) built from the archive with `e^{+ik·T}` equals
`pbc_intor('int1e_ovlp', kpts=...)` to 1e-9. KROHF, KGHF, k-symmetry-adapted
objects, a non-integer spin per cell, k-points along a non-periodic axis and
k-point densities and populations are refused. `add_wavefunction_bloch` gains
optional `ecp` and `n_electrons`. The Γ-point and molecular output is
unchanged.
