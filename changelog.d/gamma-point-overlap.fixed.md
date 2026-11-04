Clarify spec § 5.3 for the Γ-point orbitals of a periodic structure
(proposed, pending maintainer acceptance). `overlap` stays the AO overlap
without lattice images and is checked against Appendix A.1 unchanged. The
orthonormality of the orbitals, the electron count `tr(P S)` and AO
populations refer to the lattice-summed overlap `S_Γ`, the S(k) of
`wavefunction.bloch` at k = 0, which consumers build from `basis` and the
lattice. No validator, schema or archive changes.
