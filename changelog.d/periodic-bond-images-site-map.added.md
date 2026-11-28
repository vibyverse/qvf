Periodic bond images and the symmetry site map (spec § 5.1, § 5.7). `bonds`
and `bond_orders` pairs gain the optional integer triple `image`, the lattice
image of atom `j` (absent: unspecified, `[0, 0, 0]`: the same cell), so that
bonds to neighbouring cells, several images of one pair and self-image bonds
can be written. `structure.symmetry` gains `sites`, `atom_sites` (the site,
operation and lattice translation of each atom), `site_tolerance` and
`source` (`deposited` or `detected`). The Python, C++, JavaScript and Rust
validators check both (`E_MEMBER_MALFORMED`); the corpus gains positive and
negative archives. Additive; existing archives are unchanged.
