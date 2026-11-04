Add Quantum ESPRESSO and FHI-aims converters to `qvf.convert` and
`qvf convert` (WP-48). `espresso-in` imports pw.x input structures for every
INPUT_PW `ibrav` (celldm or A, B, C) and `ibrav = 0` with `CELL_PARAMETERS`
in alat, bohr or angstrom, with `ATOMIC_POSITIONS` in alat, bohr, angstrom
or crystal, and exports `ibrav = 0` inputs. `espresso-xml` imports the final
structure, total energy, forces and stress of `data-file-schema.xml`, and
`espresso-cube` imports pp.x cube output as a periodic grid without the
endpoint plane. `aims` imports and exports FHI-aims `geometry.in` with
`lattice_vector`, `atom` and `atom_frac`. Fully fixed atoms map to
`constraints.frozen_atoms`. Everything else is reported. The tests compare
against ASE's readers and writers when ASE is installed. pp.x `filplot`
files and QE `space_group` input remain open. No format, schema or validator
change.
