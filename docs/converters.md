# Convert common files through QVF

The `qvf.convert` Python API and `qvf convert` command started in
toolkit **0.5.0**; XSF and PDB geometry joined in **0.6.0**. The development
checkout also supports the Molden, core-CIF, mmCIF, POSCAR and CHGCAR subsets
below. They use the standard library;
NumPy is optional. Install the `validate` extra for JSON Schema checks as well
as the always-on reference checks.

Every conversion writes a new file, refuses an existing destination and
returns a JSON-compatible report. QVF imports pass the reference validator
before the completed file is published. Exports validate their QVF input and
verify member hashes. Keep the original: text exports select geometry, a
grid or orbitals, not a complete calculation. The report lists unselected
sections and explains the metadata that the text format cannot carry. This is not a
lossless archive edit; use `qvf repack` or `Editor` for that.

## XYZ and geometry-only extended XYZ

```sh
qvf convert molecule.xyz --from xyz --to qvf -o molecule.qvf
qvf convert molecule.qvf --from qvf --to extxyz -o molecule-export.xyz
qvf convert frames.xyz --from extxyz --to qvf -o trajectory.qvf
qvf convert trajectory.qvf --from qvf --to extxyz --section trajectory -o frames-export.xyz
```

Coordinates are angstrom. One frame becomes a canonical `structure`; several
frames also become a `trajectory` referring to the first structure. Atom
count, symbol order and per-axis periodicity must stay constant. Cells may
vary by frame. No molecular charge, multiplicity, energy, job or calculation
status is inferred. A geometry import alone is not a runnable queue job;
add an explicit request as described in the [family workflow](family_workflow.md).

extXYZ `Properties` is read as `name:type:count` triplets (types `S`, `R`,
`I`, `L`) and must include `species:S:1` and `pos:R:3`, in any order. Optional
`Lattice` has nine quoted, space-separated row-vector components in angstrom,
and optional `pbc` has three booleans. With a cell and no `pbc`, all axes are
periodic; with neither, all are nonperiodic. This follows the [extXYZ
conventions](https://github.com/libAtoms/extxyz#reading-aseatomsatoms-from-this-format).

The usual calculator results are kept (Unreleased, #39), converted from the
ASE units extXYZ uses to atomic units:

* the header `energy` (eV) becomes, for one frame, the root
  `energies.total` (hartree), and for several frames, the trajectory metadata
  `energies` when every frame has one;
* the `forces:R:3` column (eV/Å) becomes the energy gradient (hartree/bohr,
  negated): for one frame a `forces` section for the structure, and for
  several frames the trajectory `gradient` member when every frame has
  forces;
* the header `stress` (nine values, eV/Å³, the ASE sign convention) becomes
  that section's `stress` (hartree/bohr³) for one periodic frame.

A `Z:I:1` column must agree with the species. Other per-atom columns, such
as `tags` or `magmoms`, and other header keys, such as `config_type`, are
parsed and reported in `omitted`, not refused or stored. So are energies or
forces that only some frames carry, and the stress of a trajectory. Plain XYZ
comments are counted as omitted, not copied into archive metadata.
Unknown/ghost symbols, nonfinite values and degenerate periodic cells are
refused. Plain XYZ export refuses cells or periodicity; extXYZ retains them.
extXYZ export also writes `energy`, a `forces` column and `stress` in eV
units: for a structure, from the root energies when it is the archive's
first structure and from the one `forces` section assigned to it; for a
trajectory, from its metadata `energies` and `gradient` member. ASE reads
them back with the same values.

An archive with more than one structure/trajectory candidate requires
`--section ID`, including an imported multi-frame archive. This avoids
silently exporting its initial structure in place of its full trajectory.

## Scalar Gaussian cube

```sh
qvf convert density.cube --from cube --to qvf --kind volume.density -o density.qvf
qvf convert density.qvf --from qvf --to cube --section volume -o density-export.cube
```

The supported cube subset has a positive atom count, positive grid counts,
one scalar per voxel, bohr atom/origin/voxel coordinates and X/Y/Z loop order
(Z fastest). Fortran `D` exponents are accepted. Grid vectors can be skewed.
The implementation follows the scalar layout used by
[ASE's cube reader/writer](https://docs.ase-lib.org/_modules/ase/io/cube.html).
Atom positions convert to QVF angstrom; grid coordinates remain bohr and
scalar values are unchanged. The two comment lines and atom-charge column
are reported as omitted. A cube's atom-charge column is not assumed to be a
population analysis. Cube export writes that column as zero.

Cube does not reliably specify the field's physical meaning or periodicity.
The default is `volume.generic` and a nonperiodic structure. `--kind` is an
explicit assertion that the values already use the selected canonical
kind's units; it does not rescale the data. Multi-orbital cubes, multiple
values per voxel, negative grid counts and nonstandard loop orders are
refused rather than guessed.

For a known periodic grid, declare the axes explicitly:

```sh
qvf convert cell.cube --from cube --to qvf --periodic true true true -o cell.qvf
```

The default QVF cube endpoint convention is `includes_endpoint: true`: the
cell row is the voxel vector times `n - 1`. Add `--no-endpoint` for grids that
span the cell in `n` steps. No slab is dropped or duplicated. Export retains
voxel values and vectors verbatim and reports kind, periodic flags and the
endpoint choice, which must be supplied again on import because cube cannot
carry those QVF declarations. Never infer a crystal cell just from a cube's
bounding box.

## Static XSF geometry and scalar grids (toolkit 0.6.0)

```sh
qvf convert crystal.xsf --from xsf --to qvf -o crystal.qvf
qvf convert field.xsf --from xsf --to qvf --kind volume.density --periodic true true true -o field.qvf
qvf convert field.qvf --from qvf --to xsf --section volume -o field-export.xsf
```

The static subset follows the [XCrySDen XSF specification](https://web.mit.edu/xcrysden_v1.5.60/www/XCRYSDEN/doc/XSF.html).
It accepts one primitive geometry (`PRIMVEC`/`PRIMCOORD`) with `CRYSTAL`,
`SLAB` or `POLYMER`, or a molecular `ATOMS` geometry, optionally preceded by
`MOLECULE`. Periodic structure axes must be the leading cell axes (three,
two or one respectively); the converter never permutes axes. Atom identities
may be symbols or atomic numbers. Positions and cell rows remain angstrom.
Atom rows may carry forces (`Fx Fy Fz` in hartree/Å, the XSF convention) on
every row or none; they become a `forces` section whose gradient is in
hartree/bohr (negated). Structure export writes the forces of the one
`forces` section assigned to the structure back, and ASE reads them with the
same values. `CONVVEC` and `CONVCOORD` (the conventional cell) are accepted,
checked and reported in `conventional_cell_omitted`; the primitive geometry
is the one imported.

Every named scalar `DATAGRID_3D`, in one or several `BEGIN_BLOCK_DATAGRID_3D`
blocks, becomes its own volume section (`volume`, `volume_2`, …) labelled
`<block>/<grid>`; `--kind` and `--periodic` apply to all of them, and the
report's `grids` lists them. Grid origin and
spanning vectors convert from angstrom to bohr. Each voxel row is its span
divided by `n - 1`, so both boundary planes remain present. Values move from
XSF's X-fastest order to QVF's Z-fastest order without changing samples.
No plane is added, removed or resampled. Grid shape needs at least two
points per axis. A periodic QVF grid without endpoint planes is refused on
export; it must first be explicitly transformed by the caller.

A periodic structure does not prove that a particular scalar field repeats.
Import defaults to a finite grid, even with `CRYSTAL`; `--periodic` explicitly
asserts the **field's** repeating axes, which must match the structure axes
and cell spans. The default kind is `volume.generic`. `--kind` asserts the
values already use that kind's canonical units, with no scalar rescaling.
Export reports kind and field periodicity because XSF cannot encode those
QVF declarations; supply them again when reimporting. On nonperiodic axes,
QVF's endpoint flag has no effect; XSF always retains the sampled extent.

Animation (`ANIMSTEPS`), 2D and band grids, empty grid blocks, degenerate
boxes and unknown records are refused. Comments between sections are
reported as omitted. A nonperiodic
structure with a cell cannot be exported by this molecular subset. As with
other converters, a file containing both a structure and a grid requires
`--section` on export. No charge, multiplicity or job is inferred.

## PDB structures (toolkit 0.6.0; atom-site identity Unreleased)

```sh
qvf convert molecule.pdb --from pdb --to qvf -o structure.qvf
qvf convert crystal.pdb --from pdb --to qvf --expand-symmetry -o crystal.qvf
qvf convert virus.pdb --from pdb --to qvf --expand-ncs -o virus.qvf
qvf convert structure.qvf --from qvf --to pdb -o structure-export.pdb
```

This bounded subset follows the [wwPDB coordinate records](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
and [fixed-column conventions](https://www.wwpdb.org/documentation/file-format-content/format33/sect1.html).
Import reads `ATOM` and `HETATM` records in source order: explicit elements
(columns 77–78) and Cartesian angstrom coordinates (columns 31–54,
`Real(8.3)`). Names such as `CA` never determine the element. ASCII records,
unique positive decimal serials within a model and recognized element
symbols are required; isotope labels such as `D`, hybrid-36 serials or
residue numbers, and missing elements are refused.

Toolkit 0.6.0 kept only elements and coordinates and refused alternate
locations and ensembles. The development checkout also keeps each atom's
atom-site identity in the typed fields of spec § 5.1 (*Biomolecular atom
identity*):

| PDB columns | Per-atom field |
|---|---|
| 13–16 | `atom_name` |
| 17 | `alt_loc` |
| 18–20 | `residue_name` |
| 22 | `chain_id` |
| 23–26 | `residue_seq` (signed decimal) |
| 27 | `insertion_code` |
| 55–60 | `occupancy` (0 to 1; outside is refused) |
| 61–66 | `b_factor` |
| `ATOM` / `HETATM` | `hetero` (`false` / `true`) |

`HELIX` and `SHEET` records become the section's `secondary_structure` ranges
(spec § 5.1, *Secondary-structure ranges*), in file order: `helix` ranges with
their `helix_class` (columns 39–40, 1–10) and `sheet` ranges with `sheet_id`
(columns 12–14) and `strand` (columns 8–10). Each end keeps its insertion code
as `start_insertion_code` / `end_insertion_code`, so a helix from 10A starts at
10A. A range must start and end on one chain. A range whose end is not among
the imported residues, such as a disordered loop, is kept as written; the
report counts these in `secondary_structure_unresolved`. Helix identifiers,
comments and lengths, and strand sense and registration, are not kept.

`CONECT` records become the structure section's own `bonds` member
(`$defs/BondsPayload`, spec § 5.1): each distinct pair of bonded atoms is one
bond, as 0-based atom positions, sorted. PDB defines no bond order, so
`order` stays absent (meaning 1). A pair listed again from the same atom, a
convention some tools use for multiple bonds, is counted in
`conect_repeated_listings` and not read as an order. A `CONECT` that names a
serial with no atom, or bonds an atom to itself, refuses the import. The
obsolete hydrogen-bond and salt-bridge columns 32–61 are counted in
`conect_hydrogen_bond_fields_omitted`. In an ensemble written as one structure
per model, serials resolve within each model, and a serial missing from any
model refuses the import. In the structure-plus-trajectory case the bonds sit
on the structure, which names every frame's atoms.

`SSBOND` and `LINK` records add their bonds to the same member, merged with
the `CONECT` pairs, again without an order. Which record listed a bond is not
kept; the report's `ssbond` and `link` give each record count, the distinct
bonds it resolved to, and how many of those `CONECT` also lists. Atoms
resolve by atom-site identity: chain, residue number, insertion code, residue
name and atom name, with `SSBOND` naming the `SG` atom of each cysteine. A
blank alternate location names every alternate of the atom, and a given one
names that alternate, or the atom if it has no alternates; alternates bond
only to the same alternate or to an atom without alternates, so a disulfide
between two two-conformer cysteines is two bonds, A–A and B–B. An atom that
is absent or ambiguous, a self-bond, or two ends with no common alternate
refuses the import. A bond pair cannot point to a periodic image, so a record
whose symmetry operator is not `1555` (blank means `1555`) is skipped and
counted, with its operators, in `symmetry_mate_bonds_skipped`. Serial
numbers and the distance column are not kept. Ensembles resolve the records
within each model, as for `CONECT`.

A `CRYST1` cell makes the structure periodic: `lattice_vectors` in the PDB
standard orientation (a along x, b in the xy-plane, c completing a
right-handed cell) and `pbc` `[true, true, true]`. A `structure.symmetry`
section (`symmetry`, or `symmetry_<section>` per model) names the structure
and records the CRYST1 space group as `space_group_symbol`. PDB coordinates
are usually the asymmetric unit: symmetry mates are **not** generated, and the
report's `cell` says so and gives Z. The wwPDB placeholder
`1.000 1.000 1.000 90.00 90.00 90.00 P 1` (Z 1 or blank), written by NMR,
cryo-EM and model files, means no cell: the structure stays nonperiodic and
`placeholder_cell` is true. `SCALE1`–`SCALE3`, when present, must equal the
standard orthogonalization of the cell (relative tolerance 1e-4) with zero
translation. A different matrix puts the coordinates in another frame, so it
refuses the import, as do SCALE records without CRYST1, a second CRYST1, or
angles that describe no cell. In an ensemble every model gets the cell; the
structure-plus-trajectory case also writes a per-frame `lattice` member.

Each `REMARK 350` biomolecule becomes its own nonperiodic `structure` section,
`assembly_<n>`, labelled `Biological assembly <n>`; the deposited structure is
left unchanged. QVF has no field for assembly operators, so the assembly is
written out: every `BIOMT` operator of an `APPLY THE FOLLOWING TO CHAINS`
block (with its `AND CHAINS` continuations) is applied to the atoms of those
chains, block by block and operator by operator, with coordinates rounded to
1e-6 Å. A copy under an identity operator keeps its chain id; any other copy
is named `<chain>-<BIOMT serial>`, such as `A-2`, so each atom's origin stays
readable and residues stay distinct. Atom-site identity is copied, bonds are
copied when both atoms are in the same copy, and secondary-structure ranges
are copied to the renamed chains. Assemblies are built from the first model
of an ensemble. The report's `assemblies` gives each section, its atom count
and its copies with chains and the 3×4 operator `[R | t]`; the operators are
not stored in the archive. An operator that is not a proper rotation, a
missing `BIOMT` row, a chain with no atoms, or a chain copied twice under one
name refuses the import. Other `REMARK 350` lines, such as how the unit was
determined, are not kept.

`MTRIX` records (non-crystallographic symmetry) are counted in the report's
`mtrix`. With `iGiven` 1 the copy is already in the file. When a copy is not
given, the file holds only part of the asymmetric unit: the structure stays
as deposited, NCS copies are not generated, no assembly is built, and
`assemblies_not_built` says why. An `MTRIX` operator that is not a proper
rotation, or whose rows disagree on `iGiven`, refuses the import.

Generating the missing copies is opt-in: `import_pdb(..., expand_ncs=True)`
or `--expand-ncs`. Each `MTRIX` operator whose copy is not given (skipping an
identity) is applied to the deposited atoms of the first model,
x′ = R x + t. The deposited atoms followed by the copies form a new
`structure` section, `asymmetric_unit`, labelled
`Asymmetric unit (MTRIX copies generated)`, with the deposited cell and its
own `structure.symmetry` section; the deposited structure is unchanged.
Copied chains are named `<chain>-<MTRIX serial>`, and bonds and
secondary-structure ranges are copied with them. `REMARK 350` assemblies
and the opt-in unit-cell expansion are then built from this completed
asymmetric unit instead of being blocked. Assemblies still apply their
`BIOMT` operators only to the chains `REMARK 350` lists. The report's
`ncs_expansion` gives the section, its atom count and each operator. The
copies agree with gemmi's `expand_ncs`.

Crystal symmetry expansion is opt-in: `import_pdb(..., expand_symmetry=True)`
or `qvf convert --from pdb --to qvf --expand-symmetry`. It adds a periodic
`structure` section, `unit_cell`, labelled `Unit cell (<group> expanded)`,
with the deposited lattice and every copy of the first model under the
space group's operators, plus a `structure.symmetry` section
`symmetry_unit_cell` recording `P 1`; the deposited structure is unchanged.
Each copy is shifted by whole lattice vectors so its fractional centre lies
in [0, 1), so molecules stay intact. The operators come from a built-in Hall
symbol parser (Hall, 1981) and cover the 65 chiral space groups under their
wwPDB names, with the settings `I 1 2 1`, `A 1 2 1`, the other axis orders of
`P 2 2 21` and `P 21 21 2`, `H 3`/`H 3 2` (hexagonal axes) and `R 3`/`R 3 2`
(rhombohedral axes), and the short `P 2`, `P 21`, `C 2` and `I 2`. Every
table entry is checked against gemmi. Copy 1 is the identity and keeps its
chains; copy *n* is named `<chain>-<n>` (`-<n>` for a blank chain). Identity,
bonds within each copy and secondary-structure ranges are copied. A copy of
an atom that lands within 0.1 Å of an earlier copy of the same atom, modulo
the lattice (an atom on a symmetry element), is merged and counted in
`special_position_atoms_merged`; a bond to it is kept when the two copies
coincide in place and otherwise counted in `bonds_dropped`. The report's
`symmetry_expansion` gives the group, Hall symbol, operator and atom counts.
Expansion is refused without a real `CRYST1` cell, for a space group outside
the table, and when an `MTRIX` copy is not given unless `expand_ncs`
generates it first. A large cell can hold
hundreds of thousands of atoms: 5,000 atoms in `F 41 3 2` (96 operators)
take about 30 s and 21 MB.

A blank column is left absent, never written as a placeholder. When every
atom has a residue name and number, the section also lists `chains` (first
appearance order, blank chain ids excluded) and `residues` in file order,
keyed by (chain, residue number, insertion code). So residues 10 and 10A are
two residues, and a site whose alternate locations carry different residue
names is one. Every alternate location is its own atom with its own
occupancy. Nothing is selected or dropped by `alt_loc` or occupancy; select a
conformer in the consumer, as a view.

One `MODEL`/`ENDMDL` block, or none, gives one `structure`. A multi-model
ensemble is imported whole. When every model has the same ordered atoms with
the same per-atom metadata, model 1 is the `structure` and a `trajectory`
(`structure_ref: "structure"`) holds every model's coordinates, with the model
serials in its metadata `model_serials`. Otherwise each model is its own
`structure` section, `model_<serial>`, labelled `MODEL <serial>`. A `NUMMDL`
record must agree with the model count, and atoms outside `MODEL`/`ENDMDL`
are refused once models are used.

The report gives the sections written, the model serials, the number of
alternate-location atoms, whether residues were listed, counts of omitted
metadata records, atom-charge fields, and the secondary-structure ranges.
Every model's structure section carries the same ranges. Standard metadata
such as `ORIGX`, `SEQRES`, `CISPEP` and other `REMARK` records is counted but
not interpreted, and unknown records fail. No crystal symmetry expansion
(unless requested), NCS copy, bond order, protonation, charge, multiplicity or job is inferred. Atom serials, segment identifiers and atom charges are not kept.
Keep the original PDB for crystallographic and deposition information.

Export selects one structure, with `--section` required if a trajectory or
another structure is also present; for an imported ensemble, export one
model at a time. Trajectories are refused. A structure periodic on all three
axes is exported with `CRYST1` and `SCALE1`–`SCALE3` when its lattice is in
the PDB standard orientation (off-diagonal upper entries within 1e-4 Å, a
positive diagonal). The space group comes from a `structure.symmetry` section
assigned to it, else `P 1`. The report gives the space group and the largest
cell rounding error. A rotated lattice, partial periodicity, or a cell
without periodicity is refused; coordinates are never rotated, and extXYZ
remains the export for those. When every
atom has `atom_name`, `residue_name` and `residue_seq`, export writes the
identity back: `HETATM` for atoms marked `hetero` (otherwise `ATOM`), the
standard atom-name alignment, chain, insertion code, alternate location, and
occupancy and temperature factor, left blank where the archive records none.
The section's `helix` and `sheet` ranges become `HELIX` and `SHEET` records
before the atoms, with residue names taken from the atoms (blank for an end
the structure does not contain). Strand 1 gets sense 0 and other strands a
blank sense; a sheet range without `sheet_id` becomes its own one-strand
sheet. `coil` ranges have no PDB record and are counted in
`coil_ranges_omitted`.

Both forms of export write the structure's bonds as `CONECT` records after
the atoms. The bonds are the structure's own `bonds` member plus every
`bonds` section that § 5.1 assigns to it. Each bond is listed from both
atoms, four partners per record, by sequential serial. Bond orders are not
written: the report counts orders other than 1 in `bond_orders_not_written`,
and the remark says `CONECT, NO BOND ORDERS`. A bond whose atom index is
outside the structure is refused. With atom-site identity, export also writes
a bond between different residues as `SSBOND` when it joins the `SG` atoms
of two cysteines (one record per residue pair, so the two conformer bonds
of a disulfide become one record) and as `LINK`, with alternate locations,
when another such bond involves a hetero atom. Both come before the atoms,
with operators `1555` and the bond distance as the length; the report counts
them in `ssbond_records` and `link_records`. Peptide and other bonds between
standard residues stay `CONECT` only. gemmi reads both back as connections.
It writes up to 99,999 atoms with sequential serials. Values that do not fit
their columns are refused rather than truncated: names over four characters,
residue names over three, multi-character chain ids (so an assembly with
copied chains such as `A-2` is refused; export it to extXYZ instead), or
residue numbers outside −999..9999. `REMARK 350` is not written. Importing that file again reproduces the structure
payload and its `chains`/`residues`. A structure without that identity is
written as before: up to 9,999 synthetic `HETATM` records, one `UNL` residue
per atom, element-based names, occupancy 1.00 and temperature factor 0.00.
These are placeholders, not biological or experimental metadata. Coordinates
round to 0.001 angstrom, and the report gives the maximum rounding error and
whether atom-site identity was written. Neither output is a deposition-ready
PDB.

## Crystal CIF structures (Unreleased)

```bash
qvf convert crystal.cif --from cif --to qvf -o crystal.qvf
qvf convert several.cif --from cif --to qvf --block rutile -o rutile.qvf
qvf convert crystal.qvf --from qvf --to cif -o crystal-p1.cif
```

Import reads core crystallographic CIF 1.1, the small-molecule and inorganic
dialect: the syntax with comments, quoted values and semicolon text fields,
data blocks and loops. `save_` frames and `global_` blocks are refused, and
mmCIF (PDBx) categories are not interpreted. One data block is imported: the
only one with an `_atom_site_label` loop, or the one named by `block=` /
`--block` (case-insensitive, without `data_`). Several candidates without a
choice refuse the import.

The block needs the six `_cell_length_*` / `_cell_angle_*` values and an
atom-site loop with `_atom_site_label`, `_atom_site_type_symbol` and
`_atom_site_fract_x`/`y`/`z`; sites given only in Cartesian coordinates are
refused. A number's standard uncertainty, as in `4.5937(3)`, is dropped. The
element comes from `_atom_site_type_symbol` with any charge stripped (`Fe3+`
is Fe, counted in `charges_omitted`); a missing or non-element type symbol
refuses the import, because no element is guessed from labels. The label
becomes the atom's `atom_name`, `_atom_site_occupancy` its `occupancy` (in
(0, 1]; partial sites are counted in `partial_occupancy_sites`), and
`_atom_site_B_iso_or_equiv`, or 8π² × `_atom_site_U_iso_or_equiv`, its
`b_factor`. Every site is imported; disordered alternatives are not resolved.

The structure is the whole unit cell. The symmetry operations come from
`_space_group_symop_operation_xyz` (or the older `_symmetry_equiv_pos_as_xyz`).
They must be distinct, include `x,y,z` and form a group. Otherwise they come
from `_space_group_name_Hall` (or `_symmetry_space_group_name_Hall`) through
the built-in Hall-symbol parser, which reproduces gemmi's operators for all
564 settings in its table. Without either, only a group named `P 1` (or
number 1) is accepted. Each operation is applied to each fractional site and
the result is wrapped into [0, 1). An image within 0.1 Å (minimum image) of
an earlier image of the same site, a site on a special position, is merged and
counted in `special_position_atoms_merged`. Different sites at one place are
kept. The lattice is in the standard orientation (a along x, b in the
xy-plane) with `pbc` `[true, true, true]`, and Cartesian positions are rounded
to 1e-6 Å. A `structure.symmetry` section `symmetry` keeps the space-group
symbol, number and Hall symbol the block gives. The report gives the block,
the operator count and source (`listed`, `hall` or `p1`), the atom count and
the number of other tags not read (`unused_tags`). Bonds, anisotropic
displacement, geometry tables and experimental metadata are not kept.

Export writes one structure periodic on all three axes as a P 1 CIF: the
cell lengths and angles (6 decimals), the identity operation and fractional
coordinates (8 decimals) for every atom, with `_atom_site_occupancy` and
`_atom_site_B_iso_or_equiv` columns when any atom has them (`?` where one
does not). Labels are the atoms' `atom_name` values when they are unique and
CIF-safe, else element plus a running number (`Ti1`, `O1`, …); the report
says which. Any lattice orientation is accepted, but a reader places the cell
in its own standard orientation, so Cartesian coordinates may come back
rotated; the report gives the largest coordinate rounding error. A
structure without full periodicity is refused. gemmi reads the export back,
and importing it again reproduces the positions.

## mmCIF (PDBx) structures (Unreleased)

```bash
qvf convert 1abc.cif --from mmcif --to qvf -o 1abc.qvf
qvf convert 1abc.qvf --from qvf --to mmcif -o 1abc-export.cif
```

Import reads one mmCIF data block with an `_atom_site` category: the only
one, or the one named by `block=` / `--block`. It uses the CIF syntax rules
of the core-CIF import. Each `_atom_site` row becomes an atom with the PDB
atom-site identity of spec § 5.1. The author columns are preferred, as in
PDB files: `auth_atom_id`, `auth_comp_id`, `auth_seq_id` and `auth_asym_id`,
falling back to the `label_*` columns. With them come `pdbx_PDB_ins_code`,
`label_alt_id`, `occupancy`, `B_iso_or_equiv` and `hetero` from `group_PDB`.
The element comes only from `type_symbol`, and positions from
`Cartn_x`/`y`/`z`. Rows are grouped by `pdbx_PDB_model_num` and follow the
PDB ensemble rules: one structure, a structure plus a `trajectory` when every
model has the same atoms and metadata, or one `model_<n>` section per model.
`_cell` with `_symmetry.space_group_name_H-M` (or
`_space_group.name_H-M_alt`) becomes the periodic cell and
`structure.symmetry` section exactly as `CRYST1` does, including the
placeholder rule. `_struct_conf` helices (`HELX_*`) and `_struct_sheet_range`
strands become `secondary_structure` ranges with insertion codes, helix
class, sheet id and strand number; other `_struct_conf` types, such as turns,
are counted. `_struct_conn` `disulf`, `covale` and `metalc` connections become
bonds resolved like `SSBOND` and `LINK`. `hydrog` rows are counted in
`hydrogen_bonds_omitted`, and connections to a symmetry mate (an operator
other than `1_555`) in `symmetry_mate_bonds_skipped`. Assemblies
(`pdbx_struct_assembly_gen`), entities, label numbering and the rest of the
dictionary are not interpreted. A gemmi-written mmCIF of a PDB file imports
to the same structure, residues and symmetry as the PDB import of that file.

Export writes one structure as an `_atom_site` loop. The label columns repeat
the author identity, ids are sequential, the model number is 1, and Cartesian
coordinates have 4 decimals. A structure periodic on all three axes, with its
lattice in PDB orientation, also gets `_cell` and
`_symmetry.space_group_name_H-M` (from its `structure.symmetry`, else
`P 1`). Values are quoted as CIF requires. Unlike PDB, there is no limit on
chain-id length or atom count, so an assembly with copied chains such as
`A-2` can be exported. Bonds and secondary structure are not written. gemmi
reads the export back, and importing it again reproduces the structure.

## VASP POSCAR/CONTCAR structures (Unreleased)

```bash
qvf convert POSCAR --from poscar --to qvf -o cell.qvf
qvf convert POSCAR-vasp4 --from poscar --to qvf --species Ga N -o cell.qvf
qvf convert cell.qvf --from qvf --to poscar -o POSCAR.out
```

Import reads one VASP POSCAR or CONTCAR as a periodic `structure`
(`pbc` `[true, true, true]`), labelled with the comment line. The scale line
is one factor, a negative value giving the cell volume in Å³, or three
per-axis factors (VASP 6); the lattice rows, and Cartesian positions, are
multiplied by it. Elements come from the VASP 5 species line, where POTCAR
labels such as `Ga_d` or `Si/1a2b` reduce to their element. A VASP 4 file,
with counts only, needs `species=` / `--species`, one element per count; it is
never guessed from the comment. A species line that disagrees with
`species=` refuses the import. Positions are Direct (fractional) unless the
mode line starts with C or K (Cartesian), as in VASP. Trailing labels on a
position line are ignored.

`Selective dynamics` flags have no per-axis QVF field. An atom fixed on all
three axes (`F F F`) goes into the root `constraints.frozen_atoms` (spec
§ 3.3), which optimizers read; an atom fixed on only some axes is counted in
`partially_fixed_atoms_omitted` and not stored. A CONTCAR's velocity block,
with or without a coordinate-mode line, must have a row of three numbers per
atom. It is reported in `velocities_omitted` and not stored. A
`Lattice velocities and vectors` block is reported too, and the
predictor-corrector lines after the velocities are counted. Any other
content after the positions refuses the import, as do zero or missing
counts, a zero scale, linearly dependent lattice vectors, malformed flags and
too few position lines. Magnetic moments, charge and the job are not
inferred.

Export writes one structure periodic on all three axes in VASP 5 layout:
the section label (or id) as comment, scale `1.0`, the archive's lattice
rows unrotated, and Direct coordinates with 10 decimals in archive order.
The species and count lines have one entry per consecutive run of an
element, so `Si O Si` is written as three runs rather than reordered. When
the archive's root `constraints.frozen_atoms` applies, because the exported
structure is the archive's first, the file gets `Selective dynamics` with
`F F F` for those atoms and `T T T` for the rest. Constraints that cannot be
written are listed in `constraints_not_written`. The report gives the
species runs and the largest coordinate rounding error. ASE reads the export
back with the same cell, positions and fixed atoms, and ASE-written POSCAR
files import to the same structure.

## VASP CHGCAR charge densities (Unreleased)

```bash
qvf convert CHGCAR --from chgcar --to qvf -o density.qvf
qvf convert density.qvf --from qvf --to chgcar -o CHGCAR.out
```

CHGCAR, CHG, AECCAR0/1/2 and PARCHG share one layout, and all import the
same way. The POSCAR head becomes the periodic `structure` under the POSCAR
rules above, including `species=` / `--species` for a VASP 4 head and
`frozen_atoms` from Selective dynamics. Each grid follows a line with its
three point counts; its values are ρ·V_cell (V in Å³), x fastest. The first
grid becomes the `volume.density` section `density` in e/bohr³ (divided by
the cell volume in bohr³). Its grid is periodic on all three axes, starts at
the origin and has `includes_endpoint` false, so `voxel_vectors[a]` is
lattice vector a over `shape[a]` (spec § 5.2). A spin-polarized file has a
second grid, the magnetization density ρ↑ − ρ↓, which becomes the
`volume.spin` section `magnetization` in the same unit. A noncollinear file
has three more grids (mx, my, mz), which become one `volume.vector` section
`magnetization` of shape `[nx, ny, nz, 3]`. Any other number of grids, grids
of different shapes, or too few values refuse the import. PAW augmentation
occupancies and per-atom moment lines between grids are skipped and counted
in `lines_skipped`. The report gives the spin case, the cell volume, the
integrated electron count (`electrons`, the mean of ρ·V over the grid) and
the integrated magnetization. Grids are limited to five million points.

Export selects one `volume.density` whose grid spans its periodic structure
CHGCAR-style: periodic on all three axes, origin at 0, `includes_endpoint`
false and `voxel_vectors[a] × shape[a]` equal to lattice vector a. Other
grids are refused. The structure is written as in the POSCAR export, then
the density as ρ·V_cell, five values per line. A single `volume.spin` or
`volume.vector` section on the identical grid is appended as the
magnetization; with several candidates none is written, and the report
lists them. No augmentation occupancies are written, so the file is for
analysis and visualization tools, not a VASP restart. ASE's
`VaspChargeDensity` reads the export back with the same densities, and
ASE-written CHGCAR files import to the same grid.

## Molecular Molden orbitals (Unreleased)

```sh
qvf convert orbitals.molden --from molden --to qvf --spin restricted -o orbitals.qvf
qvf convert orbitals.qvf --from qvf --to molden --section wavefunction -o orbitals-export.molden
```

This subset follows the [Molden format description](https://www.theochem.ru.nl/molden/molden_format.html)
for molecular geometry, atom-centered Gaussian contractions and real orbitals.
It supports `s`, `p`, combined `sp`, `d`, `f` and `g` shells, each spherical
or Cartesian. The default is Cartesian; `[5D]` and `[5D7F]` make d and f
spherical, `[5D10F]` makes d spherical and keeps f Cartesian, `[7F]` makes f
spherical, and `[9G]` / `[15G]` make g spherical or Cartesian (Unreleased,
#41 for f and g). Conflicting flags are refused. Spherical components in
Molden order (m = 0, +1, −1, …, +l, −l) map to QVF's m = −l … +l without sign
changes. Cartesian components in Molden's `10F`/`15G` orders map to QVF's
lexicographic order, rescaled from Molden's unit-normalized components to
QVF's axial normalization (Appendix A.1).
Atoms need sequential indices and matching element/nuclear-charge values.
Coordinates require `AU` or `Angs`; basis exponents and orbital energies use
atomic units. Every orbital needs `Sym`, `Ene`, `Spin`, `Occup` and one real
coefficient for every AO. Sparse or repeated AO indices are refused.

Import requires `--spin restricted`, `restricted_open` or `unrestricted`.
An Alpha-only file cannot distinguish the first two; the caller must choose.
Unrestricted import requires both Alpha and Beta sets and occupations in
0–1; a single set uses 0–2. The resulting `structure` and `wavefunction.gto`
use existing QVF v1 members. Orbital kind is `custom`: Molden does not prove
that orbitals are canonical, localized or natural. No method, charge,
multiplicity, electron count or job is inferred.

The normalized-AO convention is checked independently with
[PySCF's Molden reader/writer](https://pyscf.org/_modules/pyscf/tools/molden.html).
Input contractions must have squared norm within `1e-5` of one; small text
rounding drift is normalized and the original norms are reported. Ordering
and Cartesian component scaling convert to QVF Appendix A. Export normalizes
contractions and compensates MO coefficients so orbital values are preserved;
it does not orthogonalize orbitals or claim a wavefunction is physically valid.
Export writes coordinates in bohr to avoid a second angstrom conversion.

ECP/core sections, ghosts, Slater functions, h and higher shells, complex or
spinor coefficients, periodic geometry, transition-weight occupations and
unknown sections are refused. Non-unit shell scale factors are refused, and
unnormalized input contractions are refused unless an explicit producer
dialect (below) accounts for them.
Export requires a basis for every atom, contiguous shells per atom and one
convention each for d, f and g, and writes the matching flags. PySCF-written
f and g files give the same orbital values at off-axis points after import
and after export. More than one wavefunction requires `--section`.

### Molden dialects (Unreleased)

Programs that write Molden files deviate from the format description in
documented ways. Import never guesses the producer from orbital norms or the
title: name it with `dialect=` (Python) or `--molden-dialect` (CLI). The
default is `standard`. Each dialect maps the written values onto standard
Molden before the Appendix A conversion above, and the report records it.
The ORCA, Turbomole, CFOUR and Psi4 rules are the documented corrections in
[IOData](https://github.com/theochem/iodata) 1.0.1 (`iodata/formats/molden.py`).
In this table N_i is the Appendix A.1 primitive norm and Π is
(2l_x−1)!!·(2l_y−1)!!·(2l_z−1)!! for one Cartesian component.

| Dialect | Producers | Written values, mapped back on import |
|---|---|---|
| `standard` | Molden format description, PySCF `tools.molden`, Psi4 1.4 and later, Molpro, Multiwfn | Normalized contractions, unit-normalized AOs, Molden phases. Squared contraction norms must be within 1e-5 of one. |
| `renormalize` | Any producer with unnormalized contractions and standard AOs | Every contraction is rescaled to unit norm, as the Molden program does when it reads a file. Norms are reported. |
| `orca` | ORCA `orca_2mkl -molden` | Each contraction coefficient is multiplied by the norm of the Cartesian primitive (0,0,0), (1,0,0), (1,1,0), (1,1,1) or (2,1,1) for s, p, d, f, g. The spherical f and g components m = ±3 and g components m = ±4 have the opposite sign. Cartesian d, f and g shells are refused. |
| `turbomole` | Turbomole `tm2molden` | Cartesian d, f and g contraction coefficients are divided by √3, √15 and √105. |
| `cfour` | CFOUR up to 2.1 | Cartesian d, f and g MO coefficients are multiplied by 1/√Π. |
| `psi4-1.3` | Psi4 up to 1.3.2 | Cartesian d, f and g MO coefficients are multiplied by √((2l−1)!!/Π), so they refer to axially normalized AOs. |
| `xtb` | xTB `printmold` (used by `qvf export-xtb`) | Contraction coefficients include N_i, and odd-l MO coefficients have the opposite sign (xTB evaluates R_atom − r). |

After the mapping, the squared contraction norm check applies to every
dialect except `renormalize`, so a wrong dialect choice is refused when the
norms disagree. A refusal names the contraction and its norm, and suggests
`orca` when the title names ORCA. Psi4 1.0 and earlier, which multiplied
primitive norms into the contraction coefficients, has no dialect.

These text forms are accepted in every dialect:

| Item | Accepted forms |
|---|---|
| Header | `[Molden Format]` in any letter case, followed by optional `[Title]` text |
| Units | `[Atoms] AU`, `Angs`, `(AU)` or `(Angs)` in any letter case |
| Atom labels | An element symbol in any case, optionally followed by a tag such as `O1` or `HE12` (Molpro style). The atomic number must match the element |
| Angular flags | `[5D]`, `[5D7F]`, `[5D10F]`, `[7F]`, `[9G]`, and `[6D]`, `[10F]`, `[15G]`, alone or combined, before or after `[GTO]`, one per line, in any letter case. Absent flags mean Cartesian |
| `[GTO]` | Atom lines `n 0` or `n`, with or without blank lines between atoms and at the end, shell scale `1.00`, `1.0` or absent |
| Numbers | Fortran exponents with `D` or `d`, as in `0.2000D+01`, everywhere |
| `[MO]` | `Sym=`, `Ene=`, `Spin=` and `Occup=` in any order and letter case, each once per orbital. `Sym` labels are kept as written (`1a1`, `1.1`, `a'`). Occupations may be fractional |

These inputs are refused with a message saying what QVF cannot represent:

- negative occupations (for example natural-orbital occupations from a
  correlated calculation) and occupations above 2 (restricted) or 1
  (unrestricted), with no clipping,
- nuclear charges that differ from the element (ECP effective cores, ghosts,
  point charges),
- `[core]` and `[Pseudo]` sections, whose ECP parameters are not in the file,
- missing `Sym`, `Ene`, `Spin` or `Occup` keys and sparse coefficient lists,
- h and higher shells, Slater functions and unknown sections.

Import reports omitted comments/title text. Export reports its selected spin,
contraction norms, omitted sections/members and generated `A` symmetry labels
when labels were absent. It retains geometry, basis, real coefficients,
energies, occupations and labels; QVF identities, bonds, method metadata,
jobs and overlap/density matrices are not retained. Supply the reported spin
on reimport and keep the original archive.

## Gaussian formatted checkpoint files (Unreleased)

```sh
qvf convert water.fchk --from fchk --to qvf -o water.qvf
qvf convert water.qvf --from qvf --to fchk --section wf -o water-export.fchk
```

```python
from qvf.convert import import_fchk, export_fchk

report = import_fchk("water.fchk", "water.qvf")
report = export_fchk("water.qvf", "water-export.fchk")
```

The file layout follows Gaussian's
[formchk documentation](https://gaussian.com/formchk/): a title line, a
line with the job type, method and basis name, then named fields. A field
header gives the name, the type `I` (integer) or `R` (real) and either one
value or `N=` and a count, followed by the values. Text (`C`, `H`) and
logical (`L`) fields are skipped and listed in the report. Reals may use `D`
exponents and the Fortran form without `E` for a three digit exponent.

| fchk field | QVF |
|---|---|
| `Atomic numbers`, `Current cartesian coordinates` (bohr) | `structure` atoms, positions in Å |
| `Charge`, `Multiplicity`, `Number of electrons` | root `provenance` `charge`, `multiplicity`, `n_electrons`, and `basis.n_electrons` |
| `Shell types`, `Number of primitives per shell`, `Shell to atom map`, `Primitive exponents`, `Contraction coefficients`, `P(S=P) Contraction coefficients` | `wavefunction.gto` `basis` |
| `Alpha MO coefficients`, `Alpha Orbital Energies` and the Beta pair | `mo_coefficients`, or `mo_coefficients_alpha` and `mo_coefficients_beta`, and the `mo_metadata` energies |
| `Total SCF Density`, `Spin SCF Density` | `density_matrix`, or `density_matrix_alpha` and `density_matrix_beta` |
| `Total Energy` | root `energies.total`, with the method of the second line in lowercase |
| `Mulliken Charges` | an `atom_properties` section with `mulliken_charge` |
| basis name of the second line | `provenance.basis` |

**Shells.** Shell type 0 is s, 1 is p, `-1` is SP, a positive type l ≥ 2 is
a Cartesian shell and a negative type `-l` a pure (spherical) one, up to
l = 7. An SP shell becomes an s and a p shell on the same exponents. The
`Pure/Cartesian d shells` and `Pure/Cartesian f shells` flags (0 pure,
1 Cartesian, the f flag for f and higher) must agree with the signs. The AO
order within a shell maps onto Appendix A.2:

| Shell | Gaussian order | QVF order |
|---|---|---|
| p | x, y, z | y, z, x when the file has pure shells (s and p then carry `pure: true`), else x, y, z |
| pure, l ≥ 2 | m = 0, +1, −1, +2, −2, …, +l, −l | m = −l, …, +l |
| Cartesian d | xx, yy, zz, xy, xz, yz | xx, xy, xz, yy, yz, zz |
| Cartesian f | xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz | xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz |
| Cartesian g and higher | the reverse of the QVF order: zzzz, yzzz, yyzz, …, xxxx | xxxx, xxxy, …, zzzz |

**Normalization.** Every Gaussian basis function is normalized.
Contraction coefficients refer to normalized primitives, and a squared
contraction norm must be within 1e-5 of one: an unnormalized contraction is
ambiguous and refused. Text rounding drift is normalized away and the norms
are reported. Each Cartesian component is normalized on its own, so the
Gaussian d_xy is √3 times the QVF d_xy of Appendix A.1, f_xyy is √5 and
f_xyz √15 times the QVF function. MO coefficients and density matrices are
divided by these factors. Pure functions are the real solid harmonics of
Schlegel and Frisch, *Int. J. Quantum Chem.* **54**, 83 (1995), the
convention of Appendix A.2, so no signs change. These conventions are the
ones [IOData](https://github.com/theochem/iodata) 1.0.1 documents for fchk
(`iodata/formats/fchk.py`, `CONVENTIONS`, with L2 normalized primitives).
No IOData, cclib or Gaussian code is used.

**Spin and occupations.** fchk stores no occupations. A file with
`Beta MO coefficients` is `unrestricted`. Otherwise equal alpha and beta
electron counts give `restricted`, and unequal counts `restricted_open`
(ROHF). Occupations follow the file's orbital order: the first
`Number of beta electrons` orbitals hold two electrons in a single set, the
next ones up to `Number of alpha electrons` one. When the file has a
`Total SCF Density`, the import checks that these occupations reproduce it
within 1e-5 of its largest element and refuses the file otherwise (an SCF
that did not occupy the lowest orbitals). The report gives the deviation.
Orbital kind is `custom`: the file does not say whether its orbitals are
canonical. A `Spin SCF Density` of a single set has no QVF member, and an
unrestricted `Total SCF Density` without the spin density cannot be split:
both are reported as omitted.

**Overlap.** The import computes the AO overlap from the basis (Obara and
Saika recursion, standard library only) and writes the `overlap` member, so
the validator checks its diagonal against Appendix A.1 (`E_NORMALIZATION`).
`import_fchk(..., overlap=False)` skips it for a very large basis: the
pure Python cost grows with the square of the number of shells.

**Refused on import.** Effective core potentials (`Nuclear charges` below the
atomic numbers, or `ECP-` fields), ghost and dummy atoms, shells away from
their atom, l > 7, and counts, charge, multiplicity or electrons that
disagree. Other numeric fields (gradients, Hessians, post SCF densities, an
`SCF Energy` that differs from `Total Energy`), the title and the job type
are listed in the report and not imported.

**Export** writes one molecular `wavefunction.gto` section with real
`restricted`, `restricted_open` or `unrestricted` orbitals and no ECP.
Because fchk has no occupations, the occupations must be integral and in
aufbau order in the row order (2, then 1, then 0 for a single set, 1 then 0
per spin channel), and both channels need the same number of orbitals.
Charge and multiplicity are derived from them and checked against
`provenance` and `basis.n_electrons`. Shells are grouped by atom (stable
order), as Gaussian writes them and as cclib expects. Each shell type carries
its own `pure` flag, with one convention for d shells and one for f and
higher. Contractions are written normalized and the MO coefficients
compensate, so orbital values are preserved. The SCF densities come from the
archive's density members or, without them, from the orbitals and
occupations (`density_source` in the report). A `restricted_open` spin
density is built from the singly occupied orbitals. `Total Energy` comes from
root `energies.total`, `SCF Energy` from `provenance.scf_energy` and
`Mulliken Charges` from a single `atom_properties` section with
`mulliken_charge` for the structure. Array values use Gaussian's E16.8 fields
(9 significant digits), scalars keep full precision. The second line names
the method and basis from `energies` or `provenance`, or `unknown`. Symmetry
labels, orbital kind, the overlap, identities, bonds and jobs are not
retained, and the report lists the omitted sections and members.

**Independent checks.** Gaussian was not run and no Gaussian output files
are used. `python/tests/fchk_fixtures/generate.py` writes synthetic files
from the documented layout without the converter: RHF with SP and pure d, f
and g shells, UHF with Cartesian d, f and g shells, and ROHF with pure d and
Cartesian f (`5D 10F`). The tests compare orbital values at random points
with a direct evaluation of Gaussian's documented AOs, and the MO overlap
C S Cᵀ with IOData's own reading and overlap of the same file. PySCF 2.14.0
RHF, UHF and ROHF calculations, spherical and Cartesian up to g shells, are
exported and read back by cclib 1.8.1 (s to f shells, cclib's limit), by
IOData 1.0.1 (orthonormal orbitals in its overlap) and by the importer,
whose orbitals match PySCF's AO evaluation and satisfy C S Cᵀ = 1 with the
archive overlap, which agrees with PySCF's. Every archive validates with no
findings. Remaining: ECPs, natural orbitals and post SCF densities.

## Python API and limits

```python
from qvf.convert import import_xyz, export_xyz, import_cube, export_cube
# Since toolkit 0.6.0:
from qvf.convert import import_xsf, export_xsf, import_pdb, export_pdb
# In the development checkout (Unreleased):
from qvf.convert import import_molden, export_molden, import_cif, export_cif
from qvf.convert import import_poscar, export_poscar, import_chgcar, export_chgcar
from qvf.convert import import_mmcif, export_mmcif

report = import_xyz("molecule.xyz", "molecule.qvf")
report = export_xyz("molecule.qvf", "geometry.xyz", extended=True)
report = import_cube("field.cube", "field.qvf", kind="volume.generic")
report = export_cube("field.qvf", "field-export.cube", section="volume")
report = import_molden("orbitals.molden", "orbitals.qvf", spin="restricted")
report = import_molden("orca.molden", "orca.qvf", spin="restricted", dialect="orca")
```

Import text is capped at 64 MiB. Geometry/grid importers limit coordinate/scalar
values to five million; Molden limits coefficient matrices to five million
entries and contractions to 512 primitives each.
These initial converters parse bounded inputs in memory; they are not a
replacement for the streaming QVF reader/writer for large datasets. QVF
reader limits also apply on export. Remaining WP-48 work: Molden h shells,
ECP sections and natural-orbital occupations outside 0 to 2, XSF animation, bonds to symmetry mates, centrosymmetric
space groups in the PDB expansion, mmCIF assemblies and writing `REMARK 350`
back. No QVF schema, canonical kind or finding code changes.
