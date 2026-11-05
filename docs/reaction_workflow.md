# Reaction workflow interchange profile

This profile connects editable reaction drawings, pending calculations, accepted
and rejected results, and reaction analysis using **QVF format 1**. It reuses
the [format specification](../spec/qvf-format-spec.md) and existing extensions;
it does not introduce a schema revision. The [single-job procedure](family_workflow.md)
explains submission and retrieval. A supported field is not evidence that a
particular released program implements this entire profile.

## Compatibility and ownership

The family currently pins toolkit **v0.6.0**. The canonical
`reaction.scheme.steps[].energies` contract from QVF #16 is in toolkit 0.7.0;
producers must not write it under the v0.6.0 pin. Existing Sketch reports use
`step.x_vibesketch.energies`. View's reaction analysis uses
`step.x_vibeview.energies`, containing the #16 terms/deltas shape and selection
metadata described below. A reader refuses competing representations on the
same step rather than choosing one silently. Move to the canonical field only
after a containing QVF release and explicit consumer adoption.

This is the **28 October 2026 source integration agreement**, tracked under
QVF #15/#16, QC #702 and View #164 with the Sketch and Queue owners. It is
not a claim that the pending fixes have shipped. The critical-content
preparation fix and audit example are in toolkit 0.7.0, not v0.6.0. Record the actual
containing release and component test evidence when the release owner
qualifies them; do not infer release capability from this document.

The component owners maintain execution and acceptance semantics. QVF owns
the shared format, schema, registry and reference validator. Queue transports
QVF and binds attempts; it does not infer scientific success. These profile
requirements are consumer checks beyond generic QVF validation. Unsupported
profile versions or required information must be refused by the operation
that needs them. An archive may still be opened for read-only inspection with
that limitation stated.

## Information carried at each edge

| Information | Existing QVF carrier | Consumer obligation |
|---|---|---|
| Editable canvas | `x_vibesketch.document` and its original members; declared `x_vibesketch` extension | Preserve drawing bytes and species index through calculation; distinguish original drawing from optimized geometry. |
| Species and geometry | `structure`, chemical identity, `bonds`, `reaction.scheme` participants | Resolve species and structure IDs, explicit atom ordering/maps, roles and positive stoichiometry. Check charge and multiplicity against the authoritative job request. Positions are angstrom, not canvas coordinates. |
| Preparation | Sketch document/stage members and run attachments | Retain ETKDG version, seed, attempts, selected conformer, MMFF termination and warnings. A failed preparation cannot become an accepted parent. |
| Pending stage | `job.spec` with explicit `structure_ref`, `options`, `provenance.run_status=pending`; Sketch stage member | Carry requested and resolved settings, stage/species identity, accepted parent result and geometry identity. Independent `job.spec` sections are not themselves a dependency graph. |
| Execution | `run.record` with `job_ref`, sequence, status, input, log, attachments | Preserve the request verbatim; record actual package/source/native identity and runtime settings. Queue runtime selection is not proof of the runtime that executed. |
| Results | Separate result `structure`, `trajectory`, `forces`, `hessian`, `vibrations`, root `energies`; section `run_ref`/`structure_ref` | Retain full trajectories, same-surface Hessian settings and convergence evidence. Root energy selection must resolve the current result run and geometry. |
| Acceptance or rejection | Hash-bound Sketch stage acceptance or View selection review; original run verdict and attachments | Preserve rejected evidence and reasons. Acceptance is an explicit operation, not a translation of exit code or archive validity. |
| Reaction analysis | `reaction.scheme` and versioned energies extension below | Select one property from each accepted source, check method/state/geometry links, balance and coefficients; retain sources and reproduce the sum. |

For the molecular demonstration, the requested sequence is ETKDGv3/MMFF94,
MACE-OFF23 medium optimization, PBE/def2-SVP optimization without dispersion,
same-surface Hessian, and TightPNO DLPNO-CCSD(T)/cc-pVDZ endpoint energy.
Record actual runtime identities; final paper chemistry uses a released
vibe-qc **>=0.18,<0.19**, without requiring a particular patch. Do not alter an
in-flight environment or substitute a new computation for retained failures.

Resolved settings include model checkpoint hash, dtype/device and optimizer
criteria for MACE; explicit functional, basis, dispersion, fitting/COSX,
grid and SCF tolerances for PBE; displacement and stencil for the Hessian;
auxiliary basis, frozen-core policy/count, numerical PNO thresholds,
reference contract and actual triples algorithm for DLPNO. A preset or mode
label alone is insufficient. The Hessian acceptance prerequisite and the
PBE geometry parent for a correlated endpoint are separate dependency edges.

Atom maps must survive independently of normalized SMILES. Stable drawing
atom IDs and mappings belong to the preserved document/atom extensions when
there is no canonical carrier. Reordering atoms requires an explicit verified
permutation, including hydrogens; a matching molecular formula is insufficient.

## Immutable identities and evolving containers

Use distinct identities for distinct objects:

1. **Archive SHA-256** hashes the exact complete ZIP file bytes at a recorded
   handoff. Keep that snapshot immutable. Repacking, accepting, adding a run,
   or editing a drawing creates a different archive identity even if geometry
   is unchanged. A UUID or `metadata.derived_from` is not a content checksum.
2. **Member SHA-256** hashes the exact uncompressed payload bytes named by a
   section/member descriptor. Check it before parsing. Record the owning
   archive, section ID and member role so equal hashes cannot hide a wrong
   species or wrong result link.
3. **Geometry identity** may use a producer's semantic digest, but must name
   that algorithm and serialization, units, atom order/identity and state.
   It is not interchangeable with the structure-member or archive hash.

The QC container profile records `run.record.x_input_archive_sha256` as the
digest of the immutable archive submitted to **that invocation**. A forced
rerun records the preceding result archive's hash in the new run. Keep the
first run and its original pending hash. `attachment.input_manifest` and
`attachment.result_manifest` preserve the submitted and generated manifests;
neither snapshot alone reconstructs the original ZIP bytes. Retain the
original submitted archive as a separate immutable artifact.

An accepted derivative refers to the hash of the original raw result and the
selected structure member. It must not claim that hash identifies itself.
Its own complete-file hash is recorded by its consumer or artifact receipt
after finalization. The next pending archive refers to that accepted
derivative and retains its acceptance evidence inside QVF. JSON geometry or
acceptance sidecars may be caches only: consumers must reconstruct and verify
the scientific information from QVF and reject a conflicting cache.

Sketch's `vibesketch_document` vendor section carries an `acceptance` JSON
member with `schema: "vibe-sketch.stage-acceptance/1"`, `raw_result_sha256`,
the full stage `verdict` (including accepted/reasons/recipe/stage/species),
and `geometry: {section_id, member_path, member_sha256}` or null. The stage
member's `start` keeps `result_sha256` as the original raw result hash and
adds `acceptance_artifact_sha256` for the accepted derivative, plus
`structure_ref` and `member_sha256`. Its existing `geometry_sha256` hashes
the original geometry JSON cache bytes, explicitly labelled
`sha256-file-bytes`; it is not a canonical scientific geometry hash. A
separate `prerequisite` carries the Hessian artifact/raw-result hashes,
stage and acceptance. In an accepted stage directory, `result.qvf` remains
the original raw result and `acceptance.qvf` is the derivative; a rejected
`assessment.qvf` remains inspectable without publishing `acceptance.qvf`.
Consumers verify both archive identities and cannot feed rejected derivatives
to descendants. These are versioned vendor payloads; the generic schema
does not validate their scientific acceptance semantics.

Queue's `vq.campaign/1`, `vq.campaign-artifact/1` and attempt journal remain
control records: node/output links, artifact hashes, runtime and attempt
identities. They complement the QVF scientific parent edges and never replace
them. A failed parent blocks only its descendants. Retrying creates a new
attempt and preserves the original rejected result.

## Energy selection and reaction analysis

SCF reference energy stays at `provenance.scf_energy`, with its reference
`provenance.method`. For a correlated endpoint, `energies.method` names
`dlpno-ccsd(t)` and `energies.total` is `{value, units: "hartree"}`. Components
retain reference, CCSD correlation and triples contributions. QC's
`components.correlation` includes CCSD correlation plus triples;
`x_reference`, `x_ccsd_correlation` and `x_triples` distinguish its pieces.
Do not add those pieces twice or add an unrequested dispersion correction.

The QC profile's open provenance fields `result_structure_ref` and
`result_run_ref` bind current root energies to the generated result. The
container finalizer remaps them to its preserved graph's IDs. A consumer
checks those links, the run's `job_ref`, method/settings, charge,
multiplicity and selected geometry. A current run must not inherit stale
root energy or correlation metadata from a preceding stage. Its retained
manifest attachments keep historical values attributable to their own runs.

QC's `correlation_converged`, `correlation_settings` and reference verdict
carry producer evidence; they do not establish a minimum or authorize an
override of a failed reference. View additionally records an explicit
source-hash-bound review. Missing legacy correlated totals are refused;
neither SCF energy nor an external JSON summary is a fallback.

View's v0.6-compatible `step.x_vibeview.energies` object contains:

- `version: 1`, `deltas.energy.value_eh`, and complete signed `terms` using
  the #16 shape. Terms name participant structures, selected absolute
  `energy_eh`, method/program/version/specification, source archive SHA-256,
  parser and exact property location, and a source `run_ref`.
- `x_vibeview.version: 1`, `quantity: "correlated_total"`, and
  `conversion.hartree_to_kj_mol: 2625.499639216`, with the definition
  `627.509474 kcal/mol per hartree × 4.184 kJ/kcal` recorded explicitly.
- `x_vibeview.selections`, keyed by participant structure ID, with the source
  structure ID, exact JSON property pointer (such as `/energies/total`),
  quantity, method and `review` containing status, source SHA-256 and note.
- A `run.record` attachment containing each original source QVF, permitting
  reopening to recheck its hash, energy pointer, run/geometry linkage,
  settings, structures and derived value without an external sidecar.

Declare `extensions.x_vibeview = {"version": "1", "critical": false}`.
Generic readers may preserve this optional analysis without interpreting it;
an analysis operation must understand its version and reject ambiguity.
All three energy representations (canonical, Sketch and View) are mutually
exclusive on a step for this profile.

For butadiene + ethene → cyclohexene, coefficients are −1, −1, +1. Recompute
the weighted sum in hartree, then apply the recorded conversion constant.
Compare complete method definitions, while retaining species-dependent
frozen-core counts and checking their requested policy. Mixed methods,
missing terms, duplicate selections, wrong state/geometry, mismatched reviews
and nonfinite values are errors. This is electronic ΔE; thermal corrections,
barriers and fictitious connecting trajectories are outside this analysis.

## Reproducible checks and evidence

Record three separate outcomes: structural/checksum validation, tested
producer/consumer integration, and scientific acceptance. Full schema checks
must actually be available; `W_SCHEMA_SKIPPED` is not a full validation pass.
A structurally valid archive may be scientifically rejected or unsupported.

The read-only [audit example](../python/examples/audit_interchange.py) checks
strict validation, exact preservation of every input section declaration and
member byte, extension declarations, and whole-archive identities:

```sh
python python/examples/audit_interchange.py pending.qvf result.qvf > handoff.json
python python/examples/audit_interchange.py result.qvf reopened.qvf > reopen.json
```

Use the actual Sketch producer and queued QC result, inspect/analyze in View,
then reopen and save in Sketch using its supported preservation path. A
geometry-only conversion or a drawing edit intentionally fails the additive
audit: retain the original archive and label the derivative accurately. The
audit never asserts scientific validity, launches work, or modifies inputs.

Reuse `conformance/corpus/reaction_scheme.qvf`, `reaction_scheme_energies.qvf`
and their existing negative counterparts for role/reference/sum/method tests.
The canonical energies fixtures test toolkit 0.7.0, not the v0.6.0 pin.
`test_prepare_job.py` checks drawing preservation, corrupt payloads and
unsupported critical content; `test_interchange_audit.py` checks additive
handoffs and structurally valid lossy saves. Component tests additionally
must cover wrong species/geometry, missing accepted prerequisites, ambiguous
energy selection, changed source hashes and unsupported required settings.

Do not copy publication calculations into the format corpus. Read actual
paper archives in place, record their real paths and hashes in external
evidence, and preserve both accepted and rejected originals. Synthetic
round-trip fixtures exercise interfaces only and must be labelled as such.
Only the paper replay owner launches the final integrated chemistry. The
release owner records which released versions contain the changes; a main
commit or passing synthetic test is not a released scientific demonstration.

## User procedure after component qualification

1. In Sketch, draw or reopen the mapped butadiene/ethene/cyclohexene scheme,
   prepare explicit-hydrogen conformers and review termination, state and atom
   maps. Save the native QVF and the selected staged recipe. Record the
   actual component versions and hashes before launch. Check that their
   released versions implement this profile, including the numerical fixes;
   do not treat the minimum version range alone as qualification.
2. Compile the reviewed recipe into a new campaign directory with Sketch's
   `campaign-compile` command. Use the operator's existing runtime settings
   and model-license acknowledgement where required. Inspect the generated
   pending QVFs, stage graph and settings. Every scientific dependency must
   read the accepted QVF, not silently fall back to a script/XYZ/JSON route.

   ```sh
   vibe-sketch campaign-compile reaction.qvf RECIPE --out campaign
   vq campaign validate campaign/campaign.json
   ```

3. The replay owner launches the compiled campaign using the exact command
   printed by Sketch, normally `vq campaign run campaign/campaign.json
   campaign/run`. Record queue attempts and immutable retrieved artifacts.
   Inspect each stage's QVF acceptance/rejection and both Hessian stencils
   before accepting the endpoint. A terminal queue state does not waive any
   scientific criterion. Preserve failures; retry explicitly with a new
   attempt rather than overwriting the original.
4. In View, use `vibe-view reaction inspect ENDPOINT.qvf` for each endpoint.
   Review the actual geometries, trajectories, modes, convergence, reference
   verdict and method. Prepare `selections.json` with the exact scheme
   section/step and one endpoint per structure ID; each selects its path,
   result structure, `quantity: "correlated_total"`, method,
   `property: "/energies/total"` and explicit hash-bound review. It contains
   selection controls, not replacement scientific results. Then run:

   ```sh
   vibe-view reaction build reaction.qvf selections.json -o analysis.qvf
   vibe-view reaction export analysis.qvf --step STEP \
     --csv reaction.csv --figure endpoints.svg > analysis.json
   vibe-view open analysis.qvf
   ```

5. Reopen the saved analysis in View to recompute from its embedded sources.
   View's analysis is a new derivative retaining whole original QVFs in
   attachments, not an additive update of their top-level sections. Check
   embedded-source hashes with the consumer's replay operation. Reopen the
   original native drawing in Sketch for editing; result/analysis imports
   may offer geometry conversion only. The supported Sketch save guard
   refuses overwrite when it cannot retain foreign data and requires Save As.
   Do not call that conversion a lossless analysis round trip.

Use `audit_interchange.py` on additive pending-to-result updates. Record
archive and member hashes, exact source/runtime identities, commands,
limitations and the containing releases alongside the scientific review.
The commands above are a reproducible procedure, not a record that the
complete publication replay has already run.
