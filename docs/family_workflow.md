# A QVF job across the vibe family

Use one `.qvf` archive to carry a calculation request, send it to a managed
vibe-qc runtime with vibe-scheduler, and open the resulting structure and supported
results in vibe-viewer. A sketch document is another use of the same format;
its export is not automatically a runnable calculation.

This workflow uses the QVF **v0.6.0 toolkit** and **format version 1**. The
[specification](../spec/qvf-format-spec.md) defines the archive contract; each
consumer implements a subset. Conformance does not imply that every consumer
can execute, display or preserve every section.

The [reaction workflow profile](reaction_workflow.md) defines the continuous
scientific interchange requirements, explicit energy selection and immutable
parent hashes. It distinguishes current source contracts from released and
scientifically demonstrated capability.

The frozen QC 0.18, Sketch 0.6, View 3.0 and Queue 1.0 interfaces use the
**v0.6.0 tag's bytes**, including its schema, rather than the current toolkit
development tree or an optional dependency's minimum version. Keep those pins
when qualifying these releases. Later toolkit additions need their own
compatible release and explicit consumer adoption.

The optional `prepare-job` step below requires toolkit v0.6.0. The existing
writer example also works with v0.5.0.

## Check the packaged QVF runtime

For each packaged platform, match the bundled QVF version and source-archive
SHA256 to the chosen release, alongside the source and CI pin. Sketch 0.6
requires QVF 0.6: a package bundling 0.5 cannot save a native QVF document even
when its source checkout tests with 0.6. Updating a cache does not change the
version selected by its runtime lock.

Use the consumer's existing installed-runtime check to save and reopen a native
QVF document and compare the complete document, including drawing, species,
geometry and recipe. Record the actual package and runtime identities with that
result. Reuse accepted checks for unchanged inputs; qualify a corrected package
through its existing release procedure. Keep archive validity, installed
compatibility and scientific acceptance recorded separately.

## Install the relevant components

The Python distribution remains `qvf-writer`. It supplies the `qvf` package
and compatible `qvf_writer` and `qvf_reader` imports. From a QVF checkout
at the v0.6.0 tag, install
into the Python environment that will run the example:

```sh
python -m pip install './python[numpy,validate]'
```

The example uses the legacy NumPy array API. The `validate` extra adds JSON Schema checks to the
reference validator. This source installation does not depend on the release
already being available on a package index or mirror. Use your normal managed
installation for vibe-qc, and install/configure vibe-scheduler on the submission
host. The queue's `--program NAME` must name a registered vibe-qc runtime;
`NAME` and `HOST` below are placeholders for that configuration.

vibe-viewer's core installation can inspect QVF metadata and structures. The
interactive `vibe-view open` command also needs its `viewer` extra. Sketch
import needs vibe-sketcher's chemistry dependencies and `qvf-writer` in its own
Python environment; its installer supports an explicit QVF checkout path.
These applications need not share a Python environment.

The toolkit supplies the unified **`qvf`** command and retains
**`qvf-validate`**. Use `qvf info job.qvf` to inspect an archive and
`qvf validate job.qvf` for reference validation.

## Write and inspect a pending job

Save
[write_pending_job.py](../python/examples/write_pending_job.py) in your working
directory, or use its `python/examples/` path in a checkout that includes it:

```sh
python write_pending_job.py job.qvf
qvf-validate --json job.qvf
vibe-view info job.qvf
```

[The example](../python/examples/write_pending_job.py) writes an illustrative
H2 geometry in angstrom and requests an RHF/STO-3G single point. It performs
no calculation and refuses to overwrite an existing output. The archive has
exactly one `structure`, one `job.spec` referring to it, and
`provenance.run_status = "pending"`. The explicit `functional: null` says
Hartree-Fock takes no density functional; omitting the field would leave the
choice to the engine. There are no calculated energies or orbitals yet.

For a read-only check in the vibe-qc environment:

```python
from vibeqc.qvf_job import load_job_container

job = load_job_container("job.qvf")
assert job.run_status == "pending"
assert job.spec["method"] == "rhf"
assert len(job.structure["atoms"]) == 2
```

`load_job_container` reads the request and verifies the input-member hashes;
it does not run the job and is not a substitute for the reference validator.
It imports the installed vibe-qc runtime, so a missing or incompatible native
library can prevent this read-only check. Restore a usable managed runtime
before claiming the QC reader check passed; `qvf validate` remains available
independently.
The development vibe-qc runner requires **exactly one job.spec**. It resolves
that section's `structure_ref`, or the payload's reference when the section
has none. A job with neither reference requires exactly one structure;
otherwise the runner refuses the ambiguity. Its supported task names are
`single_point`, `optimize` and `hessian`. Older runtimes require exactly one
structure even when a reference is present. QVF itself permits job sets;
executing several jobs remains outside this runner's current input surface.

## Submit, follow and fetch

These commands execute the requested calculation through the configured
runtime. Record the job ID printed by submission:

```sh
vq submit HOST job.qvf --program NAME
vq status HOST JOBID --json
vq fetch HOST JOBID --name job.qvf -o results/
```

vibe-scheduler recognizes the single-file payload and selects the managed
`vibeqc run` entry point. Do not pass `--python` or `--branch` for this path.
Its request parser recognizes the payload shape; it does not validate the
archive's scientific request. On scheduler hosts, the development queue uses
the launcher resolved from the verified target runtime pin. It refuses a
missing launcher or a conflicting configured wrapper before staging input;
no extra wrapper is needed when the runtime registry supplies that launcher.
The optional `--vibeqc-preflight` executes the
configured command and is not a read-only archive validation command.
The `job.spec` is declarative data; neither the queue nor this workflow
executes a stored `run.record.input` as a script. The producer updates the
working archive with results and a sequenced `run.record`. A runtime pin may
be required by site policy; `--expected-sha FULL_SHA` requests an explicit
source identity at submission and dispatch.

Read `qvf_lifecycle` in JSON status. `done` requires both a terminal queue job
and a terminal QVF status with a complete, sequenced input/log record.
`chemistry_failed` distinguishes a recorded failed calculation from
`queue_process_failed`. A queue process exiting zero alone is insufficient.

Development QC runtimes preserve an explicitly referenced input structure,
the original job specification and drawing members when adding results. The
result geometry is separate and selected through `viewer_defaults`; View
can display it while Sketch reopens the original canvas. Earlier QC runtimes
may discard the drawing during finalization, so use a runtime containing the
container-preservation fix for this round trip. A legacy job without a
structure reference retains one geometry and stores the submitted geometry
in its run record.

Fetch after completion; artifact-only fetch preserves archive bytes and
excludes queue sidecars. A retry accepts an identical file and refuses to
overwrite a different result.

```sh
qvf-validate --json results/job.qvf
vibe-view open results/job.qvf
```

A clean validation report establishes the checks performed, not scientific
accuracy or the availability of every visualization. For a deliberate rerun,
submit the result with `--qvf-force`; this appends another sequenced run record.
Use a fresh fetch directory for the new result:

```sh
vq submit HOST results/job.qvf --program NAME --qvf-force
vq fetch HOST NEW_JOBID --name job.qvf -o rerun-results/
```

## Correlated results and campaign acceptance

A correlated calculation can retain `provenance.method: "rhf"` and
`provenance.scf_energy` for its reference while the executed `job.spec` requests
`dlpno-ccsd(t)`. The reference SCF verdict and convergence describe that
reference. They do not establish the correlated method's convergence or total
energy. Preserve these fields and the original archive hash.

The Sketch 0.6 campaign stage reader supports the QC producer's existing
`run.record` attachment layout: `attachment.structured` carries the
program-specific `dlpno_ccsd_done` event. It checks the unique latest successful
run, its `job_ref`, identical hash-verified job/input bytes, the structured
member's hash, a single converged matching event, the requested triples mode
and a finite `e_total` in hartree. A conflicting root total or malformed
structured member is refused. The retained event describes the actual triples
algorithm; a requested mode label alone does not establish that algorithm.
This attachment is existing program output, not a new canonical energy kind.

Reading a total does not accept the calculation. Preserve the original raw
reference rejection and its hash. A permitted explicit operator review lives
in a separate immutable continuation with its evidence; it must not change
the original verdict, computation, protocol or artifact identities. The
campaign path consumes verified accepted stage artifacts. The generic
document `import-result` path instead requires a recorded job/run submission;
it is not a substitute for campaign import or reporting. An incomplete or
unreviewed campaign remains ineligible for an accepted reaction-energy report.

The frozen Sketch campaign export keeps its report in each reaction step's
`x_vibesketch.energies` extension. That layout is allowed by v0.6.0; the
canonical step `energies` field on later toolkit main is outside this pin.

Keep source identity at each layer: the executed runtime's Git revision,
campaign/source-directory digest, artifact manifest and result-file digest,
member digests, and the parent geometry/result digests. A correlation event's
`source_identity` is not the runtime's Git revision. An input structure and a
result structure retain distinct IDs and hashes even when their coordinates
agree to numerical precision. Resolve the job's `structure_ref` and the
result's `run_ref`; do not replace either with whichever structure is last.

View 3.0's normal energy summaries read the SCF/reference energy beside its
provenance method; they do not select a correlated total from a structured
attachment. Queue 1.0 verifies transport and lifecycle identities, not
scientific acceptance. Its released campaign execution also does not include
the later opt-in retained-runtime recovery mapping. Keep these limits separate
from archive validity and from later consumer fixes. None requires changing
the QVF v1 schema or rewriting accepted artifacts.

<a id="open-the-geometry-in-vibe-sketch"></a>

## Open the geometry in vibe-sketcher

A foreign QVF contributes molecular structures and supported reaction-scheme
content to a sketch. Keep the original result archive and write a separate
sketch archive:

```python
from pathlib import Path
from vibesketch.document.qvf import from_qvf, to_qvf

imported = from_qvf(Path("results/job.qvf").read_bytes())
for warning in imported.warnings:
    print(warning)
Path("drawing.qvf").write_bytes(to_qvf(imported.document, preview=False))
```

This is a geometry/document conversion, **not a lossless archive update**.
Foreign job specifications, run logs, wavefunctions, other result sections and
vendor sections may be absent from the new sketch archive. An existing
`x_vibesketch.document` is the sketch's own document representation; generic
QVF readers can still read the canonical structures it exports.

A sketch document may contain several species and jobs, and `to_qvf` does not
mark it as pending. Its queue submission API currently emits engine input
files. Do not send a general sketch document directly to the single-job
vibe-qc runner. Select one molecular geometry and construct an explicit
pending job with the writer pattern above, the helper below, or use the sketch's supported
engine-input submission path. Periodic structures are outside sketch's
molecular drawing scope.

### Preserve a one-species drawing when preparing a job (toolkit 0.6.0)

For a sketch export with exactly one structure and no existing jobs, or a
geometry-only converter output, the toolkit provides an explicit
preparation step:

```sh
qvf prepare-job drawing.qvf --output job.qvf --job-type molecular \
  --method rhf --task single_point --basis sto-3g --no-functional
qvf validate job.qvf --json
```

For an executable starting point, run
[write_sketch_drawing.py](../python/examples/write_sketch_drawing.py) in a
Sketch environment with its chemistry dependencies and the development QVF
toolkit installed, then use the commands above:

```sh
python python/examples/write_sketch_drawing.py drawing.qvf
```

The example exports a real one-species Sketch document with a 2D H–H drawing,
an explicitly supplied illustrative 3D geometry, bonds and the Sketch vendor
document. It creates no job or preview and refuses to overwrite a file.
The `prepare-job` step adds the explicit scientific request. Canvas positions
are drawing coordinates; the separate conformer supplies angstrom coordinates
for the calculation.

Install the development toolkit's `validate` extra. The helper needs no NumPy.
It preserves the original structure, bonds, previews, vendor document bytes
and declared hashes, adds one job linked to that structure, and publishes a
new archive only after full validation. It does not run or submit a calculation.
Use the normal `vq submit` step above with a runtime that supports the chosen
method. Retaining a drawing does not require a result viewer to render it or
certify that its displayed coordinates describe a later optimized structure.
Opening this prepared archive in Sketch restores the preserved original
document. The added canonical `job.spec` does not edit that document's own job
list. Keep the prepared archive when saving a new drawing, because a Sketch
re-export can omit the request or later calculation results.

Specify `--charge` and `--multiplicity` for a charged species or radical;
their job defaults are 0 and 1, and conflicting structural identity is refused.
Omitted basis/functional leave the choice to the engine, while `--no-basis`
and `--no-functional` explicitly say the method takes none. Multiple structures,
existing jobs/runs/results and unsupported source content are refused without
dropping sections. See [command details and refusal rules](cli.md#prepare-one-pending-job-toolkit-060).

The example's reference validation, vibe-qc request load, vibe-viewer structure
read and sketch geometry round trip are format and integration checks. They
make no claim that a quantum-chemistry calculation has been validated.

## Importing other formats with toolkit 0.6.0

The [text converters](converters.md) bring XYZ/extXYZ structures and trajectories,
scalar cube fields, static XSF geometry/grids and single-model PDB geometry into
canonical QVF sections. They
create new, validated archives and report conversion limits. Imported
geometry has no job request: use the explicit pending-job pattern above or
the `prepare-job` helper
before submitting it to vibe-qc through vibe-scheduler. Keep result archives
when exporting a selected geometry or grid to a text format.
