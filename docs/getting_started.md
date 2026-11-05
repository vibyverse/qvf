# Getting started with QVF

Use QVF to carry geometry, calculation requests and results in one archive.
Start with a small archive you can read back and validate, then add the kinds
your producer and consumer support. No quantum-chemistry engine is required
to use the toolkit.

## Install a known toolkit

The latest release is **QVF toolkit v0.7.0**, using **format version 1**.
Install its tagged source with Python 3.9 or later and Git:

```sh
python -m venv .venv
# POSIX shells; on Windows use .venv\Scripts\Activate.ps1 in PowerShell.
. .venv/bin/activate
python -m pip install "qvf-writer[numpy,validate] @ git+https://github.com/vibyverse/qvf.git@v0.7.0#subdirectory=python"
python -m qvf --version
```

The released baseline for QC 0.18, Sketch 0.6, View 3.0 and Queue 1.0 remains
**v0.6.0**. To match it, use `@v0.6.0` in the same command. Adopting 0.7.0 in
those consumers is a separate step.

PyPI publication is still planned. The distribution name is `qvf-writer`,
while the import and command name are `qvf`. The core reader, editor and
raw-buffer writer use the standard library; `numpy` adds numerical writer
helpers and array views, and `validate` adds JSON Schema checks.

If you already have a checkout, install its contents from the repository root:

```sh
python -m pip install './python[numpy,validate]'
git rev-parse HEAD
python -m qvf --version
```

Do not switch a shared checkout merely to run an example. A development
checkout may report the previous release version while containing newer
code. Record the tag/commit or package hash with the validation result;
see [validator identity](cli.md#choose-and-record-the-validator).

| Need | Use |
|---|---|
| Frozen family interchange, job preparation, standard readers/writers | The consumer's pinned v0.6.0 toolkit |
| PySCF adapter, broader Molden conversion, new reaction-energy or biomolecular fields | Toolkit v0.7.0 and its [release notes](../CHANGELOG.md#070---2026-11-05); adopt in consumers only after they pin a compatible release |
| Build or embed C++, C or Fortran | [Library guide](library_guide.md), [C ABI](c_api.md), [Fortran](fortran.md) |
| Inspect an archive in a browser or JavaScript application | [JavaScript reader](library_guide.md#javascript-reader) |

## Create, read and validate an archive

Save this as `example.py` and run `python example.py`. It writes illustrative
H₂ geometry, with coordinates in Å. It does not run a calculation or invent
result energies.

```python
from qvf import Reader, Writer

with Writer("h2.qvf", program="example", version="1.0",
            calculation="hydrogen geometry") as writer:
    writer.add_structure([
        {"symbol": "H", "atomic_number": 1, "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "atomic_number": 1, "position": [0.0, 0.0, 0.74]},
    ])

with Reader("h2.qvf") as reader:
    report = reader.validate()
    assert report["ok"], report
    atoms = reader.json("structure", "structure")["atoms"]
    assert len(atoms) == 2
    print(atoms)
```

The writer context closes the archive and writes its manifest. Keep the
reader open while using lazy section views. Opening a reader checks the
container; `validate()` also reads and checks the payloads.

```sh
python -m qvf info h2.qvf
python -m qvf cat h2.qvf structure structure
python -m qvf validate h2.qvf --json
```

A **section** has an ID and a kind. Its **members** hold the JSON or binary
payloads. Here `structure` is both the section ID and the role of its JSON
member; those are separate names, not an archive filename.

## Prepare a calculation request

The same geometry can become a pending request without changing its original
members:

```sh
python -m qvf prepare-job h2.qvf --output h2-job.qvf --job-type molecular \
  --method rhf --task single_point --basis sto-3g --no-functional
python -m qvf validate h2-job.qvf --json
```

`--no-functional` records JSON null because Hartree–Fock takes no density
functional. Charge and multiplicity default to 0 and 1; supply them explicitly
for ions or radicals. Preparation does not execute or submit the request.
The [family workflow](family_workflow.md) covers submission, result retrieval,
Sketch drawing preservation, View inspection and Queue transport limits.

The CLI requires a new output path. Preserve accepted source archives and
their hashes; create a derivative file when changing a document. The library
`Writer` creates its target and can replace an existing file, while `Editor`
updates an archive intentionally; see [editing](library_guide.md#python-package-api).

## Add results deliberately

Use the [kind catalogue](kinds.md) to choose existing kinds and the
[integration guide](integration_guide.md) to map your engine's data. Set
explicit structure and wavefunction references when an archive contains
several geometries or orbital sets. Store the actual source program, method,
basis and convergence status; keep an SCF reference separate from a correlated
or excited-state result.

For wavefunctions, check Appendix A's primitive normalization, AO order and
spherical signs. In the QVF coefficient array `[n_mo, n_ao]`, AO indices are
columns. Validation cannot prove the orbitals reproduce an engine's values,
nor can it establish scientific correctness.

| Finding or symptom | Next step |
|---|---|
| `E_SHA256_MISMATCH` | Compare with the original producer artifact; do not repair an accepted archive by silently changing its hash. |
| `E_REF_UNRESOLVED` | Check section IDs, kinds and state indices; add the referenced section. |
| `E_NORMALIZATION` | Check primitive coefficients, AO ordering and the overlap convention in Appendix A. |
| Valid archive, missing display or unsupported job | Check that consumer's supported kinds and operations; format validity alone does not promise support. |
| A legacy validator rejects a field accepted by the producer's pin | Re-run with the pinned validator and record both identities; do not alter source data to fit the older schema. |

See the [command-line guide](cli.md) for exit statuses, resource limits,
comparison and explicit-output archive maintenance. The
[specification](../spec/qvf-format-spec.md) remains the normative contract.
