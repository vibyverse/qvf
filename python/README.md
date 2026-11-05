# qvf-writer (Python)

The `qvf` package writes, reads, validates and updates **QVF** (`.qvf`, the
*Quantum Visualization Format*), a ZIP container for quantum-chemistry data
with a typed, checksummed manifest. The distribution name stays `qvf-writer`.

The package and its raw-buffer writer, reader and editor use only the standard
library. Install optional dependencies for numeric convenience writers and
full JSON Schema validation. Python 3.9+ is supported by the base package.

Toolkit 0.7.0 also includes `qvf export-xtb` and
`qvf.adapters.xtb.write_xtb`: package molecular GFN0/1/2-xTB JSON, geometry
and explicitly selected trajectory, bond-order, Molden, normal-mode and
cube outputs as one validated QVF. See the
[xTB exporter guide](https://github.com/vibyverse/qvf/blob/main/docs/adapters/xtb.md).
This adapter needs neither xTB nor NumPy at export time.

## Python release candidate

`0.7.0` is the first PyPI release candidate, the same version as toolkit
0.7.0, and includes the PySCF adapter. It has not been uploaded yet. Until
then, install the tagged source:

```sh
python -m pip install "qvf-writer[pyscf,validate] @ git+https://github.com/vibyverse/qvf.git@v0.7.0#subdirectory=python"
```

The following index command works **after the package is uploaded**:

```sh
python -m pip install "qvf-writer[pyscf,validate]==0.7.0"
```

For the PySCF tutorial use Python 3.11–3.13 on a platform with a PySCF wheel
(Linux or macOS). The QVF base package does not require PySCF or NumPy.
To test a built wheel before publication, install it directly:

```sh
python -m pip install "./qvf_writer-0.7.0-py3-none-any.whl[pyscf,validate]"
```

Toolkit 0.7.0 does not change the archive format version or move consumers
pinned to toolkit v0.6.0. That older tag has no PySCF adapter. To match
their pin, install its tagged source:

```sh
python -m pip install "qvf-writer[numpy,validate] @ git+https://github.com/vibyverse/qvf.git@v0.6.0#subdirectory=python"
```

`qvf_writer` and `qvf_reader` remain importable with a `DeprecationWarning`;
existing `QvfWriter`, `QvfReader` and `qvf-validate` callers keep working.

```python
import struct
from qvf import Reader, Section, Writer

with Writer("density.qvf", program="example", version="1.0") as writer:
    writer.add_section(Section(id="rho", kind="volume.density", members={}))
    writer.add_json("rho", "grid", {
        "origin": [0, 0, 0], "voxel_vectors": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "shape": [1, 1, 2],
    })
    with writer.begin_binary("rho", "data", "float64", [1, 1, 2]) as member:
        member.write(struct.pack("<2d", 0.5, 0.25))

with Reader("density.qvf") as reader:
    tensor = reader.tensor("rho", "data")
    print(tensor.tolist())                # [[[0.5, 0.25]]]
    raw = tensor.data                     # read-only byte memoryview
    assert reader.validate()["ok"]
```

Members are written as they are added, and the manifest is written last.
`Reader` checks resource limits at open and loads payloads on demand.
`Editor` appends updates, optionally commits atomically, and repacks orphans.
`qvf.validate.run` shares finding codes with the C++ and JavaScript validators.
The installed package contains the normative schema for structural checks.

Start with the runnable
[PySCF-to-QVF tutorial](https://github.com/vibyverse/qvf/blob/main/docs/tutorials/pyscf-to-qvf.md).
The optional adapter is included in this candidate; see the
[adapter guide](https://github.com/vibyverse/qvf/blob/main/docs/adapters/pyscf.md)
for supported methods and limitations.

```sh
qvf info density.qvf
qvf validate density.qvf
qvf cat density.qvf rho grid
qvf-validate density.qvf  # compatibility entry point
```

The `qvf` command also provides `extract`, numeric `diff`, `repack`, explicit-output
`scrub`, `conformance` and `schema`. `python -m qvf --help` lists the commands.
Validation has cumulative container, manifest, payload, semantic and reference
stages; full checks are the default. Output commands refuse existing files.
Scrub offers a dry-run change report and keeps the original archive intact.

See the [format specification, integration guide, and library
guide](https://vibe-qc.com/docs/qvf/). The C++ library and language-independent
conformance suite live in the same toolkit. Licensed Apache-2.0.
