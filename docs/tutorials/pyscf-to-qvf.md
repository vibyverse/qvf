# From a PySCF calculation to a validated QVF archive

This tutorial runs a small molecular Hartree–Fock calculation, exports its
geometry and wavefunction, and checks the archive using the installed QVF
package. No vibe-qc installation is required. The numerical calculation is
performed by PySCF; QVF stores and validates the result.

## Install in a fresh environment

Use Python 3.11–3.13 on Linux or macOS with a compatible PySCF wheel. The QVF
base reader/writer supports Python 3.9+, but the optional chemistry engine
has its own platform and Python requirements.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
```

The adapter is part of toolkit 0.7.0, and the Python package `qvf-writer` has
the same version. **Before PyPI publication**, install it from the tagged source:

```sh
python -m pip install "qvf-writer[pyscf,validate] @ git+https://github.com/vibyverse/qvf.git@v0.7.0#subdirectory=python"
```

**After publication**, the equivalent command will be:

```sh
python -m pip install "qvf-writer[pyscf,validate]==0.7.0"
```

The explicit version pin makes this example reproducible. Toolkit v0.6.0
predates the adapter; installing that tag does not provide this example.
`pyscf` installs the optional engine and NumPy; `validate` enables full JSON
Schema checks in addition to the always-enabled semantic checks.

## Run the example

Save [pyscf_to_qvf.py](../../python/examples/pyscf_to_qvf.py) in your working
directory (it is also included in the source distribution), then run:

```sh
python pyscf_to_qvf.py h2.qvf
qvf info h2.qvf
qvf validate h2.qvf
qvf-validate h2.qvf
```

The script defines H₂ with a 0.74 Å bond and the STO-3G basis, runs restricted
Hartree–Fock with one thread, refuses an unconverged result, and exports with:

```python
from qvf.adapters.pyscf import write_scf

# mf is the completed PySCF SCF object.
ids = write_scf(mf, "h2.qvf")
```

The full example uses exclusive file creation to protect previous results.
For a second run choose another output filename. Expected output includes:

```text
qvf-writer 0.7.0
RHF/STO-3G total energy: -1.1167593074 hartree
MO coefficient shape: (2, 2) (orbitals, AOs)
Maximum orbital orthonormality error: ...
Validation: PASS (zero findings); wrote h2.qvf
```

Tiny floating-point differences are normal. The energy is illustrative for
this specific method, basis and geometry, not a high-accuracy reference for
H₂. The script compares the stored energy with the calculation to 1e-12
hartree and checks `C @ S @ C.T` against the identity to 1e-10. Coefficient
rows are molecular orbitals; columns are QVF atomic orbitals. Those checks
exercise the archive readback rather than relying on a successful file write.

## Inspect and reuse the result

```python
from qvf import Reader

with Reader("h2.qvf") as reader:
    print(reader.raw_manifest["provenance"])
    coefficients = reader.get("wf").coefficients().numpy()
    metadata = reader.json("wf", "mo_metadata")
    print(metadata)
```

The archive contains geometry, basis, orbital coefficients, energies and
occupations, overlap and density matrices, population analyses, a job
specification and PySCF provenance. It uses the adapter's documented QVF AO
ordering and normalization; raw coefficients should not be compared with
PySCF arrays without that mapping. A QVF-capable viewer can open this archive;
which fields it displays depends on that viewer's version.

For DFT, open shells, excited states and NTOs, see the
[adapter guide](../adapters/pyscf.md). This example covers molecular RHF;
it does not establish support for periodic or complex wavefunctions. Cite
PySCF for the calculation and QVF for the archive format; the repository's
[CITATION.cff](../../CITATION.cff) provides the QVF citation metadata.
