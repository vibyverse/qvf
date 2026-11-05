# Python API reference

This page is generated from the docstrings of the `qvf` package in
[`python/qvf/`](../../python/qvf/) by
[mkdocstrings](https://mkdocstrings.github.io/). The sources are read
statically and the package is not imported. The
[library guide](../library_guide.md#python-package-api) explains how the
pieces fit together, and the [specification](../../spec/qvf-format-spec.md) is
normative where a docstring is shorter.

`import qvf` exposes `Writer`, `Reader`, `Editor`, `Tensor`, `Limits`,
`Manifest`, `Section`, `MemberSpec`, `Source`, `Lattice`, `GridSpec`,
`validate`, `validate_qvf`, `FINDING_CODES`, `QvfWriter`, `QvfReader` and
`__version__`. They are documented below under the module that defines them.
Names starting with an underscore are internal and not listed.

The optional PySCF adapter, `qvf.adapters.pyscf` (the `pyscf` extra), is
documented after the package modules; the [PySCF adapter guide](../adapters/pyscf.md)
shows it in use. It is part of the 0.7.0 package and is not part of the v0.6.0 package.

::: qvf.writer

::: qvf.reader

::: qvf.kinds

::: qvf.editor

::: qvf.model

::: qvf.validate

::: qvf.units

::: qvf.convert

::: qvf.jobs

::: qvf.operations

::: qvf.conformance

::: qvf.adapters.pyscf

The standard-library-only xTB companion exporter is described in the
[xTB guide](../adapters/xtb.md). It is available from the development source.

::: qvf.adapters.xtb

::: qvf.adapters.xtb_run

## Compatibility writer

`QvfWriter` is the original writer. It remains supported and is re-exported as
`qvf.QvfWriter`, and through the `qvf_writer` compatibility module.

::: qvf._legacy_writer.QvfWriter
