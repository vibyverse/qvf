# QVF media type registration draft

Draft for maintainer review, 2026-10-15. The proposed registration covers
`application/vnd.vibe-qc.qvf+zip`, the media type already specified for QVF
archives. Submission and the public contact details require the maintainer's
approval. The type remains unregistered until IANA accepts the application.

## Registration fields

| Field | Proposed value |
|---|---|
| Type | `application` |
| Registration tree | Vendor (`vnd.`) |
| Subtype | `vnd.vibe-qc.qvf+zip` |
| Required parameters | None |
| Optional parameters | None |
| Encoding | Binary |
| Intended usage | Common |
| Restrictions on usage | None |
| Provisional registration | No |
| Deprecated aliases | None |
| File extension | `.qvf` |
| Macintosh file type codes | None; the Uniform Type Identifier is `com.vibe-qc.qvf` |
| Object identifiers | None |
| Contact | QVF maintainer, `mpei@vibe-qc.com` (proposed public contact) |
| Author and change controller | The vibe-qc project, steward of QVF |

For the online form, select the vendor prefix separately and enter
`vibe-qc.qvf+zip` in the subtype text box. The form's submitter name should be
supplied by the maintainer. These fields follow the
[IANA application](https://www.iana.org/form/media-types); the vendor registration
process is defined by [RFC 6838](https://www.rfc-editor.org/rfc/rfc6838.html).

## Encoding and identification

QVF is a ZIP container containing a UTF-8 JSON manifest and JSON, binary
numerical, and optional embedded-document members. ZIP STORE and DEFLATE and
ZIP64 are supported. Binary transport encoding is required where the transport
cannot carry arbitrary octets. The `+zip` suffix follows
[RFC 6839 section 3.6](https://www.rfc-editor.org/rfc/rfc6839.html#section-3.6).

There is no QVF-specific leading magic number. A ZIP local-file header usually
starts with hexadecimal `50 4B 03 04`, which identifies ZIP, not QVF. A reader
identifies QVF by inspecting `manifest.json` for `qvf_version` and validating
the manifest and referenced members. The manifest is written last, and an
updated archive can contain superseded manifests; the final central-directory
entry named `manifest.json` is authoritative.

## Published specification

The specification is the QVF format specification, including its JSON Schema
and registry, maintained at [the QVF repository](https://github.com/vibyverse/qvf).
The intended immutable reference for this registration is
[the toolkit 0.3.0 specification](https://github.com/vibyverse/qvf/blob/v0.3.0/spec/qvf-format-spec.md),
with its [schema](https://github.com/vibyverse/qvf/blob/v0.3.0/spec/qvf_manifest.schema.json).
The maintainer should submit once those public tag links are available, or use
a later published release incorporating the same media type.

## Applications and interoperability

QVF carries quantum-chemistry structures, calculations, wavefunctions,
volumetric fields, spectra, trajectories, and related analysis data. The
vibe-qc calculation program, vibe-viewer viewer, vibe-sketcher structure editor,
and independent producers and consumers can exchange these archives.

The manifest carries an integer format version; optional additions within
version 1 do not change the media type. Readers interpret the withdrawn
version 2 as version 1 with a warning. Unsupported higher versions must be
refused. Consumers may ignore unsupported optional sections, but must refuse
an unsupported section marked critical. Numerical members use the declared
shapes, data types, fixed units, and little-endian byte order of the
specification. A generic ZIP utility can inspect the container but cannot
establish chemical or schema conformance.

## Fragment identifiers

QVF defines no fragment identifier syntax. Generic fragment handling follows
the registered `+zip` suffix; this registration adds no QVF-specific semantics.
Section identifiers in the manifest are internal references, not URI fragments.

## Security considerations

Treat every archive as untrusted. ZIP processing requires bounds on expanded
member and total sizes, JSON sizes and nesting, numerical shape arithmetic,
and processing resources. The QVF validator limits are specified in section
6.1; consumers may impose smaller resource budgets. Declared sizes must be
checked against actual expanded bytes. Compression ratio alone is not a
format validity limit: a large constant array can compress very well. Reject
invalid archive paths, duplicate member names other than superseded manifests,
corrupt compressed data, size mismatches, and checksum failures. Do not extract
members outside a selected destination or follow archive-controlled links.

The security considerations of ZIP, JSON, and any embedded image or document
format also apply. Render previews with hardened decoders and bounded image
allocation. SVG previews must be self-contained and inert; consumers must
disable scripts and external resource retrieval even if a supplied document
violates those rules. A schema-valid archive does not prove its embedded SVG
is safe to render.

An archive may contain arbitrary scripts as historical `run.record.input`
attachments. They are records, not instructions to execute, and must never be
run as a side effect of opening or executing an archive. A `job.spec` is
declarative calculation input; a runner must separately authorize computation
and enforce its resource and execution policy. Opening a document must not
silently authorize a calculation. Vendor sections and opaque attachments
require the same trust boundaries.

Structures, input files, logs, and provenance can reveal confidential research
or personal information. QVF provides no confidentiality or authenticated
sender identity. Member SHA-256 values and ZIP CRC-32 detect mismatches but
cannot authenticate an archive whose manifest an attacker can also replace.
Use authenticated transport, access controls, and external encryption or
signatures when confidentiality or origin authentication is needed. Appending
an update leaves superseded content inside the ZIP; removing a manifest
reference does not erase that content. Repack and inspect the resulting archive
before sharing data intended to exclude earlier content.

Interpretation of the core archive does not require retrieval of remote
resources. Schema identifiers and provenance links do not authorize automatic
network access. Consumers should use a trusted local schema and separately
control any explicit external retrieval.

These considerations describe handling of the existing format; they do not
relax its validation rules or add an execution mechanism. See
[the QVF security contact](../SECURITY.md) for vulnerability reports.
