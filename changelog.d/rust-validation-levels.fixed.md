The Rust validator stages `Level` as the Python and C++ reference
validators do: `Container` now reports a missing `manifest.json`, and
`Semantics` runs the cross-section checks (atom counts, orbital rows,
volume spans) and leaves out only `E_REF_UNRESOLVED`. It previously
followed the JavaScript staging, which runs them only from `References`.
Results at the default level `All` are unchanged.
