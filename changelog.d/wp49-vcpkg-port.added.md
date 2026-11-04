Add a vcpkg overlay port in `packaging/vcpkg` (WP-49). It builds the C/C++
library from the `v0.6.0` tag archive (SHA512 pinned) or from a local QVF
source tree selected with `QVF_VCPKG_SOURCE_DIR`, for static and shared
triplets, and installs the CMake package, `qvf.pc`, the single header and
the license notices. Nothing is published to a registry.
The `@qvf/core` tarball from `npm pack` now contains the compiled `dist/`
build, LICENSE and NOTICE, and no tests or scripts (`files` and `prepack`).
The Python sdist no longer ships `python/tests`, which only run from the
repository checkout, and the wheel lists its supported Python versions.
`docs/python-release.md` adds a trusted-publisher checklist for maintainers.
