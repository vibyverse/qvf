# Contributing to QVF

QVF development, issues and pull requests are hosted at
[github.com/vibyverse/qvf](https://github.com/vibyverse/qvf).
Keep patches focused, explain the expected behavior, and include the relevant
source tests or conformance cases. Contributions use this repository's
**Apache-2.0** license; see [LICENSE](LICENSE).

Keep real machine configuration, credentials, account names and private
data outside the source tree. Enable the staged-content privacy hook once
per clone:

```sh
git config --local core.hooksPath .githooks
python3 .githooks/test_privacy_hook.py -q
```

Use documented placeholders for paths and RFC 5737 address ranges for
examples. The hook checks staged text, and does not replace the filename,
binary, secret and full-publication scans required before publishing a
snapshot. Report security and privacy issues privately per [SECURITY.md](SECURITY.md).

### Continuous integration

Ordinary development pushes and merges use `[skip ci]` and focused local
source checks. Full CI belongs to authorized release candidates. The
[GitHub workflow](.github/workflows/ci.yml) runs on `release-candidate/*`
pushes or an explicit manual dispatch on a candidate branch. Main pushes
and pull requests do not start the full matrix. Do not pass GitLab's
`ci.skip` push option to GitHub.

Candidate CI checks:

- **Python tests** on the supported endpoints (3.9 and 3.14 on GitHub),
  including the privacy hook self-test. The suite already
  runs the reference conformance runner over positive, negative and legacy
  archives; there is no duplicate standalone conformance job. The 3.14 job
  installs PySCF 2.14.0 to exercise the optional
  adapter tests, with OpenMP and BLAS limited to one thread. The 3.9 job
  checks the package without PySCF.
- **JavaScript tests** on Node 20 and 24 on GitHub.
  `npm test` builds and checks the source once; CI separately checks the script
  and test TypeScript project.
- **C++ build and tests** on Linux (gcc and clang), macOS (clang) and Windows
  (MSVC), with warnings as errors: `-Wall -Wextra -Werror`, or `/W4 /WX` on MSVC,
  including the C++ conformance runner over the positive, negative and legacy
  corpora in CTest. Do not rerun that executable after CTest.
- **Sanitizers**: the C++ tests (the conformance runner included) under
  AddressSanitizer and UndefinedBehaviorSanitizer (clang), then the C++ reader
  over every `conformance/corpus/*.qvf` archive. A sanitizer finding fails the
  job.
- **Single-header amalgamation**: `cpp/qvf_single.hpp` is generated and must
  match its sources.
- **Generated tables**: the kind, unit and finding-code tables and the
  compiled JSON Schema (C++, TypeScript, Python), and the kind catalogue
  `docs/kinds.md`, are generated from the schema by `scripts/gen_kinds.py`
  and must match it.

For an authorized candidate, these are the equivalent checks from the
repository root. Run native commands under your native-build lock; keep
operational logs outside the source tree.

```sh
# Python tests (GitHub uses 3.9 and 3.14)
python .githooks/test_privacy_hook.py -q
pip install -e "python/[numpy,validate]" pytest
# Include the optional adapter on the Python 3.14 candidate job:
pip install --only-binary=pyscf -e "python/[pyscf]" "pyscf==2.14.0"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pytest python/tests -q

# Optional focused corpus check, already included in pytest above
python conformance/run_conformance.py

# C++ build and tests with warnings as errors. Pass the flags through CXXFLAGS,
# not -DCMAKE_CXX_FLAGS, so CMake keeps its MSVC defaults. On MSVC use
# /W4 /WX, and add --config Release / -C Release for Visual Studio generators.
CXXFLAGS="-Wall -Wextra -Werror" cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
# CTest includes qvf_conformance and qvf_conformance_no_schema.

# C++ under AddressSanitizer + UndefinedBehaviorSanitizer, then the reader over
# the corpus (build/ is ignored by version control)
SAN="-fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g"
CC=clang CXX=clang++ CFLAGS="$SAN" CXXFLAGS="$SAN" \
  cmake -S cpp -B build/san -DCMAKE_BUILD_TYPE=Debug
cmake --build build/san
ctest --test-dir build/san --output-on-failure
for f in conformance/corpus/*.qvf; do
  ./build/san/read_qvf "$f" || break
done

# Fuzzing the ZIP reader and the validator (clang; not run in CI yet). Where
# clang has no libFuzzer runtime, the target uses a built-in mutation driver.
CC=clang CXX=clang++ cmake -S cpp -B build/fuzz -DCMAKE_BUILD_TYPE=Debug -DQVF_FUZZ=ON
cmake --build build/fuzz --target zip_reader_fuzzer
mkdir -p build/fuzz/corpus   # libFuzzer adds new inputs to the first directory only
./build/fuzz/zip_reader_fuzzer -max_total_time=60 build/fuzz/corpus conformance/corpus conformance/negative

# Generated tables: check them, and regenerate after editing the schema,
# registry.json or spec Appendix B (then regenerate the single header too)
python scripts/gen_kinds.py --check
python scripts/gen_kinds.py

# Generated single header: check it, and regenerate after editing cpp/include or cpp/src
python scripts/amalgamate.py --check
python scripts/amalgamate.py

# Documentation site (mkdocs.yml) in a separate environment; output in site/
python -m pip install -r scripts/requirements-docs.txt
mkdocs build --strict
```

The site is built from `docs/`. `scripts/mkdocs_hooks.py` adds the
top-level documents (README, spec, governance, changelog, roadmap,
contributing, security, conformance) and points links to other repository
files at the public source tree, so documents keep their repository-relative
links. A new page under `docs/` needs a `nav` entry in `mkdocs.yml`; strict
mode fails otherwise. `docs/api/python.md` is the Python API reference:
mkdocstrings renders the docstrings of the listed `qvf` modules there, reading
the sources without importing them, so a docstring is Markdown. A literal
block (such as a usage line with `[...]` options) is indented after a `::`
line, or strict mode reports a missing cross-reference. The build also runs
Doxygen (`cpp/Doxyfile`) over `cpp/include/qvf/` and publishes the C++ and C
reference under `api/cpp-reference/`; install Doxygen first, or set
`QVF_DOCS_NO_DOXYGEN=1` to build the rest with a placeholder. A Doxygen
warning fails strict mode. The headers keep plain `//` and `/* */` comments:
`scripts/doxygen_filter.py` presents a comment block directly above a
declaration, and a trailing comment after a member, to Doxygen as
documentation, so keep a header's comment next to what it describes. A
reference to an overloaded or template member names it without `()`, since
Doxygen reads `add()` as the overload without parameters. The Doxyfile sets
`HAVE_DOT = NO`, so Graphviz is not needed. `validation.anchors` in
`mkdocs.yml` makes a link to a missing heading anchor fail strict mode too.
CI does not build or publish the site.

### Publishing the documentation site

Publication is a maintainer action. `.github/workflows/docs-pages.yml` builds
the site in strict mode, with Doxygen, and deploys it to GitHub Pages. It runs
only when a maintainer dispatches it (Actions, "Docs site", Run workflow) from
`main`. No push, tag or schedule starts it, and the release workflow does not
call it.

One-time repository settings before the first dispatch:

1. Settings, Pages, Build and deployment: set Source to **GitHub Actions**.
   This creates the `github-pages` environment the deploy job uses.
2. Settings, Environments, `github-pages`: keep the deployment branch rule
   at `main` (the default), so only `main` can publish.
3. For a custom domain, set it under Settings, Pages and keep `site_url` in
   `mkdocs.yml` equal to the final address. The workflow publishes the site
   under `docs/` of the Pages base address with a forwarding root page, which
   matches the current `site_url` (`https://vibyverse.com/qvf/docs/`) when the
   Pages base is `https://vibyverse.com/qvf/`. With the default GitHub domain
   the site is at `https://vibyverse.github.io/qvf/docs/`.

Before dispatching, run `mkdocs build --strict` locally with Doxygen installed
and check the links of the built site.

### Private operator configuration

Deployment configuration is maintained separately from the format source.
Private literal terms can be supplied to the contributor guard through
`VIBE_PRIVACY_TERMS_FILE` or clone-local `privacy.termsFile`. Use an absolute
path to a UTF-8 file outside every Git checkout and object database, one
literal per line. Matching is case-insensitive; an explicitly configured
missing, empty or in-repository file blocks the check, including symlinks
through a checkout. Never commit the private policy or resolved home paths.

Keep operator configuration and audit evidence in external private storage.
The toolset uses an absolute `VIBE_PRIVATE_ROOT`, defaulting to
`${XDG_STATE_HOME:-$HOME/.local/state}/vibe-private`; it must resolve outside
Git repositories. A terms file can live under `qvf/config/` in that root.
Set `VIBE_PRIVACY_TERMS_FILE` to its absolute path explicitly; QVF does not
create or discover private policy files automatically. Use owner-only
directories (0700) and files (0600), and keep credentials in your existing
secret store. Ordinary build outputs and caches are unaffected.


### Release and adoption checklist

This checklist is for the release/updater owner. It does not add release CI,
publication or a second review as a gate on ordinary source changes.

1. **Freeze the intended source.** Check ownership and incoming commits before
   choosing a candidate. Keep the toolkit version consistent in `VERSION`,
   language packages, CMake, citation metadata, README/roadmap and changelog.
   Record public-contract decisions in `GOVERNANCE.md`; confirm compatibility,
   canonical units and the unchanged or explicitly authorized format version.
   Keep later Unreleased work outside the candidate.
2. **Check the source inventory.** Run the affected generators' `--check`
   modes (`gen_kinds.py`, `amalgamate.py`, `gen_fortran.py`). Build source
   archives with `scripts/make_archive.sh --output-dir` naming an external
   directory. Check that tar and ZIP contain exactly committed source bytes,
   including positive, negative and legacy corpus fixtures, license and notices.
   Preserve binary fixtures intentionally designed to be unreadable; reconcile
   their exact hashes with the approved publication review before exporting.
3. **Qualify the candidate.** Run the authorized candidate CI at the frozen
   commit: language/platform checks, shared finding codes, corpus parity,
   generated files and sanitizers. A failing job holds that candidate. Fix the
   cause, identify the resulting commit, and repeat only affected checks until
   a complete qualifying candidate result exists. Do not label pending or
   unverified scientific results as passing.
4. **Check installed consumption.** Build wheel/sdist and smoke-test the
   documented imports/CLI with declared extras. Verify CMake/pkg-config
   installation from outside the source checkout and exercise static/shared
   C/C++ clients. If the [local Conan recipe](packaging/conan/README.md) is
   included, run its `test_package` under the native-build lock and retain the
   recipe revision/profile. Claim Fortran support only for compiler/platform
   combinations actually checked; installed modules are compiler-specific.
5. **Release and publish through the approved path.** Only the authorized owner
   creates the immutable tag after acceptance. Inspect the actual publication
   receipt and anonymous destination readback, including source identity,
   inventory and provenance. Tag presence alone is not publication acceptance.
   Scanner policy, schema hosting, distribution naming and media-type
   registration remain maintainer decisions. A local Conan create is not a
   package-index upload; registry publication requires its own authorization.
6. **Adopt and record.** Each family consumer pins one release tag and vendors
   byte-identical contract files. Run the affected reader/writer or handoff
   check for that consumer, preserve original archives and report unsupported
   kinds explicitly. Save source/runtime/test identities and release evidence
   in the private operational record. Remove only owned temporary outputs
   after the evidence is durable; keep reusable caches and active jobs intact.
