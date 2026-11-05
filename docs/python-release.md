# Python package release

The first PyPI release candidate is `qvf-writer==0.7.0`, the same version as
toolkit 0.7.0. It includes the optional PySCF adapter. It does not modify
format version 1 or frozen consumers' v0.6.0 pins. Python package metadata,
`qvf.__version__`, `VERSION` and the other language packages all identify
0.7.0. The candidate has not yet been uploaded. Do not advertise an index install
as available until anonymous installation succeeds.

## Prepare and verify

From the repository root:

```sh
python -m pip install build twine
python -m build python --outdir /tmp/qvf-dist
python -m twine check --strict /tmp/qvf-dist/*
```

Use a fresh output directory. The default build makes an sdist, then builds
the wheel from that sdist, proving the source distribution is self-contained.
Inspect both archives: include the package, schema, LICENSE and NOTICE;
exclude caches, local configuration, credentials and operational records.
The sdist also includes the runnable tutorial. It excludes `python/tests`,
which need the repository checkout (conformance corpus and documentation),
so run the test suite from the checkout. The wheel requires no native
compilation and no base dependencies.

Create a fresh environment, install the wheel with `--no-deps`, and run
`scripts/smoke_python_package.py` from outside the checkout with `python -I`.
Then install the wheel's `pyscf,validate` extras into a second fresh
environment and run the [tutorial](tutorials/pyscf-to-qvf.md), `qvf validate`,
`qvf-validate` and `pip check`. Record package hashes, source commit and tested
Python/PySCF versions. Existing adapter tests provide broader method coverage;
the tutorial is the installed-distribution acceptance check.

## One-time account setup

The maintainer creates/verifies a PyPI account with two-factor authentication,
then adds a **pending GitHub Trusted Publisher** at
<https://pypi.org/manage/account/publishing/>:

| Field | Value |
| --- | --- |
| PyPI project name | `qvf-writer` |
| GitHub owner | `vibyverse` |
| Repository | `qvf` |
| Workflow filename | `publish-python.yml` |
| Environment | `pypi` |

Configure the GitHub repository environment `pypi` with an appropriate release
reviewer. The workflow requests its short-lived OIDC identity only in the
publish job, which downloads the checked artifacts and runs no project code.
No API token or password needs to be shared with a contributor.
Pending publishers do not reserve names. See
[PyPI's pending publisher instructions](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/).

### Trusted publisher checklist

Only the maintainer can do these steps. They need PyPI and GitHub
repository administration rights. Nothing in the repository changes.

1. **PyPI account.** Enable two-factor authentication and keep recovery
   codes in the existing credential store. If the project will belong to a
   PyPI organization, create or join it first, so that the project is created
   under the right owner.
2. **Check the name.** Confirm that `qvf-writer` is still unregistered, or
   that the maintainer already owns it. A pending publisher does not reserve
   the name. If someone registers it first, the first upload fails and the
   pending publisher becomes useless.
3. **Add the pending publisher** with exactly the values in the table above.
   The workflow filename is the file name only, without `.github/workflows/`.
   The environment name must match the `environment.name` of the publish
   job. PyPI compares these claims with the short-lived OIDC token that
   GitHub issues to that job.
4. **Create the GitHub environment `pypi`** under the repository settings
   (Environments). Add at least one required reviewer who is a maintainer.
   Enable *Prevent self-review* when more than one maintainer exists. Limit
   deployment branches and tags to the default branch `main`, because the
   workflow is dispatched from `main` and checks out the release tag itself.
   Store no PyPI secret in this environment or in repository secrets.
5. **Protect the release refs.** Add a tag ruleset that allows only
   maintainers to create, update or delete `python-v*` tags, and keep branch
   protection on `main`. The workflow file on `main` and the tag both
   determine what is uploaded, so both must be maintainer controlled.
6. **Rehearse without upload.** Run the workflow with `publish=false`. This
   requests no OIDC token, needs no environment approval and only builds,
   checks and stores the artifacts. Download them and compare their SHA256
   hashes with a local `python -m build` of the same tag. Wheel hashes can
   differ between builds because of archive timestamps. Then compare the
   archive contents instead.
7. **First upload.** After authorization, run with `publish=true` and approve
   the `pypi` deployment. On success PyPI creates the project and converts
   the pending publisher into an ordinary trusted publisher for it. Check on
   the project's *Publishing* page that exactly one publisher is listed with
   the expected repository, workflow and environment.
8. **Verify the upload.** From a clean environment, run
   `python -m pip download --no-deps "qvf-writer==<version>"` for the wheel
   and `--no-binary :all:` for the sdist, and compare their SHA256 hashes with
   the workflow artifacts. The PyPI file pages show the provenance
   attestation that the publish action generates. It names the repository,
   workflow and run that produced each file.

To revoke publishing, remove the trusted publisher on PyPI and disable or
delete the `pypi` environment. A token leak is not possible because no
long-lived token exists. A compromised workflow run can still upload within
its token lifetime of about 15 minutes. That is why the environment reviewer
and the tag ruleset matter.

TestPyPI rehearsals need a separate trusted publisher on
<https://test.pypi.org/> and a workflow that targets its upload URL. The
current workflow does not do that. Adding it is a separate, authorized
workflow change.

## Publish after authorization

1. Freeze and qualify the exact Python candidate, including the relevant
   Python candidate checks. Preserve the existing toolkit tags and pins.
2. Land the clean candidate source on the canonical GitHub repository.
   Verify that the candidate includes this workflow and the exact qualified
   source. Do not push private pre-migration ancestry to the public repository.
3. Create the immutable public package tag `python-v0.7.0` at that verified
   snapshot. This is a Python distribution tag, not a replacement toolkit tag.
4. Manually run **Python package release** from the trusted default branch,
   with `release_tag=python-v0.7.0` and `publish=false` to build/check it.
   The workflow verifies the tag matches the package version.
5. After the maintainer authorizes upload, run with the same tag and
   `publish=true`, and approve the `pypi` environment. The workflow rebuilds
   and checks the exact tagged source before uploading those artifacts.
6. In another clean environment install
   `qvf-writer[pyscf,validate]==0.7.0` from PyPI, rerun the tutorial, and
   verify the index files and hashes. Then update installation documentation
   from “prepared candidate” to “published candidate”.

A later Python release needs its own version and artifacts. Never overwrite
an uploaded version. Base Python 3.9 support does not imply PySCF support on
every Python/platform combination. The workflow tests the optional adapter
on Linux/Python 3.13; report other verified platforms separately.
