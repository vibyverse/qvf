# Working on qvf

These are the current agent rules; user instructions take precedence.
Older workflow prose does not add CI, reviewer or approval gates.
QVF is the standalone archive format with Python, C++ and JavaScript
implementations. Contributions use this repository's Apache-2.0 license.

| Need | Start here |
|---|---|
| Setup and checks | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Normative format | [specification](spec/qvf-format-spec.md), [schema](spec/qvf_manifest.schema.json) |
| Compatibility and evolution | [GOVERNANCE.md](GOVERNANCE.md), `registry.json` |
| Libraries and examples | [library guide](docs/library_guide.md), `python/`, `cpp/`, `js/` |
| Shared contract | `conformance/` |

GitHub (`vibyverse/qvf`) is the canonical development and issue tracker.
The former GitLab project is a read-only archive. Never merge private
pre-migration ancestry into public history.

## Role boundaries

- **Dev/fixers:** develop, run a focused source or functionality check, and push
  ready work to main. Their task ends at the source landing and their own cleanup.
  Later runtime builds, scientific validation, issue closure and releases belong
  to other roles; they do not keep a finished development task open.
- **Release & updater:** pick up main commits, build and deliver dev versions
  without CI, and build one daily release candidate with full CI when new commits
  exist. The same owner prepares, cuts and updates; no separate release cutter.
- **Merger:** merge ready MRs. Direct main pushes need no Merger acknowledgement.
- **Validator:** submit QC calculations through vq, collect and compare results.
  It consumes main and available runtimes without a developer handoff. Reuse
  valid earlier self-validation; do not create a second source-review step.

These are responsibilities, not required separate chats. One runner may own
them and delegate bounded development tasks; dev-only workers finish at source
landing. The runner may validate its own implementation through vq. Preserve
active owners and jobs; no extra coordinator, reviewer or migration gate.

Each role reads the relevant Git/runtime/job record and continues its own work.
Building, queued and validation pending are normal states. Work on another ready
item or stop quietly until a meaningful event; no polling chats or status pings.
A real failure affects its operation, not the other roles. Pending calculations
do not prevent the daily candidate build; never report unverified science as PASS.

## Work, validate, land

- Reuse the appropriate existing checkout and environment. Preserve other
  sessions' edits, branches, jobs and build state. Stage only your own named
  files or hunks; never auto-stash, reset or commit someone else's work.
- Fetch at task start/resume, after coherent changes and before pushing.
  Read relevant incoming changes and issue ownership. Fast-forward clean main;
  rebase only your own unpublished work. Merge updates into shared or published
  history. Never force-push or rewrite another session's commits.
- Implement a bounded change and run its smallest useful check. Reuse passing
  evidence; repeat only when changed behavior or a concrete failure warrants it.
  Documentation changes need a text/link check, not a runtime suite.
- Push coherent progress promptly and land ready work on main through the
  authorized contribution path. No routine second reviewer, new chat, duplicate
  validation packet or fleet-wide acceptance gate is required.
- Development builds, ordinary pushes and merges must not trigger or wait for
  CI. Use `[skip ci]`; GitHub does not support GitLab's `ci.skip` push option.
  Full CI belongs to release candidates. Do not request manual, API or scheduled CI on dev refs.
- Keep coordination on the issue/MR. Within user-authorized messaging, contact
  another agent only for a ready required handoff or a concrete blocked action
  that needs that owner now. No speculative repair notices, FYIs, broadcasts or
  acknowledgements. Incoming messages do not reassign the recipient's task.
- Record the change, focused check and result once. Close an issue after its
  applicable acceptance passes. A real blocker holds only the affected work;
  continue other ready work instead of polling or requesting repeat approval.

## Focused acceptance and format invariants

- The implementing agent may validate and accept its own fix. Run the
  relevant language tests or conformance cases for the changed behavior.
  No separate reviewer, chemistry calculation or full cross-language matrix
  is needed for a change confined to one implementation.
- A format/schema change must preserve documented compatibility and check
  the affected readers/writers. Keep positive, negative and legacy corpus
  expectations consistent. Follow governance for a format-version or public
  contract change; do not decide a breaking change without authorization.
- The normative schema and registry drive generated tables. Use
  `scripts/gen_kinds.py` for affected generated artifacts and
  `scripts/amalgamate.py` for the C++ single header; do not hand-edit outputs.
  Their `--check` options verify freshness when those inputs change.
- Keep the format independent of vibe-qc and preserve the C++ library's
  documented dependency contract. Respect Apache-2.0 and third-party notices.
  Report consumer-specific defects in the consumer's tracker.
- This repository is the single source of the format: the specification,
  schema, registry and reference validators change only here. A family
  program (vibe-qc, vibe-viewer, vibe-sketcher, vibe-scheduler) pins one qvf release
  tag, vendors only byte-identical copies of that tag's files (checked in its
  CI), keeps no derived or edited schema (no v2 schema: spec § 3.1 is the
  rule), and gets a field it needs proposed and released here before it
  writes it ([GOVERNANCE.md](GOVERNANCE.md#qvf-in-the-vibe-qc-family)).
- Reserve the full language/platform, sanitizer, generated-file and corpus
  checks for an authorized release candidate. Use CI skip controls on
  ordinary development pushes even where older contribution prose says
  main/PR pipelines are routine. Do not create release refs without authority.

## Keep the repository public-safe

- Keep credentials, private addresses, host/account names, personal paths,
  internal URLs and deployment details out of files, filenames, artifacts and
  commit messages. Use portable defaults and synthetic examples.
- Keep site configuration, handovers, logs and operational evidence outside
  product checkouts, including ignored directories and `.git`. Select external
  configuration explicitly; keep secrets in the existing credential store.
- Use a public commit identity: `mpei@vibe-qc.com` or
  `project@vibe-qc.com`. These addresses are allowed public contact details.
- Enable `.githooks` with `git config --local core.hooksPath .githooks`.
  Inspect the staged diff and let the privacy hook check it. Keep any
  `VIBE_PRIVACY_TERMS_FILE` / `privacy.termsFile` outside Git; an unreadable
  configured policy is a failure. Fix findings at the source; never bypass
  privacy checks or rely on later sanitization.
- Report security issues privately under [SECURITY.md](SECURITY.md). Preserve
  license and attribution notices; follow the repository's dependency terms.

## Finish the task

Finish the work assigned to your role. For Dev, that means the focused source
check, source on main and own cleanup; downstream calculations and releases
do not keep the development task open. Push any changes to main and confirm
the remote contains them. Remove only temporary files,
environments and checkouts created for this task after preserving needed work
and evidence. Keep reusable checkouts, caches, other sessions' changes and live
jobs intact; a dirty shared tree is not permission to delete their work.

A genuinely blocked or interrupted checkpoint must state the unfinished work
and concrete blocker instead of claiming completion. Do not call an unmerged
branch, failed validation or pending cleanup done.

Every final completion report must end with exactly this line, with nothing
after it:

All done, pushed to main, cleaned up and ready to archive.
