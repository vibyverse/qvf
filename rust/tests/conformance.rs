//! The conformance runner (`qvf::conformance`, the Rust port of
//! `python/qvf/conformance.py`) over the shared corpus:
//!
//! * `conformance/corpus`: every archive validates without errors, has the
//!   expected section kinds and source program, and reproduces every check of
//!   its `expected.json` through the reader.
//! * `conformance/negative`: every archive is judged as its `expected.json`
//!   says: rejected with at least the expected finding codes, or, for a
//!   warning-only archive, valid with the warning. `limits` overrides apply.
//! * `conformance/legacy`: the same rules as `negative/`, run separately.
//!
//! Every archive of every corpus is run; none is skipped. Each test prints
//! its counts (`cargo test --test conformance -- --nocapture`).

use std::path::{Path, PathBuf};

fn conformance_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("conformance")
}

fn run(label: &str, positive: bool) {
    let dir = conformance_dir().join(label);
    let outcomes = qvf::conformance::run_dir(&dir, positive).unwrap();
    assert!(!outcomes.is_empty(), "no *.expected.json in {}", dir.display());
    let failed: Vec<String> = outcomes
        .iter()
        .filter(|o| !o.failures.is_empty())
        .map(|o| format!("{}\n     - {}", o.archive, o.failures.join("\n     - ")))
        .collect();
    println!(
        "{label} conformance (Rust): {} passed, {} failed, 0 skipped of {}",
        outcomes.len() - failed.len(),
        failed.len(),
        outcomes.len()
    );
    for f in &failed {
        println!("FAIL {f}");
    }
    assert!(failed.is_empty(), "{label} conformance: {} archive(s) failed", failed.len());
}

#[test]
fn positive_corpus() {
    run("corpus", true);
}

#[test]
fn negative_corpus() {
    run("negative", false);
}

#[test]
fn legacy_corpus() {
    run("legacy", false);
}
