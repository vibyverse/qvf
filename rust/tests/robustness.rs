//! Damaged archives are findings, never panics: every positive corpus archive
//! truncated at several lengths and with single bytes flipped goes through the
//! validator and the reader.

use std::path::Path;

use qvf::{validate_bytes, Limits, Reader, ValidateOptions};

#[test]
fn damaged_archives_do_not_panic() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/corpus");
    let mut paths: Vec<_> = std::fs::read_dir(&dir)
        .expect("corpus directory")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "qvf"))
        .collect();
    paths.sort();
    assert!(!paths.is_empty());
    let options = ValidateOptions::default();
    let mut runs = 0;
    for path in paths {
        let bytes = std::fs::read(&path).expect("archive");
        let n = bytes.len();
        let mut variants: Vec<Vec<u8>> = Vec::new();
        for k in 1..=12 {
            variants.push(bytes[..n * k / 13].to_vec());
        }
        // A deterministic spread of single-byte corruptions, the tail (central
        // directory) sampled more densely than the body.
        for k in 0..24usize {
            let at = if k < 12 { n * k / 12 } else { n - 1 - ((k - 12) * 97) % n.min(400) };
            let mut v = bytes.clone();
            v[at] ^= 0x5a;
            variants.push(v);
        }
        for v in variants {
            let report = validate_bytes(v.clone(), &options);
            assert_eq!(report.ok, !report.findings.iter().any(|f| f.code.as_str().starts_with('E')));
            if let Ok(reader) = Reader::from_bytes(v, Limits::default()) {
                for s in reader.sections() {
                    for role in s.members().keys() {
                        let _ = reader.read_bytes(s.id, role);
                        let _ = reader.read_tensor(s.id, role);
                        let _ = reader.read_json(s.id, role);
                    }
                }
            }
            runs += 1;
        }
    }
    println!("{runs} damaged archives validated and read without a panic");
}
