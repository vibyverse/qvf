//! Validate QVF archives and print one JSON report per archive.
//!
//! ```sh
//! cargo run --example validate -- archive.qvf [more.qvf ...] [--limits '{"max_ratio": 100}']
//! ```
//!
//! Exits 1 when any archive has an error finding.

use qvf::{validate_file, Limits, ValidateOptions};
use serde_json::json;

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut options = ValidateOptions::default();
    if let Some(k) = args.iter().position(|a| a == "--limits") {
        let text = args.get(k + 1).cloned().unwrap_or_default();
        let value: serde_json::Value = serde_json::from_str(&text).unwrap_or_else(|e| {
            eprintln!("--limits is not JSON: {e}");
            std::process::exit(2)
        });
        options.limits = Limits::from_json(Some(&value)).unwrap_or_else(|e| {
            eprintln!("{e}");
            std::process::exit(2)
        });
        args.drain(k..(k + 2).min(args.len()));
    }
    if args.is_empty() {
        eprintln!("usage: validate ARCHIVE.qvf [...] [--limits JSON]");
        std::process::exit(2);
    }
    let mut ok = true;
    for path in &args {
        let report = validate_file(path, &options);
        ok &= report.ok;
        let findings: Vec<_> = report
            .findings
            .iter()
            .map(|f| json!({"code": f.code.as_str(), "severity": f.severity.as_str(), "where": f.location, "message": f.message}))
            .collect();
        println!("{}", json!({"archive": path, "ok": report.ok, "n_sections": report.n_sections, "findings": findings}));
    }
    std::process::exit(if ok { 0 } else { 1 });
}
