//! Resource caps applied to untrusted archives (spec § 6.1).
//!
//! Field names and default values are those of the Python `Limits` dataclass,
//! the C++ `Limits` struct and the JavaScript `Limits`, so the `limits`
//! objects of the conformance corpus apply unchanged. Sizes are bytes. The
//! values are `f64` so that `f64::INFINITY` disables a cap.

use serde_json::Value;

const MIB: f64 = 1024.0 * 1024.0;
const GIB: f64 = 1024.0 * MIB;

/// The resource caps of spec § 6.1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Limits {
    /// Per ZIP entry, uncompressed bytes.
    pub max_member_bytes: f64,
    /// Sum over all ZIP entries, uncompressed bytes.
    pub max_total_bytes: f64,
    /// Optional application ratio policy; 0 disables it (the default).
    pub max_ratio: f64,
    /// Entries smaller than this (uncompressed) are exempt from the ratio check.
    pub ratio_exempt_bytes: f64,
    /// Largest allowed `manifest.json` or JSON member, uncompressed bytes.
    pub max_json_bytes: f64,
    /// Deepest allowed nesting of JSON arrays and objects.
    pub max_json_depth: f64,
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            max_member_bytes: 8.0 * GIB,
            max_total_bytes: 64.0 * GIB,
            max_ratio: 0.0,
            ratio_exempt_bytes: MIB,
            max_json_bytes: 256.0 * MIB,
            max_json_depth: 64.0,
        }
    }
}

impl Limits {
    /// The default limits with the overrides of a JSON object applied, such
    /// as the `limits` of a negative corpus `expected.json`. `null` means no
    /// overrides. An unknown key or a value that is not a non-negative number
    /// is an error.
    pub fn from_json(overrides: Option<&Value>) -> Result<Limits, String> {
        let mut out = Limits::default();
        let obj = match overrides {
            None | Some(Value::Null) => return Ok(out),
            Some(Value::Object(o)) => o,
            Some(_) => return Err("limits must be a JSON object".into()),
        };
        for (key, value) in obj {
            let v = value
                .as_f64()
                .filter(|x| *x >= 0.0)
                .ok_or_else(|| format!("limit {key} must be a non-negative number, got {value}"))?;
            match key.as_str() {
                "max_member_bytes" => out.max_member_bytes = v,
                "max_total_bytes" => out.max_total_bytes = v,
                "max_ratio" => out.max_ratio = v,
                "ratio_exempt_bytes" => out.ratio_exempt_bytes = v,
                "max_json_bytes" => out.max_json_bytes = v,
                "max_json_depth" => out.max_json_depth = v,
                _ => return Err(format!("unknown limit {key:?}")),
            }
        }
        Ok(out)
    }
}
