//! JSON helpers shared by the reader and the validator: the nesting-depth
//! scan that guards the parser (spec § 6.1), strict UTF-8, parsing, and the
//! value predicates the validator ports from the JavaScript reference.

use serde_json::{Map, Value};

/// Nesting depth of the JSON arrays and objects in `raw`, capped at
/// `limit + 1`. A top-level `[]` or `{}` has depth 1 and a scalar depth 0.
/// The scan is linear, skips string contents (with escapes), and stops as
/// soon as the limit is exceeded. Malformed JSON is left to the parser.
pub fn json_depth(raw: &[u8], limit: f64) -> f64 {
    let mut depth: i64 = 0;
    let mut deepest: i64 = 0;
    let mut in_string = false;
    let mut escaped = false;
    for &c in raw {
        if in_string {
            if escaped {
                escaped = false;
            } else if c == b'\\' {
                escaped = true;
            } else if c == b'"' {
                in_string = false;
            }
            continue;
        }
        match c {
            b'"' => in_string = true,
            b'[' | b'{' => {
                depth += 1;
                if depth > deepest {
                    deepest = depth;
                    if deepest as f64 > limit {
                        return deepest as f64;
                    }
                }
            }
            b']' | b'}' => depth -= 1,
            _ => {}
        }
    }
    deepest as f64
}

/// `raw` without a leading UTF-8 byte order mark.
fn strip_bom(raw: &[u8]) -> &[u8] {
    raw.strip_prefix(b"\xEF\xBB\xBF".as_slice()).unwrap_or(raw)
}

/// Decode `raw` as strict UTF-8, skipping a leading byte order mark.
pub fn decode_utf8(raw: &[u8]) -> Result<&str, std::str::Utf8Error> {
    std::str::from_utf8(strip_bom(raw))
}

/// Parse UTF-8 JSON bytes.
pub fn parse_json_bytes(raw: &[u8]) -> Result<Value, String> {
    let text = decode_utf8(raw).map_err(|e| format!("invalid UTF-8 ({e})"))?;
    serde_json::from_str(text).map_err(|e| e.to_string())
}

pub(crate) type Obj = Map<String, Value>;

pub(crate) fn as_obj(v: Option<&Value>) -> Option<&Obj> {
    v.and_then(Value::as_object)
}

/// A JSON number with an integral value (JavaScript `Number.isInteger`).
pub(crate) fn int_of(v: Option<&Value>) -> Option<f64> {
    v.and_then(Value::as_f64).filter(|x| x.fract() == 0.0)
}

pub(crate) fn num_of(v: Option<&Value>) -> Option<f64> {
    v.and_then(Value::as_f64)
}

pub(crate) fn str_of(v: Option<&Value>) -> Option<&str> {
    v.and_then(Value::as_str)
}

/// Python-style truthiness, so loosely typed flags behave as in the
/// reference validators.
pub(crate) fn truthy(v: Option<&Value>) -> bool {
    match v {
        None | Some(Value::Null) => false,
        Some(Value::Bool(b)) => *b,
        Some(Value::Number(n)) => n.as_f64() != Some(0.0),
        Some(Value::String(s)) => !s.is_empty(),
        Some(Value::Array(a)) => !a.is_empty(),
        Some(Value::Object(o)) => !o.is_empty(),
    }
}

/// JSON equality as jsonschema decides it: numbers compare by value
/// (1 equals 1.0), booleans are not numbers, containers compare deeply.
pub fn json_equal(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => x.as_f64() == y.as_f64(),
        (Value::Array(x), Value::Array(y)) => x.len() == y.len() && x.iter().zip(y).all(|(p, q)| json_equal(p, q)),
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| json_equal(v, w)))
        }
        _ => a == b,
    }
}

/// A compact rendering of a JSON value for messages.
pub(crate) fn show(v: Option<&Value>) -> String {
    let Some(v) = v else { return "nothing".into() };
    let text = v.to_string();
    if text.chars().count() > 80 {
        let head: String = text.chars().take(77).collect();
        format!("{head}...")
    } else {
        text
    }
}

/// A number as JavaScript prints it: integers without a fraction.
pub(crate) fn num(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e21 {
        format!("{}", x as i64)
    } else {
        format!("{x}")
    }
}

/// `%.10g`-like rendering used by the normalization messages.
pub(crate) fn g10(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return num(x);
    }
    let digits = 10 - 1 - x.abs().log10().floor() as i32;
    if (0..=20).contains(&digits) {
        let s = format!("{x:.*}", digits as usize);
        let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
        s
    } else {
        format!("{x:.9e}")
    }
}

pub(crate) fn escape_pointer(part: &str) -> String {
    part.replace('~', "~0").replace('/', "~1")
}

pub(crate) fn quote(s: &str) -> String {
    format!("'{s}'")
}

/// A list of integers as JavaScript joins an array: `1,2,3`.
pub(crate) fn join(v: &[f64]) -> String {
    v.iter().map(|x| num(*x)).collect::<Vec<_>>().join(",")
}

/// A JSON document that keeps the order of object keys (and nothing else):
/// `serde_json::Value` sorts them, the Python CLI walks them in file order.
enum Ordered {
    Other,
    Object(Vec<(String, Ordered)>),
    Array(Vec<Ordered>),
}

impl<'de> serde::Deserialize<'de> for Ordered {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> serde::de::Visitor<'de> for V {
            type Value = Ordered;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("JSON")
            }
            fn visit_bool<E>(self, _: bool) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_i64<E>(self, _: i64) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_u64<E>(self, _: u64) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_f64<E>(self, _: f64) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_str<E>(self, _: &str) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_unit<E>(self) -> Result<Ordered, E> {
                Ok(Ordered::Other)
            }
            fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Ordered, A::Error> {
                let mut out = Vec::new();
                while let Some(v) = seq.next_element()? {
                    out.push(v);
                }
                Ok(Ordered::Array(out))
            }
            fn visit_map<A: serde::de::MapAccess<'de>>(self, mut map: A) -> Result<Ordered, A::Error> {
                let mut out: Vec<(String, Ordered)> = Vec::new();
                while let Some((k, v)) = map.next_entry::<String, Ordered>()? {
                    // A repeated key keeps its first position and its last value, as a Python dict.
                    match out.iter_mut().find(|(key, _)| *key == k) {
                        Some(slot) => slot.1 = v,
                        None => out.push((k, v)),
                    }
                }
                Ok(Ordered::Object(out))
            }
        }
        d.deserialize_any(V)
    }
}

/// The member roles of every entry of `sections` in a raw manifest, in file
/// order (empty for an entry that is not an object with an object
/// `members`), or `None` when the bytes are not such a manifest.
pub fn member_roles_in_file_order(raw: &[u8]) -> Option<Vec<Vec<String>>> {
    let Ordered::Object(root) = serde_json::from_slice::<Ordered>(raw).ok()? else { return None };
    let (_, Ordered::Array(sections)) = root.iter().find(|(k, _)| k == "sections")? else { return None };
    Some(
        sections
            .iter()
            .map(|s| match s {
                Ordered::Object(fields) => match fields.iter().find(|(k, _)| k == "members") {
                    Some((_, Ordered::Object(m))) => m.iter().map(|(k, _)| k.clone()).collect(),
                    _ => Vec::new(),
                },
                _ => Vec::new(),
            })
            .collect(),
    )
}
