//! The table-driven JSON Schema checker.
//!
//! `scripts/gen_kinds.py` compiles `spec/qvf_manifest.schema.json` into the
//! instruction table of `generated.rs`: one node per schema object, one op per
//! assertion keyword, with shared pools of strings, numbers, JSON literals and
//! compiled patterns. This module interprets that table. It is a port of
//! `js/src/schemaCheck.ts` (itself the same interpreter as the C++
//! `cpp/src/qvf_schema_check.hpp`), so a finding's location is the same JSON
//! pointer in every validator. The messages are worded independently.
//!
//! Patterns run on a small Pike VM over code points (linear time), with
//! ECMA-262 semantics: `$` matches only at the end of the string.

use serde_json::Value;

use crate::generated::{
    op, rx_op, type_bit, DEFS, LISTS, NODES, NUMBERS, OPS, PAIRS, PAYLOAD_DEFS, ROOT, RX_OPS, RX_PROGRAMS, RX_RANGES, RX_SOURCES,
    STRINGS, TRIPLES, TYPE_NAMES, VALUES,
};
use crate::json::{json_equal, show};

/// One schema error.
#[derive(Debug, Clone)]
pub struct SchemaError {
    /// RFC 6901 JSON pointer of the failing value inside the checked document.
    pub pointer: String,
    /// The same path unescaped and `/`-joined, as the Python validator prints it.
    pub location: String,
    pub message: String,
}

fn at(table: &[i32], i: i32) -> i32 {
    table[i as usize]
}

fn s(i: i32) -> &'static str {
    STRINGS[i as usize]
}

/// The JSON literals of the table (const / enum values), parsed once.
fn literal(i: i32) -> &'static Value {
    static LITERALS: std::sync::OnceLock<Vec<Value>> = std::sync::OnceLock::new();
    &LITERALS.get_or_init(|| VALUES.iter().map(|t| serde_json::from_str(t).unwrap_or(Value::Null)).collect())[i as usize]
}

/// The `$defs` entry specifying the JSON payload of a (kind, role) member.
pub fn payload_definition(kind: &str, role: &str) -> Option<&'static str> {
    PAYLOAD_DEFS.iter().find(|(k, r, _)| *k == kind && *r == role).map(|(_, _, d)| *d)
}

fn is_integral(v: &Value) -> bool {
    v.as_f64().is_some_and(|x| x.fract() == 0.0)
}

fn type_matches(mask: i32, v: &Value) -> bool {
    match v {
        Value::Null => mask & type_bit::NULL != 0,
        Value::Bool(_) => mask & type_bit::BOOLEAN != 0,
        Value::String(_) => mask & type_bit::STRING != 0,
        Value::Number(_) => mask & type_bit::NUMBER != 0 || (mask & type_bit::INTEGER != 0 && is_integral(v)),
        Value::Array(_) => mask & type_bit::ARRAY != 0,
        Value::Object(_) => mask & type_bit::OBJECT != 0,
    }
}

fn type_names(mask: i32) -> String {
    TYPE_NAMES.iter().filter(|(_, bit)| mask & bit != 0).map(|(n, _)| format!("'{n}'")).collect::<Vec<_>>().join(", ")
}

fn class_match(code: i32, first: i32, count: i32, c: u32) -> bool {
    let hit = (0..count).any(|k| {
        let r = (2 * (first + k)) as usize;
        c >= RX_RANGES[r] && c <= RX_RANGES[r + 1]
    });
    if code == rx_op::CLASS {
        hit
    } else {
        !hit
    }
}

/// True when compiled pattern `rx` matches anywhere in `text` (JSON Schema
/// patterns are unanchored). A Pike VM: every live thread advances in lock
/// step, so the run time is linear in the text and the program.
pub fn rx_search(rx: i32, text: &str) -> bool {
    let first = at(RX_PROGRAMS, 2 * rx);
    let count = at(RX_PROGRAMS, 2 * rx + 1) as usize;
    let cps: Vec<u32> = text.chars().map(u32::from).collect();
    let mut mark = vec![-1i64; count];
    let mut step: i64 = 0;
    let mut matched = false;
    let mut stack: Vec<i32> = Vec::new();

    let mut add = |list: &mut Vec<i32>, start: i32, pos: usize, step: i64, matched: &mut bool| {
        stack.push(start);
        while let Some(pc) = stack.pop() {
            if mark[pc as usize] == step {
                continue;
            }
            mark[pc as usize] = step;
            let o = 3 * (first + pc);
            let code = at(RX_OPS, o);
            if code == rx_op::JUMP {
                stack.push(at(RX_OPS, o + 1));
            } else if code == rx_op::SPLIT {
                stack.push(at(RX_OPS, o + 2));
                stack.push(at(RX_OPS, o + 1));
            } else if code == rx_op::LINE_START {
                if pos == 0 {
                    stack.push(pc + 1);
                }
            } else if code == rx_op::LINE_END {
                if pos == cps.len() {
                    stack.push(pc + 1);
                }
            } else if code == rx_op::MATCH {
                *matched = true;
            } else {
                list.push(pc);
            }
        }
    };

    let mut cur: Vec<i32> = Vec::new();
    add(&mut cur, 0, 0, step, &mut matched);
    let mut pos = 0;
    while pos < cps.len() && !matched {
        step += 1;
        let mut next: Vec<i32> = Vec::new();
        let c = cps[pos];
        for &pc in &cur {
            let o = 3 * (first + pc);
            if class_match(at(RX_OPS, o), at(RX_OPS, o + 1), at(RX_OPS, o + 2), c) {
                add(&mut next, pc + 1, pos + 1, step, &mut matched);
            }
        }
        add(&mut next, 0, pos + 1, step, &mut matched); // a match may start anywhere
        cur = next;
        pos += 1;
    }
    matched
}

fn escape_segment(seg: &str) -> String {
    seg.replace('~', "~0").replace('/', "~1")
}

struct Evaluator<'o> {
    path: Vec<String>,
    /// Collects errors; with `None` it only decides validity and stops early.
    out: Option<&'o mut Vec<SchemaError>>,
}

fn valid(node: i32, v: &Value) -> bool {
    Evaluator { path: Vec::new(), out: None }.eval(node, v)
}

/// Whether `key` is named by the node's own Properties or PatternProperties ops.
fn declared(first: i32, count: i32, key: &str) -> bool {
    for i in first..first + count {
        let code = at(OPS, 4 * i);
        let a = at(OPS, 4 * i + 1);
        let b = at(OPS, 4 * i + 2);
        if code == op::PROPERTIES {
            if (0..b).any(|k| key == s(at(PAIRS, a + 2 * k))) {
                return true;
            }
        } else if code == op::PATTERN_PROPERTIES && (0..b).any(|k| rx_search(at(PAIRS, a + 2 * k), key)) {
            return true;
        }
    }
    false
}

impl Evaluator<'_> {
    fn collect(&self) -> bool {
        self.out.is_some()
    }

    fn eval(&mut self, node: i32, v: &Value) -> bool {
        let first = at(NODES, 2 * node);
        let count = at(NODES, 2 * node + 1);
        let mut ok = true;
        for i in 0..count {
            if !self.apply(first + i, first, count, v) {
                ok = false;
                if !self.collect() {
                    return false;
                }
            }
        }
        ok
    }

    fn fail(&mut self, message: String) -> bool {
        if let Some(out) = self.out.as_deref_mut() {
            out.push(SchemaError {
                pointer: self.path.iter().map(|p| format!("/{}", escape_segment(p))).collect(),
                location: self.path.join("/"),
                message,
            });
        }
        false
    }

    fn descend(&mut self, node: i32, v: &Value, segment: &str) -> bool {
        if !self.collect() {
            return self.eval(node, v);
        }
        self.path.push(segment.to_owned());
        let ok = self.eval(node, v);
        self.path.pop();
        ok
    }

    fn apply(&mut self, opi: i32, first: i32, count: i32, v: &Value) -> bool {
        let code = at(OPS, 4 * opi);
        let a = at(OPS, 4 * opi + 1);
        let b = at(OPS, 4 * opi + 2);
        let c = at(OPS, 4 * opi + 3);
        let collect = self.collect();
        let sv = || show(Some(v));
        match code {
            op::FALSE => self.fail(format!("False schema does not allow {}", sv())),
            op::TYPE => type_matches(a, v) || self.fail(format!("{} is not of type {}", sv(), type_names(a))),
            op::CONST => json_equal(v, literal(a)) || self.fail(format!("{} was expected", VALUES[a as usize])),
            op::ENUM => {
                if (0..b).any(|k| json_equal(v, literal(at(LISTS, a + k)))) {
                    return true;
                }
                if !collect {
                    return false;
                }
                let names: Vec<&str> = (0..b).map(|k| VALUES[at(LISTS, a + k) as usize]).collect();
                self.fail(format!("{} is not one of [{}]", sv(), names.join(", ")))
            }
            op::REF => self.eval(a, v),
            op::ALL_OF => {
                let mut ok = true;
                for k in 0..b {
                    if !self.eval(at(LISTS, a + k), v) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::ANY_OF => {
                (0..b).any(|k| valid(at(LISTS, a + k), v))
                    || self.fail(format!("{} is not valid under any of the given schemas", sv()))
            }
            op::ONE_OF => {
                let mut matches = 0;
                let mut k = 0;
                while k < b && matches < 2 {
                    if valid(at(LISTS, a + k), v) {
                        matches += 1;
                    }
                    k += 1;
                }
                if matches == 1 {
                    return true;
                }
                let msg = if matches == 0 {
                    format!("{} is not valid under any of the given schemas", sv())
                } else {
                    format!("{} is valid under more than one of the given schemas", sv())
                };
                self.fail(msg)
            }
            op::NOT => !valid(a, v) || self.fail(format!("{} should not be valid under the given schema", sv())),
            op::IF => {
                if valid(a, v) {
                    b < 0 || self.eval(b, v)
                } else {
                    c < 0 || self.eval(c, v)
                }
            }
            op::PROPERTIES => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for k in 0..b {
                    let name = s(at(PAIRS, a + 2 * k));
                    if let Some(child) = o.get(name) {
                        if !self.descend(at(PAIRS, a + 2 * k + 1), child, name) {
                            ok = false;
                            if !collect {
                                return false;
                            }
                        }
                    }
                }
                ok
            }
            op::PATTERN_PROPERTIES => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for k in 0..b {
                    let rx = at(PAIRS, a + 2 * k);
                    for (key, value) in o {
                        if rx_search(rx, key) && !self.descend(at(PAIRS, a + 2 * k + 1), value, key) {
                            ok = false;
                            if !collect {
                                return false;
                            }
                        }
                    }
                }
                ok
            }
            op::ADDITIONAL_PROPERTIES => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                let mut extras: Vec<String> = Vec::new();
                for (key, value) in o {
                    if declared(first, count, key) {
                        continue;
                    }
                    if a < 0 {
                        extras.push(format!("'{key}'"));
                        if !collect {
                            return false;
                        }
                    } else if !self.descend(a, value, key) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                if !extras.is_empty() {
                    return self.fail(format!("additional properties are not allowed ({} unexpected)", extras.join(", ")));
                }
                ok
            }
            op::REQUIRED => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for k in 0..b {
                    let name = s(at(LISTS, a + k));
                    if !o.contains_key(name) {
                        ok = self.fail(format!("'{name}' is a required property"));
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::DEPENDENT_REQUIRED => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for k in 0..b {
                    let t = a + 3 * k;
                    let prop = s(at(TRIPLES, t));
                    if !o.contains_key(prop) {
                        continue;
                    }
                    for d in 0..at(TRIPLES, t + 2) {
                        let dep = s(at(LISTS, at(TRIPLES, t + 1) + d));
                        if !o.contains_key(dep) {
                            ok = self.fail(format!("'{dep}' is a dependency of '{prop}'"));
                            if !collect {
                                return false;
                            }
                        }
                    }
                }
                ok
            }
            op::DEPENDENT_SCHEMAS => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for k in 0..b {
                    if o.contains_key(s(at(PAIRS, a + 2 * k))) && !self.eval(at(PAIRS, a + 2 * k + 1), v) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::PROPERTY_NAMES => {
                let Value::Object(o) = v else { return true };
                let mut ok = true;
                for key in o.keys() {
                    if !self.eval(a, &Value::String(key.clone())) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::MIN_PROPERTIES => match v {
                Value::Object(o) if (o.len() as i32) < a => self.fail(format!("{} does not have enough properties", sv())),
                _ => true,
            },
            op::MAX_PROPERTIES => match v {
                Value::Object(o) if (o.len() as i32) > a => self.fail(format!("{} has too many properties", sv())),
                _ => true,
            },
            op::PREFIX_ITEMS => {
                let Value::Array(arr) = v else { return true };
                let mut ok = true;
                for (k, item) in arr.iter().enumerate().take(b.max(0) as usize) {
                    if !self.descend(at(LISTS, a + k as i32), item, &k.to_string()) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::ITEMS => {
                let Value::Array(arr) = v else { return true };
                if arr.len() <= b as usize {
                    return true;
                }
                if a < 0 {
                    return self.fail(format!("expected at most {b} items, found {}", arr.len()));
                }
                let mut ok = true;
                for (i, item) in arr.iter().enumerate().skip(b as usize) {
                    if !self.descend(a, item, &i.to_string()) {
                        ok = false;
                        if !collect {
                            return false;
                        }
                    }
                }
                ok
            }
            op::CONTAINS => {
                let Value::Array(arr) = v else { return true };
                let mut matches = 0;
                for item in arr {
                    if !valid(a, item) {
                        continue;
                    }
                    matches += 1;
                    if c >= 0 && matches > c {
                        return self.fail(format!("too many items match the given schema (expected at most {c})"));
                    }
                }
                if matches >= b {
                    return true;
                }
                let msg = if matches == 0 {
                    format!("{} does not contain items matching the given schema", sv())
                } else {
                    format!("too few items match the given schema (expected at least {b}, {matches} matched)")
                };
                self.fail(msg)
            }
            op::MIN_ITEMS => match v {
                Value::Array(arr) if (arr.len() as i64) < i64::from(a) => {
                    self.fail(format!("{} is too short (at least {a} items)", sv()))
                }
                _ => true,
            },
            op::MAX_ITEMS => match v {
                Value::Array(arr) if (arr.len() as i64) > i64::from(a) => {
                    self.fail(format!("{} is too long (at most {a} items)", sv()))
                }
                _ => true,
            },
            op::UNIQUE_ITEMS => {
                let Value::Array(arr) = v else { return true };
                for i in 0..arr.len() {
                    for j in i + 1..arr.len() {
                        if json_equal(&arr[i], &arr[j]) {
                            return self.fail(format!("{} has non-unique elements", sv()));
                        }
                    }
                }
                true
            }
            op::MINIMUM | op::MAXIMUM | op::EXCLUSIVE_MINIMUM | op::EXCLUSIVE_MAXIMUM => {
                let Some(x) = v.as_f64() else { return true };
                let lim = NUMBERS[a as usize];
                match code {
                    op::MINIMUM => x >= lim || self.fail(format!("{} is less than the minimum of {lim}", sv())),
                    op::MAXIMUM => x <= lim || self.fail(format!("{} is greater than the maximum of {lim}", sv())),
                    op::EXCLUSIVE_MINIMUM => {
                        x > lim || self.fail(format!("{} is less than or equal to the minimum of {lim}", sv()))
                    }
                    _ => x < lim || self.fail(format!("{} is greater than or equal to the maximum of {lim}", sv())),
                }
            }
            op::MULTIPLE_OF => {
                let Some(x) = v.as_f64() else { return true };
                let d = NUMBERS[a as usize];
                let bad = if x.fract() == 0.0 && d.fract() == 0.0 { x % d != 0.0 } else { (x / d).fract() != 0.0 };
                !bad || self.fail(format!("{} is not a multiple of {d}", sv()))
            }
            op::MIN_LENGTH => match v {
                Value::String(t) if (t.chars().count() as i64) < i64::from(a) => {
                    let msg = if a == 1 {
                        format!("{} should be non-empty", sv())
                    } else {
                        format!("{} is too short (at least {a} characters)", sv())
                    };
                    self.fail(msg)
                }
                _ => true,
            },
            op::MAX_LENGTH => match v {
                Value::String(t) if (t.chars().count() as i64) > i64::from(a) => {
                    self.fail(format!("{} is too long (at most {a} characters)", sv()))
                }
                _ => true,
            },
            op::PATTERN => match v {
                Value::String(t) if !rx_search(a, t) => {
                    self.fail(format!("{} does not match '{}'", sv(), RX_SOURCES[a as usize]))
                }
                _ => true,
            },
            _ => self.fail(format!("unknown schema opcode {code}")),
        }
    }
}

fn check_node(node: i32, doc: &Value) -> Vec<SchemaError> {
    let mut out = Vec::new();
    Evaluator { path: Vec::new(), out: Some(&mut out) }.eval(node, doc);
    out
}

/// Every schema error of a `manifest.json` document against the root schema.
pub fn check_manifest(manifest: &Value) -> Vec<SchemaError> {
    check_node(ROOT, manifest)
}

/// Every schema error of a payload against the `$defs` entry `def` (nothing
/// for an unknown entry).
pub fn check_definition(def: &str, doc: &Value) -> Vec<SchemaError> {
    match DEFS.iter().find(|(n, _)| *n == def) {
        Some((_, node)) => check_node(*node, doc),
        None => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_pattern_compiles_and_runs() {
        for (i, src) in RX_SOURCES.iter().enumerate() {
            // A smoke test: the VM terminates on a few inputs for every program.
            for text in ["", "x_vendor", "abc123", "\u{e9}"] {
                let _ = rx_search(i as i32, text);
            }
            assert!(!src.is_empty());
        }
    }

    #[test]
    fn manifest_minimal_errors() {
        let errs = check_manifest(&serde_json::json!({}));
        assert!(errs.iter().any(|e| e.message.contains("required property")));
    }
}
