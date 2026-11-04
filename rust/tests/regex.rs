//! The scrub pattern engine (`qvf::regex`) against Python's `re.subn` on a
//! battery of patterns and texts: the same replacement and count, and an
//! error exactly where Python refuses the pattern. Skipped without Python.

use std::process::Command;

use serde_json::{json, Value};

const PATTERNS: &[&str] = &[
    "alice",
    "x*",
    "a|ab",
    "ab|a",
    "a+?",
    "a*?",
    "a{2}",
    "a{2,}",
    "a{,2}",
    "a{1,3}?",
    r"\d+",
    r"\D+",
    r"\w+",
    r"\W",
    r"\s+",
    r"\S+",
    r"\bfoo\b",
    r"\Bo",
    "^foo",
    "foo$",
    "$",
    "^",
    "(?m)^x",
    "(?m)$",
    "(?s)a.b",
    "a.b",
    "(?i)HOME",
    "(?i)[a-c]+",
    r"/home/[^/\s]+",
    r"/Users/[^/]+/",
    "[^a-z]",
    "[a-]",
    "[]a]",
    r"[\d.]+",
    r"[\w.-]+@[\w.-]+",
    "(a*)*",
    "(a|b)*c",
    "(?:ab)+",
    "(?P<user>jdoe)",
    "a(?=b)",
    "a(?!b)",
    "(?#note)x",
    "a{2",
    "{",
    r"\x41",
    r"é",
    r"\.",
    r"\\",
    "colou?r",
    r"(\w+)@example\.org",
    r"\A\w",
    r"\w\Z",
    "[A-Z][a-z]*",
    "é+",
    "(?i)é",
    "",
    "a|",
    "|b",
    r"^\s*$",
    "(a|ab)(c|bcd)",
    // refused by both
    "*a",
    "a**",
    "(",
    ")",
    "[a",
    "a{3,2}",
    r"\q",
];

const TEXTS: &[&str] = &[
    "",
    "alice and Alice",
    "abxd xx",
    "aaa ab abc",
    "foo food foo",
    "/home/USER/run /home/runner/x\n/Users/Shared/y/z",
    "mail jdoe@example.org, x.y-z@sub.example.org",
    "line1\nx line2\nx\n",
    "a\nb a.b aXb",
    "HOME home Home",
    "colour color colouur",
    "12 3.14 v1.2.3",
    "café éé É",
    "Abc DEF ghi",
    "  \n\t",
    "abcd abd ac",
];

#[test]
fn patterns_match_python_re() {
    let python = std::env::var("PYTHON").unwrap_or_else(|_| "python3".into());
    let script = r#"
import json, re, sys
cases = json.load(sys.stdin)
out = []
for p, t in cases:
    try:
        out.append(list(re.compile(p).subn(lambda _: "<R>", t)))
    except re.error:
        out.append(None)
print(json.dumps(out))
"#;
    let cases: Vec<Value> = PATTERNS.iter().flat_map(|p| TEXTS.iter().map(move |t| json!([p, t]))).collect();
    let child = Command::new(&python)
        .args(["-c", script])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn();
    let Ok(mut child) = child else {
        println!("SKIP: no Python");
        return;
    };
    use std::io::Write;
    child.stdin.take().unwrap().write_all(serde_json::to_string(&cases).unwrap().as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    let expected: Vec<Value> = serde_json::from_slice(&out.stdout).unwrap();
    let mut failures = Vec::new();
    for (case, want) in cases.iter().zip(&expected) {
        let (p, t) = (case[0].as_str().unwrap(), case[1].as_str().unwrap());
        let got = match qvf::regex::Regex::new(p) {
            Ok(r) => {
                let (s, n) = r.replace_all(t, "<R>");
                json!([s, n])
            }
            Err(_) => Value::Null,
        };
        if &got != want {
            failures.push(format!("{p:?} on {t:?}: {got}, Python {want}"));
        }
    }
    println!("{} cases, {} differences", cases.len(), failures.len());
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
