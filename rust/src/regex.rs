//! A small regular-expression engine for the redaction patterns of
//! `qvf scrub`: the commonly used subset of Python's `re` syntax, matched
//! with Python's leftmost-first semantics by a Pike VM in time linear in the
//! text (per lookahead), so a pattern cannot make the scrub hang.
//!
//! Supported: literals and escapes (`\n \t \r \f \v \a \0 \xhh \uhhhh
//! \Uhhhhhhhh`, escaped punctuation), `.`, classes `[...]` with ranges,
//! negation and `\d \w \s` (and their negations, also outside classes),
//! anchors `^ $ \A \Z \b \B`, groups `(...)`, `(?:...)`, `(?P<name>...)`,
//! lookaheads `(?=...)` and `(?!...)`, comments `(?#...)`, alternation, the
//! greedy and lazy quantifiers `* + ? {m} {m,} {m,n} {,n}`, and the global
//! flags `(?i)`, `(?m)`, `(?s)` at the start. `\w`, `\d` and `\s` are the
//! Unicode classes of Python `str` patterns (approximated by the Rust
//! character properties). Backreferences, lookbehinds, scoped flags and
//! possessive quantifiers are refused with an error rather than matched
//! differently.

use std::fmt;

/// Why a pattern was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegexError(pub String);

impl fmt::Display for RegexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RegexError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Perl {
    Digit,
    Word,
    Space,
}

#[derive(Debug, Clone, PartialEq)]
enum ClassItem {
    Range(char, char),
    Perl(Perl, bool),
}

#[derive(Debug, Clone, PartialEq)]
struct Class {
    items: Vec<ClassItem>,
    negated: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Look {
    /// `\A`, or `^` without `(?m)`.
    Start,
    /// `^` with `(?m)`.
    LineStart,
    /// `$` without `(?m)`: the end, or before a final newline.
    End,
    /// `$` with `(?m)`.
    LineEnd,
    /// `\Z`.
    AbsoluteEnd,
    WordBoundary,
    NotWordBoundary,
}

#[derive(Debug, Clone, PartialEq)]
enum Node {
    Empty,
    Char(char),
    /// `.`; `true` also matches a newline (`(?s)`).
    Any(bool),
    Class(Class),
    Look(Look),
    Lookahead(Box<Node>, bool),
    Concat(Vec<Node>),
    Alt(Vec<Node>),
    Repeat {
        node: Box<Node>,
        min: u32,
        max: Option<u32>,
        greedy: bool,
    },
}

#[derive(Debug, Clone, Copy, Default)]
struct Flags {
    ignore_case: bool,
    multiline: bool,
    dot_all: bool,
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    flags: Flags,
    pattern: &'a str,
}

fn is_word(c: char) -> bool {
    c == '_' || c.is_alphanumeric()
}

fn is_space(c: char) -> bool {
    c.is_whitespace() || ('\x1c'..='\x1f').contains(&c)
}

fn perl_matches(p: Perl, c: char) -> bool {
    match p {
        Perl::Digit => c.is_numeric() && (c.is_ascii_digit() || !c.is_ascii()),
        Perl::Word => is_word(c),
        Perl::Space => is_space(c),
    }
}

impl Parser<'_> {
    fn error<T>(&self, message: &str) -> Result<T, RegexError> {
        Err(RegexError(format!("{message} at position {} of {:?}", self.pos, self.pattern)))
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn eat(&mut self, c: char) -> bool {
        if self.peek() == Some(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn alternation(&mut self, depth: usize) -> Result<Node, RegexError> {
        if depth > 200 {
            return self.error("pattern nests too deeply");
        }
        let mut branches = vec![self.sequence(depth)?];
        while self.eat('|') {
            branches.push(self.sequence(depth)?);
        }
        Ok(if branches.len() == 1 { branches.pop().unwrap_or(Node::Empty) } else { Node::Alt(branches) })
    }

    fn sequence(&mut self, depth: usize) -> Result<Node, RegexError> {
        let mut items = Vec::new();
        while let Some(c) = self.peek() {
            if c == '|' || c == ')' {
                break;
            }
            let atom = self.atom(depth)?;
            let atom = self.quantified(atom)?;
            items.push(atom);
        }
        Ok(match items.len() {
            0 => Node::Empty,
            1 => items.pop().unwrap_or(Node::Empty),
            _ => Node::Concat(items),
        })
    }

    /// `{m}`, `{m,}`, `{m,n}` or `{,n}` at the current `{`, or `None` when
    /// it is a literal brace (as in Python).
    fn braces(&mut self) -> Option<(u32, Option<u32>)> {
        let start = self.pos;
        let mut text = String::new();
        self.pos += 1;
        while let Some(c) = self.peek() {
            self.pos += 1;
            if c == '}' {
                let (min, max) = match text.split_once(',') {
                    None => {
                        let n = text.parse().ok()?;
                        (n, Some(n))
                    }
                    Some((a, b)) => {
                        let min = if a.is_empty() { 0 } else { a.parse().ok()? };
                        let max = if b.is_empty() { None } else { Some(b.parse().ok()?) };
                        (min, max)
                    }
                };
                if text.is_empty() || !text.chars().all(|c| c.is_ascii_digit() || c == ',') || text.matches(',').count() > 1 {
                    break;
                }
                return Some((min, max));
            }
            text.push(c);
        }
        self.pos = start;
        None
    }

    fn quantified(&mut self, atom: Node) -> Result<Node, RegexError> {
        let mut node = atom;
        let mut repeated = false;
        loop {
            let (min, max) = match self.peek() {
                Some('*') => {
                    self.pos += 1;
                    (0, None)
                }
                Some('+') => {
                    self.pos += 1;
                    (1, None)
                }
                Some('?') => {
                    self.pos += 1;
                    (0, Some(1))
                }
                Some('{') => match self.braces() {
                    Some(q) => q,
                    None => return Ok(node),
                },
                _ => return Ok(node),
            };
            if matches!(node, Node::Look(_) | Node::Lookahead(..) | Node::Empty) {
                return self.error("nothing to repeat");
            }
            if repeated {
                return self.error("multiple repeat");
            }
            repeated = true;
            if max.is_some_and(|m| m < min) {
                return self.error("min repeat greater than max repeat");
            }
            if min > 1000 || max.is_some_and(|m| m > 1000) {
                return self.error("repeat counts over 1000 are not supported");
            }
            let greedy = !self.eat('?');
            if self.peek() == Some('+') {
                return self.error("possessive quantifiers are not supported");
            }
            node = Node::Repeat { node: Box::new(node), min, max, greedy };
        }
    }

    fn hex(&mut self, digits: usize) -> Result<char, RegexError> {
        let s: String = self.chars.get(self.pos..self.pos + digits).map(|c| c.iter().collect()).unwrap_or_default();
        if s.len() != digits || !s.chars().all(|c| c.is_ascii_hexdigit()) {
            return self.error("incomplete escape");
        }
        self.pos += digits;
        let v = u32::from_str_radix(&s, 16).map_err(|_| RegexError("bad escape".into()))?;
        char::from_u32(v).map_or_else(|| self.error("bad escape"), Ok)
    }

    /// An escape after `\`: a character, a Perl class, or an anchor.
    fn escape(&mut self, in_class: bool) -> Result<Node, RegexError> {
        let Some(c) = self.peek() else { return self.error("bad escape (end of pattern)") };
        self.pos += 1;
        let ch = |c| Ok(Node::Char(c));
        let perl = |p, neg| Ok(Node::Class(Class { items: vec![ClassItem::Perl(p, neg)], negated: false }));
        match c {
            'n' => ch('\n'),
            't' => ch('\t'),
            'r' => ch('\r'),
            'f' => ch('\x0c'),
            'v' => ch('\x0b'),
            'a' => ch('\x07'),
            'b' if in_class => ch('\x08'),
            '0' => ch('\0'),
            'x' => Ok(Node::Char(self.hex(2)?)),
            'u' => Ok(Node::Char(self.hex(4)?)),
            'U' => Ok(Node::Char(self.hex(8)?)),
            'd' => perl(Perl::Digit, false),
            'D' => perl(Perl::Digit, true),
            'w' => perl(Perl::Word, false),
            'W' => perl(Perl::Word, true),
            's' => perl(Perl::Space, false),
            'S' => perl(Perl::Space, true),
            'A' if !in_class => Ok(Node::Look(Look::Start)),
            'Z' if !in_class => Ok(Node::Look(Look::AbsoluteEnd)),
            'b' => Ok(Node::Look(Look::WordBoundary)),
            'B' if !in_class => Ok(Node::Look(Look::NotWordBoundary)),
            '1'..='9' => self.error("backreferences are not supported"),
            c if c.is_ascii_alphanumeric() => self.error(&format!("bad escape \\{c}")),
            c => ch(c),
        }
    }

    fn class(&mut self) -> Result<Node, RegexError> {
        let negated = self.eat('^');
        let mut items = Vec::new();
        let mut first = true;
        loop {
            let Some(c) = self.peek() else { return self.error("unterminated character set") };
            if c == ']' && !first {
                self.pos += 1;
                break;
            }
            first = false;
            self.pos += 1;
            let lo = if c == '\\' {
                match self.escape(true)? {
                    Node::Char(c) => c,
                    Node::Class(Class { items: mut inner, .. }) => {
                        items.append(&mut inner);
                        continue;
                    }
                    _ => return self.error("bad escape in a character set"),
                }
            } else {
                c
            };
            let range = self.peek() == Some('-') && self.chars.get(self.pos + 1).is_some_and(|c| *c != ']');
            if range {
                self.pos += 1;
                let Some(h) = self.peek() else { return self.error("unterminated character set") };
                self.pos += 1;
                let hi = if h == '\\' {
                    match self.escape(true)? {
                        Node::Char(c) => c,
                        _ => return self.error("bad character range"),
                    }
                } else {
                    h
                };
                if hi < lo {
                    return self.error("bad character range");
                }
                items.push(ClassItem::Range(lo, hi));
            } else {
                items.push(ClassItem::Range(lo, lo));
            }
        }
        Ok(Node::Class(Class { items, negated }))
    }

    fn atom(&mut self, depth: usize) -> Result<Node, RegexError> {
        let Some(c) = self.peek() else { return Ok(Node::Empty) };
        self.pos += 1;
        match c {
            '(' => {
                let node = if self.eat('?') {
                    match self.peek() {
                        Some(':') => {
                            self.pos += 1;
                            self.alternation(depth + 1)?
                        }
                        Some('P') => {
                            self.pos += 1;
                            if !self.eat('<') {
                                return self.error("named backreferences are not supported");
                            }
                            let name_start = self.pos;
                            while self.peek().is_some_and(|c| c != '>') {
                                self.pos += 1;
                            }
                            let name: String = self.chars[name_start..self.pos].iter().collect();
                            if !self.eat('>') || name.is_empty() || !name.chars().all(is_word) {
                                return self.error("bad group name");
                            }
                            self.alternation(depth + 1)?
                        }
                        Some('=') | Some('!') => {
                            let negate = self.peek() == Some('!');
                            self.pos += 1;
                            Node::Lookahead(Box::new(self.alternation(depth + 1)?), negate)
                        }
                        Some('#') => {
                            while self.peek().is_some_and(|c| c != ')') {
                                self.pos += 1;
                            }
                            if !self.eat(')') {
                                return self.error("missing ), unterminated comment");
                            }
                            return Ok(Node::Empty);
                        }
                        Some('<') => return self.error("lookbehinds are not supported"),
                        _ => return self.error("only (?:, (?P<name>, (?=, (?! and leading (?imsx) groups are supported"),
                    }
                } else {
                    self.alternation(depth + 1)?
                };
                if !self.eat(')') {
                    return self.error("missing ), unterminated subpattern");
                }
                Ok(node)
            }
            ')' => self.error("unbalanced parenthesis"),
            '[' => self.class(),
            '.' => Ok(Node::Any(self.flags.dot_all)),
            '^' => Ok(Node::Look(if self.flags.multiline { Look::LineStart } else { Look::Start })),
            '$' => Ok(Node::Look(if self.flags.multiline { Look::LineEnd } else { Look::End })),
            '\\' => self.escape(false),
            '*' | '+' | '?' => self.error("nothing to repeat"),
            '{' => {
                self.pos -= 1;
                if self.braces().is_some() {
                    return self.error("nothing to repeat");
                }
                self.pos += 1;
                Ok(Node::Char('{'))
            }
            c => Ok(Node::Char(c)),
        }
    }
}

#[derive(Debug, Clone)]
enum Inst {
    Char(char),
    Any(bool),
    Class(Class),
    Look(Look),
    /// Lookahead sub-program and whether it is negated.
    Lookahead(usize, bool),
    Split(usize, usize),
    Jmp(usize),
    Match,
}

struct Compiler {
    progs: Vec<Vec<Inst>>,
}

const MAX_INSTRUCTIONS: usize = 100_000;

impl Compiler {
    fn emit(&mut self, p: usize, inst: Inst) -> Result<usize, RegexError> {
        let prog = &mut self.progs[p];
        if prog.len() >= MAX_INSTRUCTIONS {
            return Err(RegexError("pattern is too large".into()));
        }
        prog.push(inst);
        Ok(prog.len() - 1)
    }

    fn compile(&mut self, p: usize, node: &Node) -> Result<(), RegexError> {
        match node {
            Node::Empty => {}
            Node::Char(c) => {
                self.emit(p, Inst::Char(*c))?;
            }
            Node::Any(nl) => {
                self.emit(p, Inst::Any(*nl))?;
            }
            Node::Class(c) => {
                self.emit(p, Inst::Class(c.clone()))?;
            }
            Node::Look(l) => {
                self.emit(p, Inst::Look(*l))?;
            }
            Node::Lookahead(inner, negate) => {
                self.progs.push(Vec::new());
                let sub = self.progs.len() - 1;
                self.compile(sub, inner)?;
                self.emit(sub, Inst::Match)?;
                self.emit(p, Inst::Lookahead(sub, *negate))?;
            }
            Node::Concat(items) => {
                for n in items {
                    self.compile(p, n)?;
                }
            }
            Node::Alt(branches) => {
                let mut jumps = Vec::new();
                for (k, b) in branches.iter().enumerate() {
                    if k + 1 < branches.len() {
                        let split = self.emit(p, Inst::Split(0, 0))?;
                        self.compile(p, b)?;
                        jumps.push(self.emit(p, Inst::Jmp(0))?);
                        let next = self.progs[p].len();
                        self.progs[p][split] = Inst::Split(split + 1, next);
                    } else {
                        self.compile(p, b)?;
                    }
                }
                let end = self.progs[p].len();
                for j in jumps {
                    self.progs[p][j] = Inst::Jmp(end);
                }
            }
            Node::Repeat { node, min, max, greedy } => {
                for _ in 0..*min {
                    self.compile(p, node)?;
                }
                let split = |a: usize, b: usize| if *greedy { Inst::Split(a, b) } else { Inst::Split(b, a) };
                match max {
                    None => {
                        let start = self.emit(p, Inst::Split(0, 0))?;
                        self.compile(p, node)?;
                        self.emit(p, Inst::Jmp(start))?;
                        let end = self.progs[p].len();
                        self.progs[p][start] = split(start + 1, end);
                    }
                    Some(max) => {
                        let mut holes = Vec::new();
                        for _ in *min..*max {
                            holes.push(self.emit(p, Inst::Split(0, 0))?);
                            self.compile(p, node)?;
                        }
                        let end = self.progs[p].len();
                        for h in holes {
                            self.progs[p][h] = split(h + 1, end);
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// A compiled pattern.
#[derive(Debug, Clone)]
pub struct Regex {
    progs: Vec<Vec<Inst>>,
    ignore_case: bool,
}

fn fold(c: char) -> char {
    let mut l = c.to_lowercase();
    match (l.next(), l.next()) {
        (Some(x), None) => x,
        _ => c,
    }
}

impl Regex {
    /// Compile a pattern of the supported subset.
    pub fn new(pattern: &str) -> Result<Regex, RegexError> {
        let mut parser = Parser { chars: pattern.chars().collect(), pos: 0, flags: Flags::default(), pattern };
        // Global flags, which Python accepts only at the start.
        while parser.chars.get(parser.pos) == Some(&'(') && parser.chars.get(parser.pos + 1) == Some(&'?') {
            let mut k = parser.pos + 2;
            let mut flags = parser.flags;
            let mut any = false;
            while let Some(&c) = parser.chars.get(k) {
                match c {
                    'i' => flags.ignore_case = true,
                    'm' => flags.multiline = true,
                    's' => flags.dot_all = true,
                    'u' => {}
                    _ => break,
                }
                any = true;
                k += 1;
            }
            if !any {
                break;
            }
            match parser.chars.get(k) {
                Some(')') => {
                    parser.flags = flags;
                    parser.pos = k + 1;
                }
                Some(':') => return parser.error("scoped inline flags are not supported"),
                _ => return parser.error("unknown flag"),
            }
        }
        let node = parser.alternation(0)?;
        if parser.pos < parser.chars.len() {
            return parser.error("unbalanced parenthesis");
        }
        let mut c = Compiler { progs: vec![Vec::new()] };
        c.compile(0, &node)?;
        c.emit(0, Inst::Match)?;
        Ok(Regex { progs: c.progs, ignore_case: parser.flags.ignore_case })
    }

    fn class_matches(&self, class: &Class, c: char) -> bool {
        let test = |c: char| {
            class.items.iter().any(|item| match item {
                ClassItem::Range(lo, hi) => (*lo..=*hi).contains(&c),
                ClassItem::Perl(p, neg) => perl_matches(*p, c) != *neg,
            })
        };
        let mut hit = test(c);
        if !hit && self.ignore_case {
            hit = test(fold(c)) || c.to_uppercase().any(test);
        }
        hit != class.negated
    }

    fn look(&self, l: Look, text: &[char], at: usize) -> bool {
        let n = text.len();
        let word = |i: usize| i < n && is_word(text[i]);
        let before = at > 0 && word(at - 1);
        match l {
            Look::Start => at == 0,
            Look::LineStart => at == 0 || text[at - 1] == '\n',
            Look::End => at == n || (at + 1 == n && text[at] == '\n'),
            Look::LineEnd => at == n || text[at] == '\n',
            Look::AbsoluteEnd => at == n,
            Look::WordBoundary => before != word(at),
            Look::NotWordBoundary => before == word(at),
        }
    }

    /// Add a thread at `pc` to `list` (with its start), following jumps,
    /// splits and zero-width assertions at `at`, in priority order.
    #[allow(clippy::too_many_arguments)]
    fn add(
        &self,
        prog: usize,
        list: &mut Vec<(usize, usize)>,
        seen: &mut [bool],
        pc: usize,
        start: usize,
        text: &[char],
        at: usize,
    ) {
        let mut stack = vec![pc];
        while let Some(pc) = stack.pop() {
            if seen[pc] {
                continue;
            }
            seen[pc] = true;
            match &self.progs[prog][pc] {
                Inst::Jmp(t) => stack.push(*t),
                Inst::Split(a, b) => {
                    stack.push(*b);
                    stack.push(*a);
                }
                Inst::Look(l) => {
                    if self.look(*l, text, at) {
                        stack.push(pc + 1);
                    }
                }
                Inst::Lookahead(sub, negate) => {
                    if self.run(*sub, text, at, true, None).is_some() != *negate {
                        stack.push(pc + 1);
                    }
                }
                _ => list.push((pc, start)),
            }
        }
    }

    /// The leftmost-first match of program `prog` starting at `from` or
    /// later (only at `from` when `anchored`), skipping an empty match that
    /// ends at `forbid_empty_at`.
    fn run(
        &self,
        prog: usize,
        text: &[char],
        from: usize,
        anchored: bool,
        forbid_empty_at: Option<usize>,
    ) -> Option<(usize, usize)> {
        let insts = &self.progs[prog];
        let mut clist: Vec<(usize, usize)> = Vec::new();
        let mut nlist: Vec<(usize, usize)> = Vec::new();
        let mut seen = vec![false; insts.len()];
        let mut found: Option<(usize, usize)> = None;
        let mut at = from;
        loop {
            if found.is_none() && (at == from || !anchored) {
                self.add(prog, &mut clist, &mut seen, 0, at, text, at);
            }
            if clist.is_empty() && (anchored || found.is_some() || at >= text.len()) {
                break;
            }
            seen.iter_mut().for_each(|s| *s = false);
            let c = text.get(at).copied();
            for &(pc, start) in clist.iter() {
                let step = match &insts[pc] {
                    Inst::Match => {
                        if forbid_empty_at == Some(at) && start == at {
                            continue;
                        }
                        found = Some((start, at));
                        break;
                    }
                    Inst::Char(x) => c.is_some_and(|c| c == *x || (self.ignore_case && fold(c) == fold(*x))),
                    Inst::Any(nl) => c.is_some_and(|c| *nl || c != '\n'),
                    Inst::Class(class) => c.is_some_and(|c| self.class_matches(class, c)),
                    _ => false,
                };
                if step {
                    self.add(prog, &mut nlist, &mut seen, pc + 1, start, text, at + 1);
                }
            }
            if at >= text.len() {
                // Threads that wait for a character cannot advance; only a
                // Match can still be taken, which the loop above did.
                break;
            }
            std::mem::swap(&mut clist, &mut nlist);
            nlist.clear();
            at += 1;
            if found.is_some() && clist.is_empty() {
                break;
            }
        }
        found
    }

    /// Whether the pattern matches anywhere in `text`.
    pub fn is_match(&self, text: &str) -> bool {
        let chars: Vec<char> = text.chars().collect();
        self.run(0, &chars, 0, false, None).is_some()
    }

    /// Replace every match by `replacement` (literally) and count them, as
    /// Python's `pattern.subn(lambda _: replacement, text)`: matches do not
    /// overlap, and an empty match is replaced except right after an empty
    /// match at the same place.
    pub fn replace_all(&self, text: &str, replacement: &str) -> (String, usize) {
        let chars: Vec<char> = text.chars().collect();
        let mut out = String::with_capacity(text.len());
        let mut copied = 0usize;
        let mut search = 0usize;
        let mut forbid = None;
        let mut count = 0usize;
        while search <= chars.len() {
            let Some((s, e)) = self.run(0, &chars, search, false, forbid) else { break };
            out.extend(&chars[copied..s]);
            out.push_str(replacement);
            count += 1;
            copied = e;
            search = e;
            forbid = (s == e).then_some(e);
        }
        out.extend(&chars[copied.min(chars.len())..]);
        (out, count)
    }
}

#[cfg(test)]
mod tests {
    use super::Regex;

    #[test]
    fn python_semantics() {
        let sub = |p: &str, t: &str| Regex::new(p).unwrap().replace_all(t, "-");
        assert_eq!(sub("x*", "abxd"), ("-a-b--d-".into(), 5));
        assert_eq!(sub("a|ab", "ab"), ("-b".into(), 1));
        assert_eq!(sub("a+?", "aaa"), ("---".into(), 3));
        assert_eq!(sub(r"\bfoo\b", "foo food foo"), ("- food -".into(), 2));
        assert_eq!(sub("(?i)HOME", "home Home"), ("- -".into(), 2));
        assert_eq!(sub("$", "a\n"), ("a-\n-".into(), 2));
        assert_eq!(sub(r"/home/[^/\s]+", "/home/USER/x /home/runner"), ("-/x -".into(), 2));
        assert_eq!(sub(r"a(?=b)", "ab ac"), ("-b ac".into(), 1));
        assert_eq!(sub(r"\d{2,3}", "1 12 1234"), ("1 - -4".into(), 2));
        assert_eq!(sub("(a*)*", "b"), ("-b-".into(), 2));
        assert!(Regex::new(r"(a)\1").is_err());
        assert!(Regex::new("(?<=a)b").is_err());
        assert!(Regex::new("*a").is_err());
        assert!(Regex::new("a{2").is_ok());
    }
}
