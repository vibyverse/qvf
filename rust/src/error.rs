//! Finding codes (spec Appendix B), errors and findings.

use std::fmt;

/// A stable finding code of spec Appendix B. `E_...` codes are errors and
/// `W_...` codes are warnings. A code is never renamed or reused; the test
/// `codes_match_appendix_b` keeps this list equal to the generated table.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Code {
    E_ZIP_UNREADABLE,
    E_ZIP_DUPLICATE_ENTRY,
    E_MEMBER_PATH,
    E_MEMBER_TOO_LARGE,
    E_TOTAL_TOO_LARGE,
    E_COMPRESSION_RATIO,
    W_ORPHAN_MEMBER,
    E_MANIFEST_MISSING,
    E_MANIFEST_JSON,
    E_JSON_TOO_LARGE,
    E_JSON_TOO_DEEP,
    E_SCHEMA,
    W_SCHEMA_SKIPPED,
    E_VERSION,
    W_VERSION_WITHDRAWN,
    E_SOURCE,
    E_SECTION_MALFORMED,
    E_SECTION_ID_DUPLICATE,
    E_RESERVED_ID,
    E_KIND_UNKNOWN,
    E_CRITICAL_UNDECLARED,
    E_MEMBER_MALFORMED,
    E_MEMBER_MISSING,
    E_SHA256_MISMATCH,
    E_JSON_INVALID,
    E_DTYPE,
    E_BINARY_SIZE,
    E_UTF8,
    E_PAYLOAD_SCHEMA,
    E_STRUCTURE_PERIODICITY,
    E_CHARGE_MULTIPLICITY,
    E_REF_UNRESOLVED,
    E_FILES_INDEX,
    W_UNIT_UNKNOWN,
    E_ROOT_INCONSISTENT,
    E_NORMALIZATION,
    E_UNIT_CONTRADICTED,
}

impl Code {
    /// Every code, in Appendix B order.
    pub const ALL: &'static [Code] = &[
        Code::E_ZIP_UNREADABLE,
        Code::E_ZIP_DUPLICATE_ENTRY,
        Code::E_MEMBER_PATH,
        Code::E_MEMBER_TOO_LARGE,
        Code::E_TOTAL_TOO_LARGE,
        Code::E_COMPRESSION_RATIO,
        Code::W_ORPHAN_MEMBER,
        Code::E_MANIFEST_MISSING,
        Code::E_MANIFEST_JSON,
        Code::E_JSON_TOO_LARGE,
        Code::E_JSON_TOO_DEEP,
        Code::E_SCHEMA,
        Code::W_SCHEMA_SKIPPED,
        Code::E_VERSION,
        Code::W_VERSION_WITHDRAWN,
        Code::E_SOURCE,
        Code::E_SECTION_MALFORMED,
        Code::E_SECTION_ID_DUPLICATE,
        Code::E_RESERVED_ID,
        Code::E_KIND_UNKNOWN,
        Code::E_CRITICAL_UNDECLARED,
        Code::E_MEMBER_MALFORMED,
        Code::E_MEMBER_MISSING,
        Code::E_SHA256_MISMATCH,
        Code::E_JSON_INVALID,
        Code::E_DTYPE,
        Code::E_BINARY_SIZE,
        Code::E_UTF8,
        Code::E_PAYLOAD_SCHEMA,
        Code::E_STRUCTURE_PERIODICITY,
        Code::E_CHARGE_MULTIPLICITY,
        Code::E_REF_UNRESOLVED,
        Code::E_FILES_INDEX,
        Code::W_UNIT_UNKNOWN,
        Code::E_ROOT_INCONSISTENT,
        Code::E_NORMALIZATION,
        Code::E_UNIT_CONTRADICTED,
    ];

    /// The code as Appendix B spells it.
    pub fn as_str(self) -> &'static str {
        match self {
            Code::E_ZIP_UNREADABLE => "E_ZIP_UNREADABLE",
            Code::E_ZIP_DUPLICATE_ENTRY => "E_ZIP_DUPLICATE_ENTRY",
            Code::E_MEMBER_PATH => "E_MEMBER_PATH",
            Code::E_MEMBER_TOO_LARGE => "E_MEMBER_TOO_LARGE",
            Code::E_TOTAL_TOO_LARGE => "E_TOTAL_TOO_LARGE",
            Code::E_COMPRESSION_RATIO => "E_COMPRESSION_RATIO",
            Code::W_ORPHAN_MEMBER => "W_ORPHAN_MEMBER",
            Code::E_MANIFEST_MISSING => "E_MANIFEST_MISSING",
            Code::E_MANIFEST_JSON => "E_MANIFEST_JSON",
            Code::E_JSON_TOO_LARGE => "E_JSON_TOO_LARGE",
            Code::E_JSON_TOO_DEEP => "E_JSON_TOO_DEEP",
            Code::E_SCHEMA => "E_SCHEMA",
            Code::W_SCHEMA_SKIPPED => "W_SCHEMA_SKIPPED",
            Code::E_VERSION => "E_VERSION",
            Code::W_VERSION_WITHDRAWN => "W_VERSION_WITHDRAWN",
            Code::E_SOURCE => "E_SOURCE",
            Code::E_SECTION_MALFORMED => "E_SECTION_MALFORMED",
            Code::E_SECTION_ID_DUPLICATE => "E_SECTION_ID_DUPLICATE",
            Code::E_RESERVED_ID => "E_RESERVED_ID",
            Code::E_KIND_UNKNOWN => "E_KIND_UNKNOWN",
            Code::E_CRITICAL_UNDECLARED => "E_CRITICAL_UNDECLARED",
            Code::E_MEMBER_MALFORMED => "E_MEMBER_MALFORMED",
            Code::E_MEMBER_MISSING => "E_MEMBER_MISSING",
            Code::E_SHA256_MISMATCH => "E_SHA256_MISMATCH",
            Code::E_JSON_INVALID => "E_JSON_INVALID",
            Code::E_DTYPE => "E_DTYPE",
            Code::E_BINARY_SIZE => "E_BINARY_SIZE",
            Code::E_UTF8 => "E_UTF8",
            Code::E_PAYLOAD_SCHEMA => "E_PAYLOAD_SCHEMA",
            Code::E_STRUCTURE_PERIODICITY => "E_STRUCTURE_PERIODICITY",
            Code::E_CHARGE_MULTIPLICITY => "E_CHARGE_MULTIPLICITY",
            Code::E_REF_UNRESOLVED => "E_REF_UNRESOLVED",
            Code::E_FILES_INDEX => "E_FILES_INDEX",
            Code::W_UNIT_UNKNOWN => "W_UNIT_UNKNOWN",
            Code::E_ROOT_INCONSISTENT => "E_ROOT_INCONSISTENT",
            Code::E_NORMALIZATION => "E_NORMALIZATION",
            Code::E_UNIT_CONTRADICTED => "E_UNIT_CONTRADICTED",
        }
    }

    /// The severity Appendix B gives the code.
    pub fn severity(self) -> Severity {
        if self.as_str().starts_with('E') {
            Severity::Error
        } else {
            Severity::Warning
        }
    }

    /// The code for an Appendix B spelling.
    pub fn parse(s: &str) -> Option<Code> {
        Code::ALL.iter().copied().find(|c| c.as_str() == s)
    }
}

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Severity of a finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Severity {
    Error,
    Warning,
}

impl Severity {
    /// `"error"` or `"warning"`, as Appendix B and the corpus spell it.
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Error => "error",
            Severity::Warning => "warning",
        }
    }
}

/// A failure of the reader, with the finding code that describes it.
#[derive(Debug, Clone)]
pub struct QvfError {
    pub code: Code,
    pub message: String,
}

impl QvfError {
    pub fn new(code: Code, message: impl Into<String>) -> Self {
        QvfError { code, message: message.into() }
    }
}

impl fmt::Display for QvfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

impl std::error::Error for QvfError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_match_appendix_b() {
        let generated: Vec<(&str, &str)> = crate::generated::FINDING_CODES.to_vec();
        let ours: Vec<(&str, &str)> = Code::ALL.iter().map(|c| (c.as_str(), c.severity().as_str())).collect();
        assert_eq!(ours, generated);
    }
}
