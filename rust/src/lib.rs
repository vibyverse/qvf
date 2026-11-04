//! An independent Rust implementation of QVF, the archive format of the
//! vibe-qc family: a reader with typed little-endian tensors and a validator
//! that reports the stable finding codes of the format specification
//! (spec § 6, Appendix B), certified by the shared conformance corpus
//! (`conformance/`).
//!
//! ```no_run
//! let report = qvf::validate_file("water.qvf", &qvf::ValidateOptions::default());
//! for f in &report.findings {
//!     println!("{} {} {}", f.code, f.location, f.message);
//! }
//! let reader = qvf::Reader::open("water.qvf")?;
//! for s in reader.sections() {
//!     println!("{} {}", s.id, s.kind);
//! }
//! # Ok::<(), qvf::QvfError>(())
//! ```
//!
//! The crate does not depend on the C++ library or its C ABI. The JSON
//! Schema check runs on the instruction table `scripts/gen_kinds.py`
//! compiles from `spec/qvf_manifest.schema.json` (`src/generated.rs`).

#![forbid(unsafe_code)]

pub mod error;
#[rustfmt::skip]
pub mod generated;
pub mod json;
pub mod limits;
pub mod reader;
pub mod schema_check;
pub mod validate;
pub mod zip;

pub use error::{Code, QvfError, Severity};
pub use limits::Limits;
pub use reader::{DType, Reader, Section, Tensor, TensorData};
pub use validate::{validate_bytes, validate_file, Finding, Level, Report, ValidateOptions};
