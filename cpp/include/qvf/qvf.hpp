// qvf/qvf.hpp — public API for the zero-dependency QVF writer library.
//
// QVF ("Quantum Visualization Format") is a ZIP container with a mandatory
// manifest.json plus typed JSON / binary members. This library lets any C++
// code emit conforming .qvf archives. It depends only on the C++17 standard
// library — ZIP, CRC-32, SHA-256 and JSON emission are all self-contained.
//
// See spec/qvf-format-spec.md for the format; docs/library_guide.md for usage.
// License: Apache-2.0.
#ifndef QVF_QVF_HPP
#define QVF_QVF_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "qvf/container.hpp"
#include "qvf/json.hpp"
#include "qvf/model.hpp"
#include "qvf/units.hpp"

namespace qvf {

// ---------------------------------------------------------------------------
// Value types
// ---------------------------------------------------------------------------

// Source (program, version, calculation) is defined in qvf/model.hpp.

struct Atom {
    std::string symbol;
    int atomic_number = 0;
    std::array<double, 3> position{{0, 0, 0}};  // Ångström
};

struct Bond {
    int i = 0;
    int j = 0;
    double order = 1.0;
};

// Binary member dtypes permitted by the schema.
enum class DType {
    Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64, Float32, Float64
};

// A little-endian, C-contiguous binary payload. Build via the tensor_* helpers,
// which serialize host-endian-independently.
struct Tensor {
    DType dtype = DType::Float64;
    std::vector<int64_t> shape;
    std::vector<uint8_t> bytes;
};

Tensor tensor_f64(const std::vector<double>& v, std::vector<int64_t> shape);
Tensor tensor_f32(const std::vector<double>& v, std::vector<int64_t> shape);
Tensor tensor_i32(const std::vector<int32_t>& v, std::vector<int64_t> shape);
Tensor tensor_i64(const std::vector<int64_t>& v, std::vector<int64_t> shape);

// Row-major dense matrix of doubles (rows are MOs for wavefunction.gto).
struct Matrix {
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<double> data;  // size rows*cols, row-major
};

struct Shell {
    int center = 0;                     // 0-based atom index
    int l = 0;                          // shell angular momentum
    std::vector<double> exponents;      // bohr^-2
    std::vector<double> coefficients;   // on N_i-normalized primitives
                                        // (spec Appendix A.1)
    bool pure = true;                   // true = spherical, false = Cartesian
};

struct WavefunctionGTO {
    std::vector<Shell> shells;
    std::string structure_ref = "structure";
    std::string orbital_kind = "canonical";
    bool pure = true;

    // Restricted: fill mo_coefficients. Unrestricted: set unrestricted=true and
    // fill the alpha/beta matrices instead.
    bool unrestricted = false;
    // The mo_metadata `spin` (QVF spec § 5.3). Empty means "restricted", or
    // "unrestricted" when `unrestricted` is set. "restricted_open" and
    // "generalized" fill mo_coefficients; "generalized" (two-component) needs
    // complex_coefficients and 2 n_ao columns: the α components of the AOs,
    // then their β components.
    std::string spin;
    Matrix mo_coefficients;                 // [n_mo, n_ao]
    std::vector<double> energies;           // Hartree
    std::vector<double> occupations;
    Matrix mo_coefficients_alpha;
    Matrix mo_coefficients_beta;
    std::vector<double> energies_alpha, occupations_alpha;
    std::vector<double> energies_beta, occupations_beta;
    // Complex coefficients (QVF spec § 5.3): every coefficient Matrix then
    // holds rows * cols * 2 doubles, (re, im) interleaved, and is written as
    // float64 [rows, cols, 2]; mo_metadata gets `complex: true`.
    bool complex_coefficients = false;

    // Optional Γ-point tag.
    bool has_kpoint = false;
    std::array<double, 3> kpoint{{0, 0, 0}};

    // Divide each coefficient by the primitive norm N_i before writing (set
    // this when your engine stores libint/libcint/GBW-style coefficients that
    // are pre-multiplied by N_i). See spec Appendix A.
    bool coeffs_are_libint_normalized = false;

    // --- optional, QVF spec § 5.3 ---
    // Effective core potentials: a JSON array of {center, n_core_electrons,
    // label?, projectors?: [{l, exponents, coefficients, powers, local?}]}, at
    // most one per atom; Json() omits it. A Json::raw fragment is parsed.
    Json ecp;
    // Electrons the wavefunction describes (Σ Z − Σ n_core_electrons −
    // charge); < 0 omits it. The written occupations must sum to it.
    int n_electrons = -1;
    // One irreducible-representation label per MO; empty omits them. The
    // single-set labels, or the α and β labels of an unrestricted wavefunction.
    std::vector<std::string> symmetry_labels;
    std::vector<std::string> symmetry_labels_alpha, symmetry_labels_beta;
    // AO-basis matrices, each [n_ao, n_ao] in the AO order of the basis as
    // written; an empty Matrix omits it. `overlap` is the AO overlap matrix,
    // which validators check against Appendix A.1 (producers SHOULD write
    // it); `density_matrix` (with mo_coefficients) is the spin-summed density
    // matrix, `density_matrix_alpha` and `density_matrix_beta` (with the
    // alpha and beta coefficients, both or neither) the per-spin ones.
    Matrix overlap;
    Matrix density_matrix;
    Matrix density_matrix_alpha;
    Matrix density_matrix_beta;

    std::string section_id = "wf";

    // Natural transition orbitals (QVF spec § 5.3) are two sections per
    // state, with orbital_kind "nto": nto_role "hole" or "particle",
    // state_ref the id of the excited_states section and state_index (>= 1)
    // the state's index there. Empty / 0 omits a field; state_index needs
    // state_ref, nto_role needs orbital_kind "nto". An "nto" wavefunction is
    // also written with occupation_semantics "transition_weight" (its
    // occupations are the pair weights). state_ref is not resolved here,
    // because the excited_states section is usually added after its orbitals.
    std::string state_ref;
    int state_index = 0;
    // Unrestricted NTOs still use two sections. Alpha/beta row counts may
    // differ; each channel's weights are nonnegative, non-increasing and
    // identical to its partner. The combined spin sum is one only for a
    // normalized state; full TDDFT/RPA sums may differ. Missing energies are
    // zero-filled. The validator checks pairs and forward state references.
    std::string nto_role;
    // Minimal reference basis for localized orbitals; absent omits the field,
    // an explicitly empty string is rejected.
    std::optional<std::string> localization_reference_basis;
};

// k-resolved Bloch orbitals of a periodic system (kind wavefunction.bloch,
// QVF spec § 5.3.1, provisional): complex MO coefficients at every k-point of
// a set over the atom-centred Gaussian basis of the cell, with the Bloch sum
// phi_mu,k(r) = N^-1/2 sum_T exp(+i k.T) chi_mu(r - T). The basis is the one
// of WavefunctionGTO (shells, pure, coeffs_are_libint_normalized); n_ao is
// derived from the shells, n_k from kpoints and n_spin from spin.
struct WavefunctionBloch {
    std::vector<Shell> shells;
    std::string structure_ref = "structure";   // must name a periodic structure
    std::string orbital_kind = "canonical";
    bool pure = true;
    bool coeffs_are_libint_normalized = false;
    std::string spin = "restricted";           // "restricted" (n_spin 1) or
                                               // "unrestricted" (n_spin 2)
    // k-points in fractional reciprocal-lattice coordinates.
    std::vector<std::array<double, 3>> kpoints;
    // Optional kpoints payload fields; an empty vector or a false has_* omits.
    std::vector<double> weights;               // n_k, >= 0, summing to 1 (1e-6)
    std::vector<std::string> labels;           // n_k
    std::vector<int> mesh;                     // 3 integers >= 1
    bool has_shift = false;
    std::array<double, 3> shift{{0, 0, 0}};
    bool has_time_reversal = false;
    bool time_reversal = false;
    // Orbitals per k-point and spin (axis 2 of the members below).
    int64_t n_mo = 0;
    // Optional: the defined orbitals at each k-point (n_k entries <= n_mo).
    std::vector<int64_t> n_mo_per_k;
    bool has_fermi_energy = false;
    double fermi_energy_ev = 0.0;              // written as {value, "ev"}
    // Row-major float64 [n_k, n_spin, n_mo, n_ao, 2], last axis (re, im).
    std::vector<double> mo_coefficients;
    std::vector<double> energies;              // [n_k, n_spin, n_mo], hartree
    std::vector<double> occupations;           // [n_k, n_spin, n_mo], electrons
    std::string section_id = "wf_bloch";
};

struct Spectrum {
    std::vector<double> frequencies;
    std::vector<double> intensities;
    // Optional unit names (QVF spec § 4.3) of frequencies and intensities;
    // empty omits them. Declared, never converted; add_spectrum throws on a
    // string that is not a unit name.
    std::string x_units;
    std::string y_units;
};

// A volumetric grid descriptor (bohr). shape is taken from the data tensor:
// its three spatial extents (the leading three of a volume.vector or
// volume.complex tensor, the trailing three of a volume.series tensor).
struct Grid {
    std::array<double, 3> origin{{0, 0, 0}};
    std::array<std::array<double, 3>, 3> voxel_vectors{{{{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}}};
    // Optional `units` declaration of the grid payload (QVF spec § 4.3,
    // "Declared units"): a grid field or member role -> unit name, e.g.
    // {{"origin", "bohr"}, {"voxel_vectors", "bohr"}}. Declared, never
    // converted: origin and voxel_vectors can only restate "bohr", and the
    // add_* call throws on anything else or on a string that is not a unit name.
    std::map<std::string, std::string> units;
    // Optional periodic grid semantics (QVF spec § 5.2); each is written only
    // when set. periodic[a] true: the field repeats with lattice vector a of
    // the grid's structure along data axis a. includes_endpoint (absent means
    // true, cube style): the last plane along a periodic axis repeats the
    // first, so voxel_vectors[a] * (shape[a] - 1) spans lattice vector a;
    // false (CHGCAR style): voxel_vectors[a] * shape[a] spans it. Validators
    // check the span against the structure's cell; the writer does not.
    // structure_ref (empty = omit) names the structure section whose cell the
    // grid samples; the add_* call throws if it names none.
    std::optional<std::array<bool, 3>> periodic;
    std::optional<bool> includes_endpoint;
    std::string structure_ref;
};

// One program invocation's self-contained record (kind run.record): the
// verbatim input handed to a quantum-chemistry code and/or the full
// log/output it produced. `program` names the code that ran — which may
// differ from the archive-level Source when the QVF is written by a
// converter. `input_text` / `log_text` are stored as opaque UTF-8 bytes;
// empty optional strings mean "omit the field".
struct RunRecord {
    std::string program;               // REQUIRED, e.g. "orca", "vibe-qc"
    std::string input_text;            // verbatim main input (UTF-8)
    std::string log_text;              // full log/output (UTF-8)
    std::string program_version;
    std::string command;               // launch command line, one string
    bool has_exit_status = false;
    int exit_status = 0;
    std::string started_utc;           // ISO 8601 UTC
    std::string finished_utc;
    int sequence = -1;                 // >= 0 to emit (multi-run ordering)
    Json files;                        // optional role -> {filename, ...} index
    // Auxiliary files: role suffix -> raw bytes, stored as attachment.<suffix>.
    std::vector<std::pair<std::string, std::vector<uint8_t>>> attachments;
    // In a job set (QVF spec § 5.8, § 5.9): the id of the job.spec section
    // this run executed, which must already be added, and the run's outcome,
    // "running", "converged" or "failed". Empty omits a field.
    std::string job_ref;
    std::string run_status;
};

// An opaque binary member (no dtype or shape, QVF spec § 4.2) of a vendor
// section: its role, its bytes and, optionally, the IANA media type of the
// bytes ("type/subtype" without parameters, e.g. "image/svg+xml"), which an
// opaque vendor member SHOULD carry (spec § 7.3). Empty media_type omits it.
struct OpaqueMember {
    std::string role;
    std::vector<uint8_t> bytes;
    std::string media_type;
};

// kBohrPerAngstrom (Bohr per Ångström, CODATA 2018) and the other unit
// constants and conversions are in qvf/units.hpp, included above.

// A periodic cell that converts between QVF's two lattice conventions (spec
// § 4.3, "Lattice conventions"):
//   * structure.lattice_vectors — three ROW vectors in ÅNGSTRÖM; row i pairs
//     with pbc[i];
//   * reaction.path.lattice     — a 3x3 matrix whose COLUMNS are a, b, c, in
//     BOHR (one matrix per frame in the [n_frames, 3, 3] form).
// The rows in Å are stored internally; the columns in bohr are
// rows_angstrom() transposed and scaled by kBohrPerAngstrom. Construction
// validates that every entry is finite (the 3x3 shape is the type) and throws
// std::runtime_error otherwise; a degenerate or non-physical row (an aperiodic
// axis, spec § 5.1) is accepted.
struct Lattice {
    using Mat3 = std::array<std::array<double, 3>, 3>;

    // structure.lattice_vectors: row i is the lattice vector of axis i, in Å.
    static Lattice from_rows_angstrom(const Mat3& rows);
    // reaction.path.lattice (one matrix): column j is lattice vector j, in bohr.
    static Lattice from_columns_bohr(const Mat3& columns);

    Mat3 rows_angstrom() const;     // the structure.lattice_vectors form
    Mat3 columns_bohr() const;      // the reaction.path.lattice form
    double volume_angstrom3() const;   // |a . (b x c)|
    double volume_bohr3() const;

private:
    explicit Lattice(const Mat3& rows_angstrom) : rows_(rows_angstrom) {}
    Mat3 rows_;                     // rows, Å
};

// The optional per-frame members of a trajectory (QVF spec § 5.5). Every
// pointer may be nullptr, which omits the member; each tensor is float64 and
// in frame order, with the frame count (and, for the vectors, the atom count)
// of `coords`. `reaction.path` carries the four shared members (velocities,
// gradient, times, energies); its add_reaction_path overload refuses the
// trajectory-only ones (lattice, cell_volumes, temperatures).
struct TrajectoryExtras {
    // The cell in the STRUCTURE convention: rows in Å (row i is the lattice
    // vector of axis i, paired with metadata pbc[i]), [3, 3] for a fixed cell
    // or [n_frames, 3, 3] for one cell per frame. Not the columns in bohr of
    // add_reaction_path's lattice (spec § 4.3, "Lattice conventions"). The
    // metadata must carry "pbc" when it is given.
    const Tensor* lattice = nullptr;
    const Tensor* velocities = nullptr;    // [n_frames, n_atoms, 3], Å/fs
    const Tensor* gradient = nullptr;      // [n_frames, n_atoms, 3], dE/dR in hartree/bohr
    const Tensor* times = nullptr;         // [n_frames], fs
    const Tensor* energies = nullptr;      // [n_frames], hartree; must match metadata
                                           // "energies" within 1e-10 when both exist
    const Tensor* cell_volumes = nullptr;  // [n_frames], Å^3
    const Tensor* temperatures = nullptr;  // [n_frames], kelvin
};

// Normalization factor for the *axial* Cartesian Gaussian primitive of a
// shell. Depends only on l, not on (l_x, l_y, l_z), so for `pure: false`
// shells the primitives are not individually unit-normalized (<xy|xy> = 1/3
// for a d shell). Consumers must not add a per-component correction on top
// of it -- see spec Appendix A.1.
double primitive_norm(double alpha, int l);

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

class QvfWriter;
namespace detail {
struct WriterOutput;   // the archive file a streaming or updating writer is bound to
}

// A binary member being streamed (QvfWriter::begin_binary). Push its
// little-endian, C-contiguous bytes with write() in slabs of any size, then
// close(): the member spec (path, format, dtype, shape, sha256, computed on
// the fly) joins its section. The bytes must total itemsize x prod(shape).
// Destroyed without close(), or closed short, the member is abandoned and
// never listed. The writer must outlive the sink and not move while it is
// open; only one sink of a writer is open at a time.
class MemberSink {
public:
    MemberSink(MemberSink&&) noexcept;
    MemberSink& operator=(MemberSink&&) noexcept;
    ~MemberSink();

    void write(const void* data, size_t n);
    void write(const std::vector<uint8_t>& bytes) { write(bytes.data(), bytes.size()); }
    void close();
    uint64_t bytes_written() const;
    uint64_t bytes_expected() const;
    const std::string& path() const;   // the member's ZIP entry name

private:
    friend class QvfWriter;
    struct State;
    explicit MemberSink(std::unique_ptr<State> state);
    void abandon() noexcept;
    std::unique_ptr<State> s_;
};

class QvfWriter {
public:
    explicit QvfWriter(Source source);
    // A streaming writer: the archive `path` (".qvf" appended if missing) is
    // created now, every member is written to it as it is added (memory
    // holds one member at a time), and write() writes manifest.json last and
    // finishes the file. to_bytes() is not available.
    QvfWriter(Source source, const std::string& path);
    // One-shot streaming output (memory, file, or non-seekable callback).
    QvfWriter(Source source, std::unique_ptr<container::Sink> sink);
    // Open an existing archive for an in-place update (spec § 2.3): the
    // manifest is loaded (source, root blocks, sections), add_* and
    // begin_binary() append new members, remove_section() removes a section
    // (remove and add under the same id to replace one), and write() appends
    // the new members and a new manifest.json, STOREd and last. The previous
    // manifest, and the members of removed or replaced sections, stay in the
    // archive as superseded entries (W_ORPHAN_MEMBER until repack()); readers
    // use the last manifest. An update never rewrites an existing entry: a
    // new member whose default path the archive already holds (an
    // atom_properties meta.json, or a member of a replaced section) gets a
    // fresh path "<stem>.<n><suffix>" (n = 2, 3, ...; the stem ends at the
    // first '.' of the file name, so "s/structure.json" becomes
    // "s/structure.2.json"). Root keys the writer does not model are kept;
    // auto_metadata is off, and metadata.modified_utc is stamped unless
    // set_metadata() gives one. An archive stamped with the withdrawn
    // qvf_version 2 is read as version 1, and the update writes 1 and the v1
    // schema_uri (spec § 3.1). Until write() the file is untouched
    // (to_bytes() returns the updated archive in memory) unless
    // begin_binary() starts streaming into it. Throws for an archive that
    // repeats an entry name other than manifest.json, has no manifest.json,
    // a qvf_version other than 1 or 2, or a byte_order other than "little".
    static QvfWriter open_for_update(const std::string& path);

    // --- root metadata (each takes a qvf::Json object; see docs) ---
    // set_provenance drops a top-level "hostname" key: the field is deprecated
    // (QVF spec § 3.2) because archives are meant to be shared and machine
    // names are private operational detail. Pass allow_hostname = true to
    // keep it. A Json::raw fragment is parsed for the checks (throws if it is
    // not valid JSON) and stays verbatim unless the key had to be removed.
    // The block must be a JSON object, and its typed fields are checked (the
    // call throws, and changes nothing, on a wrong one): `software` (an array
    // of {name, version?, url?} objects, `name` non-empty, further keys must
    // begin x_), `git_commit` (7 to 64 hexadecimal characters), `n_threads`
    // and `n_mpi_ranks` (integers >= 1), `n_gpus` (an integer >= 0) and
    // `started_utc` / `finished_utc` (strings YYYY-MM-DDThh:mm:ssZ, UTC; QVF
    // spec § 3.2). Every other key of the open block is written as given.
    void set_provenance(Json provenance);  // allow_hostname = false
    void set_provenance(Json provenance, bool allow_hostname);
    void set_thermochemistry(Json thermochemistry);
    void set_dipole_moment(Json dipole);
    void set_constraints(Json constraints);
    // Root energies block (QVF spec § 3.6): {"total": {value, units}, "method"?,
    // "components"?, "by_method"?, "units_note"?}; write energies in hartree.
    // An object must carry "total"; a Json::raw fragment (the C ABI path) is
    // emitted verbatim. An energy stated for the method provenance.method
    // names must agree with provenance.scf_energy (validators check it).
    void set_energies(Json energies);
    void set_extensions(Json extensions);
    void set_viewer_defaults(Json viewer_defaults);
    // The root `metadata` block (QVF spec § 3.5): title, description,
    // created_utc, modified_utc, uuid, authors, license, doi, keywords,
    // derived_from, tags and x_-prefixed keys. Replaces any earlier block.
    // Throws if `metadata` is not a JSON object (a Json::raw fragment is
    // parsed first); its fields are checked by the schema, not here.
    void set_metadata(Json metadata);
    // When true (the default), build_manifest() adds `created_utc` (now, UTC)
    // and a random version-4 `uuid` to the metadata block unless
    // set_metadata() supplied them. Both are made on the first
    // build_manifest() and reused afterwards. Set it to false for a manifest
    // that is the same on every run, such as a test vector.
    bool auto_metadata = true;

    // --- section references and grouping (QVF spec § 4.1) ---
    // Set field `key` of the already added section `id` to `value`. The
    // fields `id`, `kind` and `members` cannot be changed. A `structure_ref`
    // or `run_ref` must be a string naming an existing section of kind
    // `structure` or `run.record` respectively. Throws otherwise, changing
    // nothing.
    void set_section_field(const std::string& id, const std::string& key, Json value);
    // Set the generic references of section `id`: the `structure` section it
    // describes, the `run.record` of the run that produced it, and a display
    // `group` label. An empty string leaves that field as it is. Both
    // references are checked before anything is set. A section's
    // structure_ref takes precedence over the one in a wavefunction.gto basis
    // payload, so keep the two equal. Returns `id`.
    std::string annotate_section(const std::string& id,
                                 const std::string& structure_ref = "",
                                 const std::string& run_ref = "",
                                 const std::string& group = "");

    // --- removing and replacing sections (QVF spec § 2.3) ---
    // Remove the section `id` from the manifest and return it (its section
    // object, as it would have been written). The members this writer has
    // not yet written are dropped: a buffered writer discards their bytes,
    // and a streaming writer drops the entries it already streamed from the
    // archive's directory (their bytes stay in the file as unreferenced
    // space). Entries the archive held before (open_for_update, or an
    // earlier write() of a bound writer) are never rewritten: they stay as
    // superseded entries (W_ORPHAN_MEMBER until repack()). To replace a
    // section, remove it and add the new one under the same id; in an update
    // its members get fresh paths where the old ones are taken. Nothing else
    // changes: references to the removed section from other sections
    // (structure_ref, run_ref, job_ref, a reaction.scheme participant, ...)
    // dangle until a section with that id is added again, and validators
    // report any left dangling. Throws for an unknown id or while a member
    // sink is open, changing nothing.
    Json remove_section(const std::string& id);

    // --- structure & connectivity ---
    // `pbc` says which axes are periodic and needs no particular ordering; row i
    // of `lattice_vectors` pairs with pbc[i]. A row whose axis is aperiodic is
    // bookkeeping, not a cell edge (QVF spec § 5.1). Throws if pbc marks a
    // periodic axis with no lattice to span it.
    std::string add_structure(const std::vector<Atom>& atoms,
                              const std::array<bool, 3>& pbc = {{false, false, false}},
                              const std::vector<std::array<double, 3>>& lattice_vectors = {},
                              const std::vector<Bond>& bonds = {},
                              const std::string& section_id = "structure");
    // Same, with the cell as a Lattice: written as its rows_angstrom(), so the
    // payload is identical to passing those three rows as lattice_vectors.
    std::string add_structure(const std::vector<Atom>& atoms,
                              const std::array<bool, 3>& pbc,
                              const Lattice& lattice,
                              const std::vector<Bond>& bonds = {},
                              const std::string& section_id = "structure");
    // The chemical identity of the species of the already added structure
    // section `id` (QVF spec § 5.1, "Chemical identity"), set on its section
    // object. `identity` is a JSON object (a Json::raw fragment is parsed)
    // with any of: "charge" (an integer), "multiplicity" (an integer >= 1),
    // "name" (a non-empty string), "smiles" (a non-empty string without
    // whitespace), "inchi" (beginning "InChI=", without whitespace) and
    // "inchikey" (AAAAAAAAAAAAAA-BBBBBBBBBB-C, uppercase letters). A given
    // field replaces the section's value and a null removes it; the others
    // stay. The section's charge and multiplicity after the call must be
    // possible for its atoms (when every atomic number is a non-negative
    // integer): the electron count sum(Z) - charge is not negative, and
    // multiplicity - 1 does not exceed it and has its parity. Throws on
    // another kind, an unknown key, a malformed value or an impossible
    // charge and multiplicity, changing nothing. Use set_section_field(id,
    // "label", ...) for the species' label in a scheme (its compound number).
    // In an archive with one structure section, validators also require
    // charge and multiplicity to equal provenance.charge and
    // provenance.multiplicity. Returns `id`.
    std::string set_structure_identity(const std::string& id, const Json& identity);
    // A standalone `bonds` section (QVF spec § 5.1). `structure_ref`, when
    // not empty, names the already added `structure` section whose atoms the
    // indices address and is written on the section object; without it the
    // bonds belong to the archive's first structure, so an archive with
    // several structures should pass it. A reference that names no added
    // structure section throws before anything is added; validators check
    // the indices against that structure's atoms.
    std::string add_bonds(const std::vector<Bond>& pairs,
                          const std::string& section_id = "bonds",
                          const std::string& structure_ref = "");
    std::string add_bond_orders(const std::string& method,
                                const Json& pairs,
                                const std::string& section_id = "bond_orders");

    // --- wavefunction (GBW-equivalent) ---
    // Checks everything (shapes, complex layout, spin, ECP entries,
    // n_electrons against the occupations, symmetry_labels lengths, the
    // AO-basis matrices and their pairing) and throws before writing anything.
    std::string add_wavefunction_gto(const WavefunctionGTO& wf);
    // Add the `overlap` member, [n_ao, n_ao], to the already added
    // wavefunction.gto section `section_id` (the C ABI's
    // qvf_set_wavefunction_overlap). Throws on another kind, a shape that is
    // not the basis's [n_ao, n_ao], or a section that already has one.
    std::string set_wavefunction_overlap(const std::string& section_id,
                                         const Matrix& overlap);
    // Set the basis payload's `ecp` array of the already added
    // wavefunction.gto section `section_id`, replacing any earlier one (the C
    // ABI's qvf_set_wavefunction_ecp_json). A Json::raw fragment is parsed.
    // Throws on another kind or a malformed entry, changing nothing.
    std::string set_wavefunction_ecp(const std::string& section_id, const Json& ecp);
    // k-resolved Bloch orbitals (wavefunction.bloch, spec § 5.3.1,
    // provisional). Throws before anything is written on an unknown spin, no
    // k-points, array sizes that disagree with [n_k, n_spin, n_mo, n_ao(, 2)],
    // weights that are not n_k non-negative numbers summing to 1 within 1e-6,
    // labels / n_mo_per_k of the wrong length, or a mesh that is not three
    // integers >= 1. Returns the section id.
    std::string add_wavefunction_bloch(const WavefunctionBloch& wf);

    // --- spectra ---
    // spectrum.x_units / y_units, when set, are written as the payload's
    // x_units / y_units (unit names, QVF spec § 4.3).
    std::string add_spectrum(const std::string& kind, const Spectrum& spectrum,
                             const std::string& label = "",
                             const std::string& section_id = "");
    // NMR and other object-shaped spectra: pass the full JSON payload. For a
    // one-dimensional kind, x_units, y_units and broadening.units in the
    // payload must be unit names, and a states_ref (which makes the peaks'
    // assignment.state_index values state indices, spec § 5.4) must name an
    // added excited_states section (throws otherwise).
    std::string add_spectrum(const std::string& kind, const Json& payload,
                             const std::string& label = "",
                             const std::string& section_id = "");

    // --- analysis / provenance ---
    // atom_properties (QVF spec § 5.7). Its roles are the schema's role list:
    // the 22 scalar roles mulliken_charge, loewdin_charge, hirshfeld_charge,
    // iao_charge, spin_population, cm5_charge, chelpg_charge, resp_charge,
    // esp_charge, nbo_charge, mbis_charge, bader_charge, ddec_charge, valence,
    // bonded_valence, fukui_plus, fukui_minus, fukui_zero, dual_descriptor,
    // atomic_energy_eh, atomic_volume_bohr3 and nmr_shielding_iso_ppm, each
    // float64 [n_atoms]; and the two vector roles magnetic_moment and
    // force_eh_bohr (the energy gradient, not the force), each float64
    // [n_atoms, 3]. A role named custom.<name> (<name> is [A-Za-z0-9_]+) is
    // also accepted once set_atom_properties_meta has described it. All roles
    // of one section have the same length; every form below throws if they
    // differ or if no role is given.
    //
    // Mulliken, Löwdin and spin populations; an empty vector omits that role.
    std::string add_atom_properties(const std::vector<double>& mulliken = {},
                                    const std::vector<double>& loewdin = {},
                                    const std::vector<double>& spin_population = {},
                                    const std::string& section_id = "atom_properties");
    // The same, for a braced first argument: add_atom_properties({-0.4, 0.4})
    // would otherwise be ambiguous with the map form below.
    std::string add_atom_properties(std::initializer_list<double> mulliken,
                                    const std::vector<double>& loewdin = {},
                                    const std::vector<double>& spin_population = {},
                                    const std::string& section_id = "atom_properties");
    // Any scalar role, keyed by role name, e.g.
    // add_atom_properties({{"hirshfeld_charge", q_h}, {"iao_charge", q_iao}}).
    // A custom.<name> key writes a scalar custom role and the section's `meta`
    // member. Throws on a role outside the list, a vector role (use
    // add_atom_vector_property), a custom role without a scalar entry from
    // set_atom_properties_meta, or an empty value vector.
    std::string add_atom_properties(const std::map<std::string, std::vector<double>>& roles,
                                    const std::string& section_id = "atom_properties");
    // Append one scalar role to the atom_properties section `section_id`,
    // creating the section on first use (the C ABI's qvf_add_atom_property).
    // Throws on a role outside the list, a vector role, a custom role without
    // a scalar `meta` entry, a role already in the section, an empty `values`,
    // a length that differs from the section's other roles, or a `section_id`
    // that names a section of another kind.
    std::string add_atom_property(const std::string& role,
                                  const std::vector<double>& values,
                                  const std::string& section_id = "atom_properties");
    // The same for a vector role (magnetic_moment, force_eh_bohr or a custom
    // role whose meta entry has kind "vector"): one {x, y, z} per atom, written
    // as float64 [n_atoms, 3]. Throws under the same conditions.
    std::string add_atom_vector_property(
        const std::string& role,
        const std::vector<std::array<double, 3>>& values,
        const std::string& section_id = "atom_properties");
    // Describe custom roles for the atom_properties sections written after
    // this call: an object mapping "custom.<name>" to {units, description,
    // kind} with kind "scalar" or "vector" (strict: only those keys and keys
    // beginning with x_). Entries merge into earlier calls, a repeated role
    // replaces its entry, and a custom role cannot be written without one. A
    // Json::raw fragment is parsed; the call throws on invalid JSON or a
    // malformed entry and then changes nothing. The writer emits each
    // section's `meta` member from the entries of its custom roles.
    void set_atom_properties_meta(const Json& meta);
    std::string add_scf_history(const Json& iterations,
                                const std::string& section_id = "scf_history");
    std::string add_citations(const std::string& bibtex,
                              const std::string& section_id = "citations");
    // Self-contained record of one program invocation (input + log +
    // attachments). Requires record.program and at least one of
    // record.input_text / record.log_text. record.job_ref, when set, must
    // name an added job.spec section and record.run_status be one of
    // "running", "converged" or "failed" (QVF spec § 5.8); the call throws
    // otherwise before anything is written.
    std::string add_run_record(const RunRecord& record,
                               const std::string& section_id = "run_record");
    // Declarative specification of the calculation the archive requests
    // (kind job.spec, spec § 5.9). `spec` is the JobSpecPayload JSON object;
    // it must carry job_type = "molecular" | "periodic". Its typed optional
    // fields are checked, and the call throws, and adds nothing, on a wrong
    // one: structure_ref (a section id; one that already names a section must
    // name a `structure`, a later one is resolved by the validator),
    // dispersion and relativistic (non-empty strings), and the objects
    // solvation {model, solvent?, epsilon?}, scf {convergence?,
    // max_iterations?, guess?, damping?}, geometry {convergence?, max_steps?,
    // coordinates?: cartesian | internal | redundant} and resources
    // {n_cores?, memory_mb?, walltime_s?}, which are hints, not requirements
    // (a key outside an object's list must begin x_; engine keywords go under
    // `options`). Pair with set_provenance({... run_status: "pending" ...})
    // for an archive that describes a job not yet run. A structure_ref is
    // written in the payload and on the section object, where a reader finds
    // it in the manifest (spec § 5.9). Several job.spec sections form a job
    // set: give each its own section_id and, when the archive holds more than
    // one structure, its structure_ref; the run.record of each job names it
    // with RunRecord::job_ref.
    std::string add_job_spec(const Json& spec,
                             const std::string& section_id = "job_spec");

    // --- reaction schemes (spec § 5.14) ---
    // A reaction.scheme section relating added structure sections as the
    // steps of a scheme. `steps` is a JSON array (a Json::raw fragment is
    // parsed) of step objects, in reading order, each with "id" (unique among
    // the steps, [A-Za-z0-9_.-]+), "reactants" and "products" (at least one
    // participant each) and optionally "label", "reagents", "arrow"
    // ("forward", the default, "equilibrium", "retrosynthetic" or
    // "resonance"), "conditions" (an array of strings, the text lines at the
    // arrow), "yield_percent" (0 to 100), "yield_text" (the yield as
    // printed), "energies" (an object, $defs/ReactionStepEnergies, written
    // as given and checked by the validator) and keys beginning x_. A participant is the id of an added
    // structure section, or an object {"structure_ref": id, "coefficient"?:
    // number > 0} (further keys must begin x_); within one role of a step a
    // structure appears once. Roles are chemical: reactants are consumed and
    // products formed in the forward direction whatever the arrow, so the
    // target of a retrosynthetic step is among its products. Null values are
    // treated as absent. `label` (empty = omit) titles the scheme. Throws
    // before anything is added on a malformed step or a participant that
    // names no structure section. Returns the section id.
    std::string add_reaction_scheme(const Json& steps, const std::string& label = "",
                                    const std::string& section_id = "reaction_scheme");

    // --- preview (spec § 5.15) ---
    // A preview section, a picture of the archive for file managers and
    // quick-look panels: `thumbnail` is PNG or JPEG bytes (its longest side
    // should be 256 to 1024 pixels), `image` an optional larger or vector
    // rendering, SVG, PNG or JPEG (empty omits it, and image_media_type is
    // then ignored). The media types are sniffed from the bytes unless given
    // (empty = sniff), and given ones are checked against the bytes: a PNG
    // begins with the PNG signature, a JPEG with FF D8 FF, and an SVG is UTF-8
    // with an <svg element. An SVG image must be self-contained and inert (no
    // external resources, no scripts); that is not checked. Both members are
    // stored uncompressed at "<id>/thumbnail.png|jpg" and
    // "<id>/image.svg|png|jpg". The picture of the archive is the thumbnail
    // of the first preview section. Throws before anything is added on bytes
    // that are not an allowed image type. Returns the section id.
    std::string add_preview(const std::vector<uint8_t>& thumbnail,
                            const std::vector<uint8_t>& image = {},
                            const std::string& label = "",
                            const std::string& section_id = "preview",
                            const std::string& thumbnail_media_type = "",
                            const std::string& image_media_type = "");

    // --- vendor escape hatch (x_<vendor>.*) ---
    // json_members: object mapping role -> JSON payload.
    std::string add_vendor_section(const std::string& kind,
                                   const Json& json_members,
                                   bool critical = false,
                                   const std::string& schema_uri = "",
                                   const std::string& section_id = "");
    // The same with opaque binary members (spec § 4.2, § 7.3), each stored at
    // "<id>/<role>.bin" with its media_type when given; `json_members` may
    // then be empty (Json::object()). Throws on a role given twice, a media
    // type that is not "type/subtype", or no member at all.
    std::string add_vendor_section(const std::string& kind,
                                   const Json& json_members,
                                   const std::vector<OpaqueMember>& opaque_members,
                                   bool critical = false,
                                   const std::string& schema_uri = "",
                                   const std::string& section_id = "");
    // Add the opaque binary member `role` (bytes without dtype or shape,
    // spec § 4.2), with its IANA media type when `media_type` is not empty,
    // to the already added vendor section `section_id`, at
    // "<section_id>/<role>.bin". Throws on a section that is not an
    // x_<vendor> section, a role it already has, or a media type that is not
    // "type/subtype" (RFC 6838, without parameters), changing nothing.
    // Returns `section_id`.
    std::string add_opaque_member(const std::string& section_id, const std::string& role,
                                  const std::vector<uint8_t>& bytes,
                                  const std::string& media_type = "");

    // --- volumes & basis AO (grid shape taken from the data Tensor) ---
    // kind is a volume kind (QVF spec § 5.2). `data` is float32 or float64:
    // [nx, ny, nz] for the scalar kinds (volume.density, orbital, spin, elf,
    // generic, potential, rdg, laplacian, transition_density; for the last,
    // add_volume_transition_density also writes its section fields);
    // [nx, ny, nz, 3] for volume.vector
    // (Cartesian x, y, z on the last axis); [nx, ny, nz, 2] for volume.complex
    // (last axis [re, im]). volume.series goes through add_volume_series.
    // Throws on another rank, trailing axis or dtype, or a grid.structure_ref
    // naming no structure section, before anything is written. A non-empty
    // grid.units is written as the grid payload's `units`; grid.periodic,
    // grid.includes_endpoint and grid.structure_ref when set.
    std::string add_volume(const std::string& kind, const Grid& grid,
                           const Tensor& data, const std::string& label = "",
                           const std::string& component = "",
                           const std::string& section_id = "");
    // A volume.series section (QVF spec § 5.2): `data` is float32 or float64
    // [n, nx, ny, nz], slab s the field at axis value s; `axis` is the
    // VolumeSeriesAxisPayload object written as the `axis` member: "kind"
    // ("time", "frame", "state" or "custom"), "values" (n numbers; 0-based
    // frame indices for a frame axis), and optionally "units" (a unit name),
    // "labels" (n strings), "trajectory_ref" (an added trajectory or
    // reaction.path section) and x_-prefixed keys. A Json::raw fragment is
    // parsed. Throws on a malformed axis or data before anything is written.
    std::string add_volume_series(const Grid& grid, const Tensor& data,
                                  const Json& axis, const std::string& label = "",
                                  const std::string& section_id = "");
    std::string add_volume_difference(const Grid& grid, const Tensor& data,
                                      const std::string& operand_a,
                                      const std::string& operand_b,
                                      const std::string& description = "",
                                      const std::string& label = "",
                                      const std::string& section_id = "diff");
    // volume.transition_density (QVF spec § 5.2): the transition density of
    // one excitation, e/bohr^3. excited_states_ref (empty = omit) names the
    // excited_states section of the state and state_index (0 = omit; needs
    // excited_states_ref) its index there. Neither is resolved here: the
    // excited_states section usually names this one through density_ref and
    // is added after it; validators resolve both.
    std::string add_volume_transition_density(const Grid& grid, const Tensor& data,
                                              const std::string& excited_states_ref = "",
                                              int state_index = 0,
                                              const std::string& label = "",
                                              const std::string& section_id = "");
    // Name the wavefunction.gto section `wavefunction_ref` that the added
    // volume.density, volume.spin or volume.orbital section `volume_id` was
    // computed from (QVF spec § 5.2, Volumes of a wavefunction). For a
    // volume.orbital, mo_index (-1 = omit) is the 0-based row of the orbital
    // in that section's coefficient member and spin ("alpha" or "beta",
    // empty = omit) the member of an unrestricted wavefunction: required
    // with mo_coefficients_alpha/_beta, refused with mo_coefficients. The
    // wavefunction must already be added. Throws, changing nothing, on a
    // wrong kind, an unresolved reference or a row outside the member.
    // Returns `volume_id`.
    std::string link_volume_wavefunction(const std::string& volume_id,
                                         const std::string& wavefunction_ref,
                                         int64_t mo_index = -1,
                                         const std::string& spin = "");
    std::string add_basis_ao(const Json& ao_metadata, const Grid& grid,
                             const Tensor& data, const std::string& label = "",
                             const std::string& section_id = "");

    // --- bands & DOS ---
    // eigenvalues: float64 [n_spin, n_k, n_bands], eV. projections (optional,
    // fat bands): float64 [n_spin, n_k, n_bands, n_channels], or the legacy
    // [n_k, n_bands, n_channels], which means n_spin = 1. occupations
    // (optional): float64 [n_spin, n_k, n_bands], electrons per state. Shapes
    // must agree with eigenvalues (spec § 5.6); throws otherwise.
    std::string add_bands(const Json& kpath, const Tensor& eigenvalues,
                          const Tensor* projections = nullptr,
                          const std::string& section_id = "bands");
    std::string add_bands(const Json& kpath, const Tensor& eigenvalues,
                          const Tensor* projections, const Tensor* occupations,
                          const std::string& section_id = "bands");
    // section_meta merges optional fields onto the section object (QVF spec
    // § 5.6): n_spin (1 or 2), smearing (eV), smearing_type, fermi_energy_ev,
    // n_electrons and, for dos.projected, channels (one label object per
    // projection row: atom, l, m, spin, label).
    std::string add_dos_total(const Tensor& energies, const Tensor& dos,
                              const Json& section_meta = Json::object(),
                              const std::string& section_id = "dos_total");
    // projections: [n_spin, n_channels, n_points], or the legacy
    // [n_channels, n_points], which means n_spin = 1 (spec § 5.6). A rank-3
    // tensor's spin axis must equal an integer n_spin in section_meta.
    std::string add_dos_projected(const Tensor& energies, const Tensor& projections,
                                  const Json& section_meta = Json::object(),
                                  const std::string& section_id = "dos_projected");
    std::string add_dos_coop_cohp(const std::string& kind, const Tensor& energies,
                                  const Tensor& projections, const Tensor& integrated,
                                  const Json& meta,
                                  const std::string& section_id = "");

    // --- trajectories / reactions / scans / vibrations ---
    // A `units` object in the metadata of add_trajectory, add_reaction_path
    // and add_scan_surface (and in the payload of add_topology_qtaim) is a
    // declaration (QVF spec § 4.3, "Declared units"): each value must be a
    // unit name, and the units the spec fixes can only be restated:
    // trajectory `energies` "hartree", `coords` "angstrom", `times` and
    // `timestep_fs` "fs", `lattice` "angstrom", `cell_volumes` "a3",
    // `temperatures` "kelvin"; a reaction path's `lattice` "bohr". These calls
    // throw otherwise; nothing is converted.
    //
    // add_trajectory (QVF spec § 5.5): `coords` is [n_frames, n_atoms, 3] in
    // Å, float64, or float32 for a long MD run: float32 coords write
    // "precision": "float32" into the metadata (added when absent), and the
    // call throws if the metadata names a precision other than the dtype of
    // `coords`. `extras` adds the optional per-frame members; with a lattice
    // the metadata must carry "pbc" (three booleans), and a true pbc entry
    // requires the lattice. Everything is checked before anything is written.
    std::string add_trajectory(const Json& metadata, const Tensor& coords,
                               const std::string& section_id = "trajectory");
    std::string add_trajectory(const Json& metadata, const Tensor& coords,
                               const TrajectoryExtras& extras,
                               const std::string& section_id = "trajectory");
    // A fixed cell from a Lattice, written as its rows_angstrom() [3, 3]
    // (extras.lattice must then be nullptr).
    std::string add_trajectory(const Json& metadata, const Tensor& coords,
                               const Lattice& lattice,
                               const TrajectoryExtras& extras = TrajectoryExtras(),
                               const std::string& section_id = "trajectory");
    // One cell per frame (lattices.size() == n_frames), written as the
    // [n_frames, 3, 3] rows in Å (extras.lattice must then be nullptr).
    std::string add_trajectory(const Json& metadata, const Tensor& coords,
                               const std::vector<Lattice>& lattices,
                               const TrajectoryExtras& extras = TrajectoryExtras(),
                               const std::string& section_id = "trajectory");
    // `lattice` (optional, may be nullptr) marks the path periodic: float64
    // [3,3] (fixed cell) or [n_frames,3,3] (variable cell), columns = a,b,c, in
    // bohr. Presence — not the manifest version — signals a periodic path. This
    // is NOT the layout or unit of add_structure's lattice_vectors (rows, Å);
    // see spec § 4.3, "Lattice conventions". The Tensor is written as given.
    std::string add_reaction_path(const Json& metadata, const Tensor& coords,
                                  const Tensor* lattice = nullptr,
                                  const std::string& section_id = "reaction_path");
    // Fixed-cell periodic path from a Lattice: written as the [3,3] matrix
    // columns_bohr(), so the conversion to the on-disk convention is done for
    // you. nullptr writes a molecular path, as for the Tensor overload.
    std::string add_reaction_path(const Json& metadata, const Tensor& coords,
                                  const Lattice* lattice,
                                  const std::string& section_id = "reaction_path");
    // A literal `nullptr` argument would be ambiguous between the two pointer
    // overloads above; this keeps existing calls that pass one compiling.
    std::string add_reaction_path(const Json& metadata, const Tensor& coords,
                                  std::nullptr_t,
                                  const std::string& section_id = "reaction_path");
    // The same with the per-frame members a reaction path shares with a
    // trajectory: extras.velocities, gradient, times and energies (shapes and
    // units as for add_trajectory). extras.lattice, cell_volumes and
    // temperatures must be nullptr: a reaction path's cell is `lattice`
    // above, in columns and bohr. Coordinates stay float64.
    std::string add_reaction_path(const Json& metadata, const Tensor& coords,
                                  const Tensor* lattice, const TrajectoryExtras& extras,
                                  const std::string& section_id = "reaction_path");
    std::string add_reaction_waypoints(const std::string& trajectory_ref,
                                       const Json& waypoints,
                                       const std::string& section_id = "reaction_waypoints");
    std::string add_scan_surface(const Json& metadata, const Tensor& axis_a,
                                 const Tensor& axis_b, const Tensor& energies,
                                 const Tensor* geometries = nullptr,
                                 const std::string& section_id = "scan_surface");
    // `hessian_ref` (optional, empty = omit) names the hessian section the
    // modes were derived from (spec § 5.5); it must already be added.
    std::string add_vibrations(const Json& metadata, const Tensor& displacements,
                               const std::string& section_id = "vibrations",
                               const std::string& hessian_ref = "");

    // --- energy derivatives (spec § 5.10, § 5.11) ---
    // `gradient`: float64 [n_atoms, 3], the energy GRADIENT dE/dR in
    // hartree/bohr (the forces are its negative), rows in the atom order of
    // the structure. `stress` (optional, may be nullptr, periodic systems):
    // float64 [3, 3] stress tensor (1/V) dE/d(strain) in hartree/bohr^3.
    // `meta` (optional; Json() omits it): the ForcesMetaPayload object
    // (structure_ref, method, convention = "gradient", units). Throws on a
    // wrong dtype or shape before anything is written.
    std::string add_forces(const Tensor& gradient, const Tensor* stress = nullptr,
                           const Json& meta = Json(),
                           const std::string& section_id = "forces");
    // `matrix`: float64 [3N, 3N] Cartesian Hessian d2E/(dx_i dx_j) in
    // hartree/bohr^2, index 3 * atom + xyz; symmetric. `meta` (optional;
    // Json() omits it): the HessianMetaPayload object (structure_ref,
    // mass_weighted, projected, atom_masses_amu, units). Throws on a
    // non-square matrix, an extent not divisible by 3, a mass-weighted matrix
    // without atom_masses_amu, or masses that are not one per atom, before
    // anything is written.
    std::string add_hessian(const Tensor& matrix, const Json& meta = Json(),
                            const std::string& section_id = "hessian");

    // --- excited states (spec § 5.12) ---
    // `states`: a JSON array of state objects in ascending order of `index`
    // (an integer >= 1, unique), with the fields of $defs/ExcitedState:
    // index and excitation_energy (required), total_energy, multiplicity,
    // symmetry, oscillator_strength, rotatory_strength_cgs (10^-40 esu^2
    // cm^2), transition_dipole and transition_magnetic_dipole (3 numbers,
    // atomic units), contributions ({from_mo, to_mo, spin?, weight}), nto_ref,
    // density_ref and x_ keys. excitation_energy / total_energy are {value,
    // units} objects with a unit name, or numbers, written in "ev" /
    // "hartree". nto_ref must name an added wavefunction.gto section,
    // density_ref an added volume.transition_density or volume.difference.
    // `meta` (optional; Json() omits it): the ExcitedStatesMeta object
    // (method, reference_state, n_states, tda, spin_flip, units); n_states
    // must equal the number of states. `wavefunction_ref` (empty = omit)
    // names the added wavefunction.gto of the reference state; every
    // contribution's from_mo / to_mo must then be one of its orbital rows.
    // A Json::raw fragment is parsed. Throws on anything else before anything
    // is written.
    std::string add_excited_states(const Json& states, const Json& meta = Json(),
                                   const std::string& wavefunction_ref = "",
                                   const std::string& section_id = "excited_states");

    // --- property kinds (spec § 5.13) ---
    // Each takes the `properties` payload as a JSON object (a Json::raw fragment
    // is parsed), checks it, and writes the section's single `properties`
    // member as given. Every call throws before anything is written on a
    // payload the schema would reject that can be judged without a schema
    // engine: an unknown key (keys beginning x_ are allowed), a block or matrix
    // of the wrong shape or with a non-finite number, a missing required
    // field, a value outside an enumeration, or a `units` declaration that is
    // not a unit name or contradicts the fixed unit. Units are fixed and
    // nothing is converted.
    //
    // properties.electric: at least one of `dipole` ({"vector": [x, y, z],
    // "total"?, "units"?: "debye"}, debye), `quadrupole` ({"matrix": 3x3
    // symmetric, "traceless": bool, "convention": "buckingham" | "cartesian"},
    // debye-angstrom), `octupole` ({"tensor": 27 numbers, row-major [3, 3, 3],
    // "convention"?}, debye-angstrom^2), `polarizability` ({"matrix": 3x3,
    // "frequency_ev"?, "isotropic"?}, atomic units) and `hyperpolarizability`
    // ({"tensor": 27 numbers, "frequencies_ev"?: [w1, w2]}, atomic units);
    // optional `origin` (a name or a 3-vector in bohr) and `units` ({"dipole":
    // "debye", "origin": "bohr"}). A dipole that the root dipole_moment also
    // states must agree with it within 1e-6 debye (validators check).
    std::string add_electric_properties(const Json& properties,
                                        const std::string& section_id = "electric");
    // properties.elastic: `stiffness` (6x6, Voigt order xx yy zz yz xz xy, GPa,
    // symmetric, required); optional `compliance` (6x6, 1/GPa, symmetric),
    // `bulk_modulus_gpa`, `shear_modulus_gpa`, `youngs_modulus_gpa`,
    // `poisson_ratio`, `average` ("voigt" | "reuss" | "hill") and `units`.
    std::string add_elastic_properties(const Json& properties,
                                       const std::string& section_id = "elastic");
    // properties.dielectric: `electronic` (3x3, required); optional `ionic`,
    // `static` (3x3; static must equal electronic + ionic within 1e-6 when all
    // three are given), `born_charges` (array of 3x3, one per atom, in e),
    // `piezoelectric` (3x6, C/m^2) and `units`. That there is one Born tensor
    // per atom of the section's structure is checked by validators.
    std::string add_dielectric_properties(const Json& properties,
                                          const std::string& section_id = "dielectric");

    // --- symmetry & periodic kinds ---
    std::string add_structure_symmetry(const Json& data,
                                       const std::string& section_id = "symmetry");
    std::string add_fermi_surface(const Json& mesh, const Tensor& energies,
                                  const std::string& section_id = "fermi_surface");
    // `eigenvectors` (optional, may be nullptr): float64 [n_q, n_modes, n_atoms, 3]
    // when real, or [n_q, n_modes, n_atoms, 3, 2] with the last axis [re, im]
    // when complex; its [n_q, n_modes] must match `frequencies` (QVF spec
    // § 5.6). Throws on another rank, trailing axis or dtype, or a mismatch,
    // before anything is written.
    std::string add_phonon_bands(const Json& qpath, const Tensor& frequencies,
                                 const Tensor* eigenvectors = nullptr,
                                 const std::string& section_id = "phonon_bands");
    // `projected` (optional, may be nullptr): float64 [n_channels, n_points] on
    // the [n_points] `frequencies` grid. Label its rows with a `channels` array
    // in `meta`, one object per row; when `meta` has one, its length must be
    // n_channels. Throws otherwise, before anything is written.
    std::string add_phonon_dos(const Json& meta, const Tensor& frequencies,
                               const Tensor& dos, const Tensor* projected = nullptr,
                               const std::string& section_id = "phonon_dos");
    std::string add_equation_of_state(const Tensor& volumes, const Tensor& energies,
                                      const Json& fit, const std::string& section_id = "eos");
    std::string add_topology_qtaim(const Json& critical_points,
                                   const std::string& section_id = "qtaim");

    // Low-level typed-facade entry points. Declare before adding members;
    // member roles must be new and cannot introduce paths outside the section.
    std::string add_section(Section section);
    void add_json_member(const std::string& id, const std::string& role, Json value);
    void add_section_opaque(const std::string& id, const std::string& role,
                           const std::vector<uint8_t>& bytes,
                           const std::string& media_type = "");

    // --- streaming members (WP-40) ---
    // Stream a binary member `role` into the already added section
    // `section_id` (any kind; the role must be new). Its path is
    // "<section_id>/<role>.bin", or a fresh one when an archive being
    // updated already holds that (MemberSink::path() names it). A streaming
    // or updating writer writes the bytes straight into the archive; a
    // buffered writer keeps them until to_bytes()/write(). Close the sink
    // before adding anything else.
    MemberSink begin_binary(const std::string& section_id, const std::string& role,
                            DType dtype, std::vector<int64_t> shape);
    // Write ZIP64 structures (extra fields, end record and locator) even
    // though nothing needs them. Off by default; ZIP64 is used anyway when a
    // size, an offset or the entry count passes the 32-bit limits.
    bool force_zip64 = false;

    // --- output ---
    Json build_manifest() const;
    // The archive: members in the order added, then manifest.json, STOREd.
    std::vector<uint8_t> to_bytes();
    // Write the archive to `path` (".qvf" appended if missing), streaming
    // it to the file. A writer bound to a file accepts only its own path.
    std::string write(const std::string& path);
    // Finish a writer bound to a file (QvfWriter(source, path) or
    // open_for_update): manifest.json last, then the central directory. The
    // writer stays usable: later members are appended by the next write().
    std::string write();

private:
    friend class MemberSink;
    void open_bound_archive();                       // appends to out_->path
    void emit_archive(container::ZipWriter& zip);
    Json final_manifest() const;
    std::shared_ptr<detail::WriterOutput> out_;      // null for a buffered writer
    bool member_open_ = false;

    std::string reserve_id(const std::string& id);
    void check_new_id(const std::string& id) const;   // reserve_id's checks only
    // The entry name a new member with the default `path` gets: `path`, or,
    // in a bound writer, a fresh "<stem>.<n><suffix>" when the archive
    // already holds `path` or its manifest names it. Throws for a path a
    // member of this session already uses.
    std::string entry_name(const std::string& path) const;
    bool pending_file(const std::string& path) const;   // written this session
    // Adds the member file under entry_name(path) and returns that name;
    // `store` writes it uncompressed (an already compressed image).
    std::string add_file(const std::string& path, std::vector<uint8_t> bytes,
                         bool store = false);
    void drop_file(const std::string& path);
    Json json_member(const std::string& path, const Json& obj);
    Json binary_member(const std::string& path, const Tensor& t, bool with_shape = true);
    Json opaque_member(const std::string& path, std::vector<uint8_t> bytes,
                       const std::string& media_type = "", bool store = false);
    // Throws unless `ref` names an added section of kind `kind`.
    void require_section_kind(const std::string& what, const std::string& ref,
                              const char* kind) const;
    // atom_properties: one role as a [n_atoms] or [n_atoms, 3] float64 member.
    std::string add_atom_role(const std::string& role, const std::vector<double>& flat,
                              size_t n_atoms, bool vector, const std::string& section_id);
    Json atom_properties_meta_payload(const std::vector<std::string>& custom_roles) const;
    // trajectory / reaction.path: checks every member, then writes the section.
    std::string add_frames_section(const std::string& kind, const Json& metadata,
                                   const Tensor& coords, const Tensor* lattice,
                                   const TrajectoryExtras& extras,
                                   const std::string& section_id);
    // properties.*: the section with its single `properties` member.
    std::string add_properties_section(const std::string& kind,
                                       const std::string& section_id,
                                       const Json& payload);
    size_t section_index(const std::string& id) const;   // throws if absent
    void check_section_ref(const std::string& key, const Json& value) const;
    // A wavefunction.gto section's index; throws for an absent id or another kind.
    size_t wavefunction_index(const std::string& id) const;
    Json member_json(const Json& member) const;   // a written JSON member, parsed
    void replace_section(size_t index, Json section);

    Source source_;
    Json root_ = Json::object();          // optional root blocks
    Json sections_ = Json::array();
    std::vector<std::pair<std::string, std::vector<uint8_t>>> files_;
    std::set<std::string> ids_;
    std::set<std::string> paths_;
    std::set<std::string> stored_;        // member paths written uncompressed
    // Structure section id -> the sum of its atomic numbers, when every one is
    // a non-negative integer (for set_structure_identity, spec § 5.1).
    std::map<std::string, int64_t> structure_sum_z_;
    Json atom_properties_meta_ = Json::object();  // custom role -> description
    // The automatic metadata values, made on the first build_manifest().
    mutable std::string auto_created_utc_;
    mutable std::string auto_uuid_;
};

// Rewrite the archive `in` as `out` from its current (last) manifest (spec
// § 2.3): the entries a member path references, in archive order, then that
// manifest.json, STOREd, last. Orphans and superseded manifests are dropped;
// member bytes and the manifest are copied unchanged, so every sha256 and
// size stays valid. `in` and `out` may name the same file.
void repack(const std::string& in, const std::string& out);

}  // namespace qvf

#endif  // QVF_QVF_HPP
