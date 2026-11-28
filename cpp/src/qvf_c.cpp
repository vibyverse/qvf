// qvf_c.cpp — C ABI implementation over qvf::QvfWriter. License: Apache-2.0.
#include "qvf/qvf_c.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "qvf/qvf.hpp"

#include "qvf_c_internal.hpp"

// --- volumes: grid JSON and raw data (spec § 5.2) ----------------------------

static std::array<double, 3> vec3_of(const qvf::Json& j, const char* what) {
    if (!j.is_array() || j.size() != 3)
        throw std::runtime_error(std::string("grid_json: ") + what +
                                 " must be three numbers");
    std::array<double, 3> v{};
    for (size_t i = 0; i < 3; ++i) {
        if (!j.at(i).is_number())
            throw std::runtime_error(std::string("grid_json: ") + what +
                                     " must be three numbers");
        v[i] = j.at(i).as_double();
    }
    return v;
}

// A Grid from the C ABI's grid JSON; `spatial` is the [nx, ny, nz] of the data.
static qvf::Grid grid_from_json(const char* grid_json, const std::vector<int64_t>& spatial) {
    if (!grid_json) throw std::runtime_error("grid_json is NULL");
    const qvf::Json doc = qvf::Json::parse(grid_json);
    if (!doc.is_object()) throw std::runtime_error("grid_json must be a JSON object");
    if (!doc.has("origin") || !doc.has("voxel_vectors"))
        throw std::runtime_error("grid_json needs origin and voxel_vectors");
    qvf::Grid g;
    for (const auto& kv : doc.items()) {
        const std::string& key = kv.first;
        const qvf::Json& v = kv.second;
        if (key == "origin") {
            g.origin = vec3_of(v, "origin");
        } else if (key == "voxel_vectors") {
            if (!v.is_array() || v.size() != 3)
                throw std::runtime_error("grid_json: voxel_vectors must be 3 x 3 numbers");
            for (size_t i = 0; i < 3; ++i) g.voxel_vectors[i] = vec3_of(v.at(i), "voxel_vectors");
        } else if (key == "shape") {
            bool same = v.is_array() && v.size() == spatial.size();
            for (size_t i = 0; same && i < spatial.size(); ++i)
                same = v.at(i).type() == qvf::Json::Type::Int && v.at(i).as_int() == spatial[i];
            if (!same)
                throw std::runtime_error("grid_json: shape differs from the data's "
                                         "spatial extents");
        } else if (key == "units") {
            if (!v.is_object())
                throw std::runtime_error("grid_json: units must be an object");
            for (const auto& u : v.items()) {
                if (!u.second.is_string())
                    throw std::runtime_error("grid_json: units values must be strings");
                g.units[u.first] = u.second.as_string();
            }
        } else if (key == "periodic") {
            if (!v.is_array() || v.size() != 3)
                throw std::runtime_error("grid_json: periodic must be three booleans");
            std::array<bool, 3> flags{};
            for (size_t i = 0; i < 3; ++i) {
                if (!v.at(i).is_bool())
                    throw std::runtime_error("grid_json: periodic must be three booleans");
                flags[i] = v.at(i).as_bool();
            }
            g.periodic = flags;
        } else if (key == "includes_endpoint") {
            if (!v.is_bool())
                throw std::runtime_error("grid_json: includes_endpoint must be a boolean");
            g.includes_endpoint = v.as_bool();
        } else if (key == "structure_ref") {
            if (!v.is_string() || v.as_string().empty())
                throw std::runtime_error("grid_json: structure_ref must be a section id");
            g.structure_ref = v.as_string();
        } else {
            throw std::runtime_error("grid_json: unknown key '" + key + "'");
        }
    }
    return g;
}

// Largest element count the volume entry points accept (2^32 values).
static const int64_t kMaxVolumeValues = int64_t(1) << 32;

// A float32 / float64 Tensor of `shape` from host-order C values at `data`.
static qvf::Tensor volume_tensor(const char* dtype, int rank, const int64_t* shape,
                                 const void* data) {
    if (rank < 1 || rank > 8 || !shape)
        throw std::runtime_error("volume: rank must be 1..8 with a non-NULL shape");
    std::vector<int64_t> dims(shape, shape + rank);
    int64_t count = 1;
    for (int64_t d : dims) {
        if (d < 0 || (d > 0 && count > kMaxVolumeValues / d))
            throw std::runtime_error("volume: shape extents must be >= 0 and hold at "
                                     "most 2^32 values");
        count *= d;
    }
    if (!data && count > 0) throw std::runtime_error("volume: data is NULL");
    const std::string type = dtype && *dtype ? dtype : "float32";
    const size_t n = static_cast<size_t>(count);
    if (type == "float32") {
        const float* p = static_cast<const float*>(data);
        return qvf::tensor_f32(std::vector<double>(p, p + n), dims);
    }
    if (type == "float64") {
        const double* p = static_cast<const double*>(data);
        return qvf::tensor_f64(std::vector<double>(p, p + n), dims);
    }
    throw std::runtime_error("volume: dtype must be \"float32\" or \"float64\"");
}

extern "C" {

qvf_writer* qvf_create(const char* program, const char* version,
                       const char* calculation) {
    qvf_c_detail::global_error.clear();
    try {
        qvf::Source s{program ? program : "", version ? version : "",
                      calculation ? calculation : ""};
        return new qvf_writer(std::move(s));
    } catch (const std::exception& e) {
        qvf_c_detail::global_error = e.what(); return nullptr;
    } catch (...) {
        qvf_c_detail::global_error = "unknown error"; return nullptr;
    }
}

void qvf_destroy(qvf_writer* w) { delete w; }

const char* qvf_last_error(qvf_writer* w) {
    return w ? w->last_error.c_str() : "";
}

int qvf_add_structure(qvf_writer* w, const char* const* symbols,
                      const int* atomic_numbers, const double* positions,
                      int n_atoms) {
    return guard(w, [&] {
        std::vector<qvf::Atom> atoms;
        atoms.reserve(static_cast<size_t>(n_atoms));
        for (int i = 0; i < n_atoms; ++i) {
            qvf::Atom a;
            if (symbols && symbols[i]) a.symbol = symbols[i];
            if (atomic_numbers) a.atomic_number = atomic_numbers[i];
            if (positions) {
                a.position = {{positions[i * 3], positions[i * 3 + 1],
                               positions[i * 3 + 2]}};
            }
            atoms.push_back(std::move(a));
        }
        w->backend().add_structure(atoms);
    });
}

}  // extern "C"

// The shells of qvf_add_wavefunction_gto and qvf_add_wavefunction_bloch from
// their parallel arrays (the primitive runs concatenated in shell order). It
// returns a C++ type, so it stays outside the C-linkage block (MSVC C4190).
static std::vector<qvf::Shell> shells_from_arrays(int n_shells, const int* shell_center,
                                                  const int* shell_l, const int* shell_nprim,
                                                  const double* exponents_flat,
                                                  const double* coeffs_flat, int pure) {
    std::vector<qvf::Shell> shells;
    size_t off = 0;
    for (int s = 0; s < n_shells; ++s) {
        qvf::Shell sh;
        sh.center = shell_center[s];
        sh.l = shell_l[s];
        sh.pure = pure != 0;
        int np = shell_nprim[s];
        for (int p = 0; p < np; ++p) {
            const size_t k = off + static_cast<size_t>(p);
            sh.exponents.push_back(exponents_flat[k]);
            sh.coefficients.push_back(coeffs_flat[k]);
        }
        off += static_cast<size_t>(np);
        shells.push_back(std::move(sh));
    }
    return shells;
}

extern "C" {

int qvf_add_wavefunction_gto(
    qvf_writer* w, int n_shells, const int* shell_center, const int* shell_l,
    const int* shell_nprim, const double* exponents_flat,
    const double* coeffs_flat, int pure, int coeffs_are_libint_normalized,
    int n_mo, int n_ao, const double* mo_coeffs, const double* energies,
    const double* occupations) {
    return guard(w, [&] {
        qvf::WavefunctionGTO wf;
        wf.pure = pure != 0;
        wf.coeffs_are_libint_normalized = coeffs_are_libint_normalized != 0;
        wf.shells = shells_from_arrays(n_shells, shell_center, shell_l, shell_nprim,
                                       exponents_flat, coeffs_flat, pure);
        wf.mo_coefficients.rows = n_mo;
        wf.mo_coefficients.cols = n_ao;
        wf.mo_coefficients.data.assign(
            mo_coeffs,
            mo_coeffs + static_cast<size_t>(n_mo) * static_cast<size_t>(n_ao));
        if (energies)
            wf.energies.assign(energies, energies + n_mo);
        if (occupations)
            wf.occupations.assign(occupations, occupations + n_mo);
        w->backend().add_wavefunction_gto(wf);
    });
}

int qvf_add_wavefunction_gto_unrestricted(
    qvf_writer* w, int n_shells, const int* shell_center, const int* shell_l,
    const int* shell_nprim, const double* exponents_flat,
    const double* coeffs_flat, int pure, int coeffs_are_libint_normalized,
    int n_ao, int n_mo_alpha, const double* mo_coeffs_alpha,
    const double* energies_alpha, const double* occupations_alpha,
    int n_mo_beta, const double* mo_coeffs_beta, const double* energies_beta,
    const double* occupations_beta) {
    return guard(w, [&] {
        if (n_shells < 0 || n_ao < 0 || n_mo_alpha < 0 || n_mo_beta < 0)
            throw std::runtime_error(
                "qvf_add_wavefunction_gto_unrestricted: counts must be >= 0");
        if ((n_shells && (!shell_center || !shell_l || !shell_nprim)) ||
            !mo_coeffs_alpha || !mo_coeffs_beta)
            throw std::runtime_error(
                "qvf_add_wavefunction_gto_unrestricted: shell arrays and both "
                "coefficient matrices must be non-NULL");
        qvf::WavefunctionGTO wf;
        wf.pure = pure != 0;
        wf.coeffs_are_libint_normalized = coeffs_are_libint_normalized != 0;
        wf.shells = shells_from_arrays(n_shells, shell_center, shell_l, shell_nprim,
                                       exponents_flat, coeffs_flat, pure);
        wf.unrestricted = true;
        auto fill = [n_ao](qvf::Matrix& m, int n_mo, const double* data) {
            m.rows = n_mo;
            m.cols = n_ao;
            m.data.assign(data, data + static_cast<size_t>(n_mo) * static_cast<size_t>(n_ao));
        };
        fill(wf.mo_coefficients_alpha, n_mo_alpha, mo_coeffs_alpha);
        fill(wf.mo_coefficients_beta, n_mo_beta, mo_coeffs_beta);
        if (energies_alpha) wf.energies_alpha.assign(energies_alpha, energies_alpha + n_mo_alpha);
        if (occupations_alpha)
            wf.occupations_alpha.assign(occupations_alpha, occupations_alpha + n_mo_alpha);
        if (energies_beta) wf.energies_beta.assign(energies_beta, energies_beta + n_mo_beta);
        if (occupations_beta)
            wf.occupations_beta.assign(occupations_beta, occupations_beta + n_mo_beta);
        w->backend().add_wavefunction_gto(wf);
    });
}

int qvf_set_wavefunction_overlap(qvf_writer* w, const char* section_id, int n_ao,
                                 const double* overlap) {
    return guard(w, [&] {
        if (!overlap || n_ao < 0)
            throw std::runtime_error(
                "qvf_set_wavefunction_overlap: overlap must be non-NULL with n_ao >= 0");
        qvf::Matrix m;
        m.rows = n_ao;
        m.cols = n_ao;
        m.data.assign(overlap, overlap + static_cast<size_t>(n_ao) * static_cast<size_t>(n_ao));
        w->backend().set_wavefunction_overlap(section_id ? section_id : "wf", m);
    });
}

int qvf_set_wavefunction_ecp_json(qvf_writer* w, const char* section_id,
                                  const char* ecp_json) {
    return guard(w, [&] {
        if (!ecp_json)
            throw std::runtime_error("qvf_set_wavefunction_ecp_json: ecp_json is NULL");
        w->backend().set_wavefunction_ecp(section_id ? section_id : "wf",
                                     qvf::Json::raw(ecp_json));
    });
}

int qvf_add_volume(qvf_writer* w, const char* kind, const char* grid_json,
                   const char* dtype, int rank, const int64_t* shape,
                   const void* data, const char* label) {
    return guard(w, [&] {
        if (!kind) throw std::runtime_error("qvf_add_volume: kind is NULL");
        if (std::string(kind) == "volume.series")
            throw std::runtime_error("qvf_add_volume: use qvf_add_volume_series for "
                                     "volume.series");
        const qvf::Tensor t = volume_tensor(dtype, rank, shape, data);
        const std::vector<int64_t> spatial(t.shape.begin(),
                                           t.shape.begin() + std::min<size_t>(3, t.shape.size()));
        w->backend().add_volume(kind, grid_from_json(grid_json, spatial), t, label ? label : "");
    });
}

int qvf_add_volume_series(qvf_writer* w, const char* grid_json,
                          const char* axis_json, const char* dtype,
                          const int64_t* shape4, const void* data,
                          const char* label) {
    return guard(w, [&] {
        if (!axis_json) throw std::runtime_error("qvf_add_volume_series: axis_json is NULL");
        const qvf::Tensor t = volume_tensor(dtype, 4, shape4, data);
        const std::vector<int64_t> spatial(t.shape.begin() + 1, t.shape.end());
        w->backend().add_volume_series(grid_from_json(grid_json, spatial), t,
                                  qvf::Json::raw(axis_json), label ? label : "");
    });
}

// Upper bound on each extent of a wavefunction.bloch call, so that the element
// counts below cannot overflow.
static const int kMaxBlochExtent = 1 << 20;

int qvf_add_wavefunction_bloch(
    qvf_writer* w, int n_shells, const int* shell_center, const int* shell_l,
    const int* shell_nprim, const double* exponents_flat,
    const double* coeffs_flat, int pure, int coeffs_are_libint_normalized,
    int n_k, const double* kpoints_frac, const double* weights_or_null,
    int n_spin, int n_mo, int n_ao, const double* coeffs_re_im,
    const double* energies, const double* occupations) {
    return guard(w, [&] {
        if (n_shells < 0 || (n_shells > 0 && (!shell_center || !shell_l || !shell_nprim ||
                                               !exponents_flat || !coeffs_flat)))
            throw std::runtime_error("qvf_add_wavefunction_bloch: shell arrays must be "
                                     "non-NULL");
        if (n_k <= 0 || n_k > kMaxBlochExtent || n_mo < 0 || n_mo > kMaxBlochExtent ||
            n_ao < 0 || n_ao > kMaxBlochExtent || (n_spin != 1 && n_spin != 2))
            throw std::runtime_error("qvf_add_wavefunction_bloch: need 0 < n_k, 0 <= n_mo, "
                                     "0 <= n_ao (each <= 2^20) and n_spin 1 or 2");
        if (!kpoints_frac || !coeffs_re_im || !energies || !occupations)
            throw std::runtime_error("qvf_add_wavefunction_bloch: kpoints, coefficients, "
                                     "energies and occupations must be non-NULL");
        qvf::WavefunctionBloch wf;
        wf.pure = pure != 0;
        wf.coeffs_are_libint_normalized = coeffs_are_libint_normalized != 0;
        wf.shells = shells_from_arrays(n_shells, shell_center, shell_l, shell_nprim,
                                       exponents_flat, coeffs_flat, pure);
        wf.spin = n_spin == 1 ? "restricted" : "unrestricted";
        for (int k = 0; k < n_k; ++k)
            wf.kpoints.push_back({{kpoints_frac[3 * k], kpoints_frac[3 * k + 1],
                                   kpoints_frac[3 * k + 2]}});
        if (weights_or_null) wf.weights.assign(weights_or_null, weights_or_null + n_k);
        wf.n_mo = n_mo;
        const size_t n_orb = static_cast<size_t>(n_k) * static_cast<size_t>(n_spin) *
                             static_cast<size_t>(n_mo);
        wf.mo_coefficients.assign(coeffs_re_im,
                                  coeffs_re_im + n_orb * static_cast<size_t>(n_ao) * 2);
        wf.energies.assign(energies, energies + n_orb);
        wf.occupations.assign(occupations, occupations + n_orb);
        int64_t basis_ao = 0;
        for (const auto& sh : wf.shells)
            basis_ao += sh.pure ? (2 * sh.l + 1) : ((sh.l + 1) * (sh.l + 2) / 2);
        if (basis_ao != n_ao)
            throw std::runtime_error("qvf_add_wavefunction_bloch: n_ao " + std::to_string(n_ao) +
                                     " differs from the " + std::to_string(basis_ao) +
                                     " AOs of the shells");
        w->backend().add_wavefunction_bloch(wf);
    });
}

int qvf_add_spectrum_xy(qvf_writer* w, const char* kind, const double* freqs,
                        const double* intensities, int n) {
    return guard(w, [&] {
        qvf::Spectrum s;
        s.frequencies.assign(freqs, freqs + n);
        s.intensities.assign(intensities, intensities + n);
        w->backend().add_spectrum(kind, s);
    });
}

int qvf_add_spectrum_json(qvf_writer* w, const char* kind,
                          const char* spectrum_json) {
    return guard(w, [&] {
        w->backend().add_spectrum(kind, qvf::Json::raw(spectrum_json));
    });
}

int qvf_add_atom_properties(qvf_writer* w, const double* mulliken,
                            const double* loewdin, const double* spin_population,
                            int n_atoms) {
    return guard(w, [&] {
        auto vec = [&](const double* p) {
            return p ? std::vector<double>(p, p + n_atoms) : std::vector<double>{};
        };
        w->backend().add_atom_properties(vec(mulliken), vec(loewdin),
                                    vec(spin_population));
    });
}

int qvf_add_atom_property(qvf_writer* w, const char* role,
                          const double* values, int n_atoms) {
    return guard(w, [&] {
        if (!role) throw std::runtime_error("qvf_add_atom_property: role is NULL");
        if (!values || n_atoms <= 0)
            throw std::runtime_error(
                "qvf_add_atom_property: values must be non-NULL with n_atoms > 0");
        w->backend().add_atom_property(role, std::vector<double>(values, values + n_atoms));
    });
}

int qvf_add_atom_vector_property(qvf_writer* w, const char* role,
                                 const double* values_n_atoms_x3, int n_atoms) {
    return guard(w, [&] {
        if (!role)
            throw std::runtime_error("qvf_add_atom_vector_property: role is NULL");
        if (!values_n_atoms_x3 || n_atoms <= 0)
            throw std::runtime_error(
                "qvf_add_atom_vector_property: values must be non-NULL with "
                "n_atoms > 0");
        std::vector<std::array<double, 3>> rows(static_cast<size_t>(n_atoms));
        for (size_t i = 0; i < rows.size(); ++i)
            rows[i] = {{values_n_atoms_x3[3 * i], values_n_atoms_x3[3 * i + 1],
                        values_n_atoms_x3[3 * i + 2]}};
        w->backend().add_atom_vector_property(role, rows);
    });
}

int qvf_set_atom_properties_meta_json(qvf_writer* w, const char* meta_json) {
    return guard(w, [&] {
        if (!meta_json)
            throw std::runtime_error(
                "qvf_set_atom_properties_meta_json: meta_json is NULL");
        w->backend().set_atom_properties_meta(qvf::Json::raw(meta_json));
    });
}

int qvf_add_scf_history_json(qvf_writer* w, const char* iterations_json) {
    return guard(w, [&] {
        w->backend().add_scf_history(qvf::Json::raw(iterations_json));
    });
}

int qvf_add_citations(qvf_writer* w, const char* bibtex) {
    return guard(w, [&] { w->backend().add_citations(bibtex ? bibtex : ""); });
}

int qvf_add_run_record(qvf_writer* w, const char* program,
                       const char* input_text, const char* log_text,
                       const char* program_version, const char* command,
                       int has_exit_status, int exit_status,
                       const char* started_utc, const char* finished_utc,
                       int sequence) {
    return guard(w, [&] {
        qvf::RunRecord record;
        record.program = program ? program : "";
        record.input_text = input_text ? input_text : "";
        record.log_text = log_text ? log_text : "";
        record.program_version = program_version ? program_version : "";
        record.command = command ? command : "";
        record.has_exit_status = has_exit_status != 0;
        record.exit_status = exit_status;
        record.started_utc = started_utc ? started_utc : "";
        record.finished_utc = finished_utc ? finished_utc : "";
        record.sequence = sequence;
        w->backend().add_run_record(record);
    });
}

int qvf_add_job_spec(qvf_writer* w, const char* job_type,
                     const char* method, const char* basis,
                     const char* functional,
                     int has_charge, int charge, int multiplicity,
                     const int* kpoints,
                     const char* tasks_json, const char* options_json) {
    return guard(w, [&] {
        qvf::Json spec = qvf::Json::object();
        spec.set("job_type", job_type ? job_type : "");
        if (method && *method) spec.set("method", method);
        if (basis && *basis) spec.set("basis", basis);
        if (functional && *functional) spec.set("functional", functional);
        if (has_charge) spec.set("charge", charge);
        if (multiplicity >= 1) spec.set("multiplicity", multiplicity);
        if (kpoints)
            spec.set("kpoints", qvf::Json::from_ints(
                std::vector<int>(kpoints, kpoints + 3)));
        if (tasks_json && *tasks_json)
            spec.set("tasks", qvf::Json::raw(tasks_json));
        if (options_json && *options_json)
            spec.set("options", qvf::Json::raw(options_json));
        w->backend().add_job_spec(spec);
    });
}

int qvf_add_job_spec_json(qvf_writer* w, const char* spec_json) {
    return guard(w, [&] {
        if (!spec_json || !*spec_json)
            throw std::runtime_error("qvf_add_job_spec_json: spec_json is empty");
        qvf::Json spec;
        try {
            spec = qvf::Json::parse(spec_json);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                std::string("qvf_add_job_spec_json: spec_json is not valid JSON: ") +
                e.what());
        }
        // Throws unless this is a JobSpecPayload object with a valid job_type
        // and valid typed fields (spec § 5.9).
        w->backend().add_job_spec(spec);
    });
}

int qvf_set_provenance_json(qvf_writer* w, const char* provenance_json) {
    return guard(w, [&] {
        if (!provenance_json)
            throw std::runtime_error("qvf_set_provenance_json: provenance_json is NULL");
        // The writer parses the raw string, checks the typed fields (spec § 3.2)
        // and drops a deprecated top-level "hostname" key; malformed JSON or a
        // wrong field throws -> non-zero.
        w->backend().set_provenance(qvf::Json::raw(provenance_json));
    });
}

int qvf_set_energies_json(qvf_writer* w, const char* energies_json) {
    return guard(w, [&] {
        if (!energies_json || !*energies_json)
            throw std::runtime_error("qvf_set_energies_json: energies_json is empty");
        w->backend().set_energies(qvf::Json::raw(energies_json));
    });
}

// Largest atom count the derivative entry points accept: keeps (3 n)^2
// element counts far from size_t overflow on every platform.
static const int kMaxDerivativeAtoms = 1 << 20;

int qvf_add_forces(qvf_writer* w, int n_atoms, const double* gradient,
                   const double* stress_or_null) {
    return guard(w, [&] {
        if (!gradient || n_atoms <= 0 || n_atoms > kMaxDerivativeAtoms)
            throw std::runtime_error(
                "qvf_add_forces: gradient must be non-NULL with 0 < n_atoms <= 2^20");
        const size_t n = static_cast<size_t>(n_atoms) * 3;
        const qvf::Tensor g = qvf::tensor_f64(std::vector<double>(gradient, gradient + n),
                                              {static_cast<int64_t>(n_atoms), 3});
        if (stress_or_null) {
            const qvf::Tensor st = qvf::tensor_f64(
                std::vector<double>(stress_or_null, stress_or_null + 9), {3, 3});
            w->backend().add_forces(g, &st);
        } else {
            w->backend().add_forces(g);
        }
    });
}

int qvf_add_hessian(qvf_writer* w, int n_atoms, const double* matrix) {
    return guard(w, [&] {
        if (!matrix || n_atoms <= 0 || n_atoms > kMaxDerivativeAtoms)
            throw std::runtime_error(
                "qvf_add_hessian: matrix must be non-NULL with 0 < n_atoms <= 2^20");
        const int64_t dim = 3 * static_cast<int64_t>(n_atoms);
        const size_t count = static_cast<size_t>(dim) * static_cast<size_t>(dim);
        w->backend().add_hessian(
            qvf::tensor_f64(std::vector<double>(matrix, matrix + count), {dim, dim}));
    });
}

int qvf_add_excited_states_json(qvf_writer* w, const char* states_json,
                                const char* meta_json,
                                const char* wavefunction_ref) {
    return guard(w, [&] {
        if (!states_json || !*states_json)
            throw std::runtime_error("qvf_add_excited_states_json: states_json is empty");
        // Parsed and checked by the writer; malformed JSON throws.
        w->backend().add_excited_states(
            qvf::Json::raw(states_json),
            meta_json && *meta_json ? qvf::Json::raw(meta_json) : qvf::Json(),
            wavefunction_ref ? wavefunction_ref : "");
    });
}

int qvf_add_trajectory(qvf_writer* w, const char* metadata_json, int n_frames,
                       int n_atoms, const double* coords_f64,
                       const double* lattice_rows_or_null, int lattice_per_frame,
                       const double* velocities_or_null,
                       const double* gradient_or_null, const double* times_or_null,
                       const double* energies_or_null) {
    return guard(w, [&] {
        if (!metadata_json)
            throw std::runtime_error("qvf_add_trajectory: metadata_json is NULL");
        if (!coords_f64 || n_frames <= 0 || n_atoms <= 0)
            throw std::runtime_error(
                "qvf_add_trajectory: coords must be non-NULL with n_frames > 0 and "
                "n_atoms > 0");
        // Keep 9 * n_frames and n_frames * n_atoms * 3 doubles far from size_t
        // overflow on every platform.
        const uint64_t frames = static_cast<uint64_t>(n_frames);
        const uint64_t per_frame = static_cast<uint64_t>(n_atoms) * 3;
        if (frames * per_frame > (SIZE_MAX / 8) / 2)
            throw std::runtime_error("qvf_add_trajectory: n_frames * n_atoms is too large");
        const size_t nf = static_cast<size_t>(n_frames);
        const size_t nv = nf * static_cast<size_t>(per_frame);
        const std::vector<int64_t> vec_shape = {n_frames, n_atoms, 3};
        // "precision": "float32" in the metadata stores the coordinates as float32.
        const qvf::Json parsed = qvf::Json::parse(metadata_json);
        const bool f32 = parsed.is_object() && parsed.contains("precision") &&
                         parsed.at("precision").is_string() &&
                         parsed.at("precision").as_string() == "float32";
        const std::vector<double> xyz(coords_f64, coords_f64 + nv);
        const qvf::Tensor coords =
            f32 ? qvf::tensor_f32(xyz, vec_shape) : qvf::tensor_f64(xyz, vec_shape);
        qvf::Tensor lattice, velocities, gradient, times, energies;
        qvf::TrajectoryExtras extras;
        if (lattice_rows_or_null) {
            const size_t n = lattice_per_frame ? 9 * nf : 9;
            lattice = qvf::tensor_f64(
                std::vector<double>(lattice_rows_or_null, lattice_rows_or_null + n),
                lattice_per_frame ? std::vector<int64_t>{n_frames, 3, 3}
                                  : std::vector<int64_t>{3, 3});
            extras.lattice = &lattice;
        }
        if (velocities_or_null) {
            velocities = qvf::tensor_f64(
                std::vector<double>(velocities_or_null, velocities_or_null + nv), vec_shape);
            extras.velocities = &velocities;
        }
        if (gradient_or_null) {
            gradient = qvf::tensor_f64(
                std::vector<double>(gradient_or_null, gradient_or_null + nv), vec_shape);
            extras.gradient = &gradient;
        }
        if (times_or_null) {
            times = qvf::tensor_f64(std::vector<double>(times_or_null, times_or_null + nf),
                                    {n_frames});
            extras.times = &times;
        }
        if (energies_or_null) {
            energies = qvf::tensor_f64(
                std::vector<double>(energies_or_null, energies_or_null + nf), {n_frames});
            extras.energies = &energies;
        }
        w->backend().add_trajectory(qvf::Json::raw(metadata_json), coords, extras);
    });
}

int qvf_add_electric_properties_json(qvf_writer* w, const char* properties_json) {
    return guard(w, [&] {
        if (!properties_json || !*properties_json)
            throw std::runtime_error("qvf_add_electric_properties_json: properties_json is empty");
        w->backend().add_electric_properties(qvf::Json::raw(properties_json));
    });
}

int qvf_add_elastic_properties_json(qvf_writer* w, const char* properties_json) {
    return guard(w, [&] {
        if (!properties_json || !*properties_json)
            throw std::runtime_error("qvf_add_elastic_properties_json: properties_json is empty");
        w->backend().add_elastic_properties(qvf::Json::raw(properties_json));
    });
}

int qvf_add_dielectric_properties_json(qvf_writer* w, const char* properties_json) {
    return guard(w, [&] {
        if (!properties_json || !*properties_json)
            throw std::runtime_error(
                "qvf_add_dielectric_properties_json: properties_json is empty");
        w->backend().add_dielectric_properties(qvf::Json::raw(properties_json));
    });
}

int qvf_set_metadata_json(qvf_writer* w, const char* metadata_json) {
    return guard(w, [&] {
        if (!metadata_json) throw std::runtime_error("metadata_json is NULL");
        // Parsed by the writer; malformed JSON or a non-object throws.
        w->backend().set_metadata(qvf::Json::raw(metadata_json));
    });
}

int qvf_annotate_section(qvf_writer* w, const char* id,
                         const char* structure_ref, const char* run_ref,
                         const char* group) {
    return guard(w, [&] {
        if (!id) throw std::runtime_error("section id is NULL");
        w->backend().annotate_section(id, structure_ref ? structure_ref : "",
                                 run_ref ? run_ref : "", group ? group : "");
    });
}

int qvf_add_vendor_section(qvf_writer* w, const char* kind, const char* role,
                           const char* payload_json, int critical) {
    return guard(w, [&] {
        qvf::Json members = qvf::Json::object();
        members.set(role, qvf::Json::raw(payload_json));
        w->backend().add_vendor_section(kind, members, critical != 0);
    });
}

int qvf_set_extensions_json(qvf_writer* w, const char* extensions_json) {
    return guard(w, [&] {
        w->backend().set_extensions(qvf::Json::raw(extensions_json));
    });
}

int qvf_write(qvf_writer* w, const char* path) {
    return guard(w, [&] {
        qvf_c_detail::require(w->legacy_mode, "use qvf_writer_finish for a streaming/editor handle");
        w->backend().write(path ? path : "out.qvf");
    });
}

}  // extern "C"
