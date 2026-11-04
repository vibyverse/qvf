# Fortran binding (toolkit 0.5.0)

The optional Fortran modules use the existing C ABI 1. `use qvf` supplies
checked convenience writers and re-exports the raw bindings from `qvf_c`.
No format, schema, registry or finding code changes are involved.

## Build and link

```sh
cmake -S cpp -B build -DQVF_FORTRAN=ON \
  -DQVF_BUILD_TESTS=OFF -DQVF_BUILD_EXAMPLES=ON
cmake --build build --parallel 2
./build/write_sections_fortran example.qvf
cmake --install build --prefix "$PWD/install"
```

`QVF_FORTRAN` defaults to `OFF`; ordinary C/C++ consumers need no Fortran
compiler. Enabling it requires a Fortran compiler with `ISO_C_BINDING` and
Fortran 2008 allocatable-character support. See
[checked compilers](#checked-compilers) for the combinations that have run.
No CI is needed for an ordinary development push.
Fortran examples follow `QVF_BUILD_EXAMPLES`; tests follow `QVF_BUILD_TESTS`.
With tests enabled, choose a `Python3_EXECUTABLE` with NumPy and jsonschema
installed, then run `ctest --test-dir build -R '^qvf_fortran_'`.

## Checked compilers

Only the combinations below have been run. Other compilers, versions and
platforms (macOS, Windows, Intel ifx, NVIDIA nvfortran, Cray, AArch64) are
unchecked. They may work, but no claim is made.

| Fortran compiler | Platform | Checks |
|---|---|---|
| GNU gfortran 11.5.0 | Ubuntu 24.04, x86_64 | focused development checks |
| GNU gfortran 12.4.0 | Ubuntu 24.04, x86_64 | focused development checks |
| GNU gfortran 12 | Linux | toolkit 0.5.0 candidate (all 12 CTests) |
| GNU gfortran 13.3.0 | Ubuntu 24.04, x86_64 | focused development checks, ASan/UBSan |
| GNU gfortran 14.2.0 | Ubuntu 24.04, x86_64 | focused development checks |
| GNU gfortran 16 | Linux | focused implementation checks for 0.5.0 |
| LLVM flang 18.1.3 | Ubuntu 24.04, x86_64 | focused development checks |
| LLVM flang 19.1.1 | Ubuntu 24.04, x86_64 | focused development checks |
| LLVM flang 20.1.2 | Ubuntu 24.04, x86_64 | focused development checks |

The focused development checks used the Ubuntu packages (`gfortran-N`,
`flang-N`) with GCC 13.3.0 as the C/C++ compiler, CMake 3.28.3 and Python
3.13 with NumPy and jsonschema. For each compiler they cover:

- Release and Debug builds of `qvf::fortran`, the example and the tests with
  warnings as errors: `-Wall -Wextra -pedantic -std=f2008 -fcheck=all -Werror`
  for gfortran, `-pedantic -Werror` for flang. No diagnostics were reported.
- Both `qvf_fortran_` CTests: descriptor sizes against the C ABI, ownership,
  refusals, exact error message copies, and Python validation (zero
  findings) with array layout read back through the Python `Reader`.
- Installation and an external `find_package(qvf)` consumer built with the
  same compiler, whose output validates with zero findings.

The gfortran 13 sanitizer run builds the C++ library and the Fortran sources
with `-fsanitize=address,undefined`. LeakSanitizer reports the allocatable
variables of the test main programs, which are implicitly saved and never
freed before exit. Those reports are not library leaks, so that run sets
`ASAN_OPTIONS=detect_leaks=0`.

For an installed or embedded consumer:

```cmake
project(example LANGUAGES CXX Fortran)
find_package(qvf CONFIG REQUIRED)  # or add_subdirectory(qvf/cpp)
add_executable(example example.f90)
target_link_libraries(example PRIVATE qvf::fortran)
```

The target carries the module directory and C++/Fortran link dependencies.
Set `CMAKE_PREFIX_PATH` to the installation prefix. Compiled `.mod` files are
compiler-specific: build the library and consumer with compatible compilers.
Installed modules live beneath `include/qvf/fortran/<compiler>-<version>`;
module sources are also installed under `share/qvf/fortran`. `qvf.pc` describes
the C/C++ library; Fortran consumers should use the CMake target.

## Ownership and errors

```fortran
use qvf
implicit none
type(c_ptr) :: writer = c_null_ptr
character(len=:), allocatable :: message
integer(c_int) :: status
status = qvf_open(writer, 'example.qvf', 'my-code', '1', 'example', message)
! Check status after every operation; print message on failure.
if (status /= QVF_OK) error stop 'open failed'
status = qvf_structure(writer, [1_c_int], &
                      reshape([0.0_c_double,0.0_c_double,0.0_c_double],[1,3]), message)
if (status == QVF_OK) status = qvf_finish(writer, message)
! Closing releases ownership; it never implicitly finishes an archive.
status = qvf_close(writer)
```

Initialize the handle to `c_null_ptr`. `qvf_open` refuses an already owned
handle, so a second open cannot leak it. `qvf_close` destroys the handle and
sets it to null; repeated close is harmless. Do not copy owning handles,
destroy borrowed handles, or call one handle concurrently. Always finish
successfully before close to obtain a complete archive. The destination is
opened by the C writer and may replace an existing file; choose it explicitly.
A failed write may leave an incomplete file which must not be used as QVF.

Convenience functions return `integer(c_int)`, zero (`QVF_OK`) on success.
An optional allocatable `message` receives a copied error string, including
wrapper shape/text errors. Its length is exactly that of the copied text.
`qvf_error_message(handle)` copies the last C
writer error; without a handle it copies the C global error. Wrapper refusal
messages are returned through `message`, not written into the C error buffer.
Text arguments trim trailing blanks and reject embedded NUL bytes. Supply
UTF-8 bytes in default-kind Fortran strings; no encoding conversion occurs.

## Arrays and units

The wrappers accept assumed-shape arrays, including strided slices, and copy
values explicitly into canonical row-major order. Input arrays remain owned
by the caller and are borrowed only for that call. Temporary buffers require
memory proportional to the supplied arrays; use the raw streaming API for
larger data. Empty arrays, non-finite values, inconsistent extents and arrays
with more than `huge(0_c_int)` elements are refused before crossing the ABI.
Non-finite values are refused without raising an IEEE invalid exception, so
callers that enable halting (for example `-ffpe-trap=invalid`) receive a
refusal rather than a signal.

| Function | Array dimensions in the Fortran call | Meaning |
|---|---|---|
| `qvf_structure` | `positions(atom,xyz)`, atomic numbers | Molecular positions in angstrom; symbols derived from atomic numbers 1–118 |
| `qvf_volume` | `values(x,y,z)` | Scalar float64 volume; Z varies fastest on disk; grid JSON uses the established bohr convention |
| `qvf_trajectory` | `coords(frame,atom,xyz)` | Angstrom coordinates; metadata JSON describes the same atoms |
| `qvf_gto` | `mo(orbital,ao)` | Restricted GTO coefficient rows in canonical AO order |

`qvf_trajectory` accepts optional `lattice(1,3,3)` for a fixed cell or
`lattice(n_frames,3,3)` for changing cells, with lattice **row vectors in
angstrom**. `velocities` and `gradients` use the same shape as `coords`:
angstrom/fs and the energy gradient dE/dR in hartree/bohr, respectively.
`times` and `energies` have one value per frame, in fs and hartree. Metadata
must declare `pbc` when a lattice is supplied. No force-sign or unit conversion
is implicit.

`qvf_gto` takes parallel `centers`, `angular_momenta`, `nprim` arrays and
concatenated primitive `exponents`, `coefficients`. Centers are **zero-based**
structure atom indices even in Fortran. Angular momenta are 0–6 in this
wrapper. `pure` defaults true (spherical); `libint_normalized` defaults false
and has the same convention as the C API. Optional `energies` (hartree) and
`occupations` (0–2) have one entry per MO. A primitive-count mismatch is
rejected before C can read any pointer. Neither wrappers nor examples perform
or certify a scientific calculation.

Convenience methods use the C API's existing default section IDs, so repeated
same-kind writes may fail for duplicate IDs. The raw facade supports explicit
IDs, periodic structures, all canonical kinds, unrestricted wavefunctions,
readers, validation, streaming and editing.

## Raw ABI

`scripts/gen_fortran.py` generates `fortran/qvf_c.f90` from the public C
headers. Run it after C ABI changes; `--check` detects stale bindings. All 161
current exports are declared as `c_<C function name>`, with positional
arguments in header order. `bind(C)` structures have a `_t` suffix to avoid
Fortran's case-insensitive type/function name collisions. Constants retain
the C names.

Every C pointer parameter is `type(c_ptr), value`, including `T**` output
parameters: pass `c_loc(output)` where `output` is a `target` variable of the
appropriate interoperable type. `const` ownership/lifetime rules remain
exactly those in [the C API guide](c_api.md) and the public headers. Pass
`c_null_ptr` for absent optional data. Callbacks use `type(c_funptr), value`
and `c_funloc` of a matching `bind(C)` procedure. Raw calls require correctly
sized, contiguous buffers in C row-major order and NUL-terminated strings.
`qvf_c_string` creates such a string array; give its receiving allocatable
array `target` and keep it alive through the call.

C enums/status and the small unsigned ABI version use `integer(c_int)`;
`size_t` uses `integer(c_size_t)`. The 64-bit unsigned limit/count fields use
`integer(c_int64_t)` because Fortran has no portable unsigned kind: use only
nonnegative values representable in that kind. Raw bindings add no range,
shape, lifetime or byte-order checks. Prefer convenience wrappers where they
cover the operation, and the C facade documentation elsewhere.
