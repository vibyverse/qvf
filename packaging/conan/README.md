# Local Conan 2 package (toolkit 0.6.0)

This recipe packages the C/C++ library into a local Conan cache. It does not
upload packages or imply a ConanCenter listing. Use Conan 2.20 or later
within major version 2, Git, CMake 3.16 or later and a C++17 toolchain.
Conan 2.33 is the development verification version.

The recipe uses the [Conan package lifecycle](https://docs.conan.io/2/tutorial/creating_packages/create_your_first_package.html),
with [custom source export](https://docs.conan.io/2/reference/conanfile/methods/export_sources.html)
and the toolkit's [installed CMake configuration](https://docs.conan.io/2/reference/tools/cmake/cmakedeps.html#disable-cmakedeps-for-installed-cmake-configuration-files).
It reads the version from `VERSION`. Export requires a QVF Git checkout with
committed library/recipe files; it copies those exact Git blobs, excluding
untracked files, caches, checkout filters and Git history. Source tarballs
without Git metadata are built with CMake directly. A different `--version`
is refused. A local package revision identifies the exported source bytes;
record the Git commit and Conan revision when adopting development source.

## Create and check locally

From the repository root, use an external cache and test-build directory:

```sh
QVF_CONAN_ROOT="$(mktemp -d)"
export CONAN_HOME="$QVF_CONAN_ROOT/cache"
conan profile detect --name qvf-local
conan create packaging/conan --no-remote -pr:h qvf-local -pr:b qvf-local \
  -s:h compiler.cppstd=17 -c tools.build:jobs=2 \
  -c tools.cmake.cmake_layout:test_folder="$QVF_CONAN_ROOT/test-static"
conan create packaging/conan --no-remote -pr:h qvf-local -pr:b qvf-local \
  -s:h compiler.cppstd=17 -o 'qvf/*:shared=True' -c tools.build:jobs=2 \
  -c tools.cmake.cmake_layout:test_folder="$QVF_CONAN_ROOT/test-shared"
```

Review the detected compiler profile before building. Run native commands
under your environment's native-build lock; the two-job limit also applies
to the small test consumers. Keep the cache while consumers use it. When
finished, remove only this explicitly created temporary directory.

Static linking is the default. `shared=True` uses the existing ABI-versioned
shared library; `fPIC` controls static position-independent code on supported
platforms and is removed for shared libraries or Windows. Toolkit examples,
unit tests, fuzzing and Fortran are disabled in these packages. The included
`test_package` compiles and runs a C99 consumer and a C++ consumer that writes
and reads a small archive. Cross builds compile those consumers and run them
only when Conan's `can_run` permits execution. Local package tests do not
establish the full language/platform or scientific validation matrix.

## Consume the package

Require `qvf/<VERSION>` (pin its recipe revision for reproducibility) in your
Conan consumer and generate `CMakeToolchain` plus `CMakeDeps`. Then use:

```cmake
find_package(qvf CONFIG REQUIRED)
target_link_libraries(your_app PRIVATE qvf::qvf)
```

Conan points CMake to QVF's installed configuration, preserving its C++17
requirement and shared C-ABI import definition. A C-only target still needs
the C++ runtime; set its `LINKER_LANGUAGE` to `CXX`, as the test package does.
The package also retains relocatable `lib/pkgconfig/qvf.pc`, the single
header, `licenses/LICENSE` and `licenses/NOTICE`. With that `pkgconfig`
directory on `PKG_CONFIG_PATH`, compile C++ clients with `-std=c++17` and
`pkg-config --cflags --libs qvf`; compile C translation units as C and use
the C++ driver for linking. Shared-library runtime discovery follows the
consumer's platform or Conan run environment.

The recipe does not package compiler-specific Fortran modules. Fortran users
should follow the [separate binding guide](../../docs/fortran.md). A
[vcpkg overlay port](../vcpkg/README.md) is available separately. Registry
publication and broader platform acceptance remain WP-49 work.
