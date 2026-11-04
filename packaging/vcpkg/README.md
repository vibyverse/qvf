# vcpkg overlay port (toolkit 0.6.0)

`ports/qvf` packages the C/C++ library as a
[vcpkg overlay port](https://learn.microsoft.com/vcpkg/concepts/overlay-ports).
It does not publish anything and does not imply a listing in the vcpkg
curated registry or in any other registry. It needs a current vcpkg
checkout and tool, Git, CMake 3.16 or later and a C++17 toolchain. The
development verification used the vcpkg tool `2026-10-29` with a vcpkg
checkout from 2026-11-04 on x86_64 Linux and GCC.

The port version in `ports/qvf/vcpkg.json` matches `VERSION`. It builds
from one of two sources:

| Source | Selection | Identity |
| --- | --- | --- |
| Tag archive | `QVF_VCPKG_SOURCE_DIR` unset | GitHub archive of tag `v<VERSION>`, checked against the SHA512 in `portfile.cmake` |
| Local source tree | `QVF_VCPKG_SOURCE_DIR=<qvf root>` | A QVF Git checkout or extracted source archive whose `VERSION` matches the port |

A local tree with a different `VERSION` is refused. The port copies only
`LICENSE`, `NOTICE`, `VERSION` and the library inputs under `cpp/`
(`CMakeLists.txt`, `cmake/`, `include/`, `src/`, `qvf_single.hpp`) into its
buildtree, so the checkout is never configured or written and local build
directories stay out. Unlike the [Conan recipe](../conan/README.md), local
mode reads the working tree, including uncommitted edits. vcpkg does not hash
local sources into the package ABI, so always pass `--binarysource=clear`
in local mode and record the Git commit and the state of the working tree
when adopting development source.

Static or shared linkage follows the triplet (`VCPKG_LIBRARY_LINKAGE`).
Shared builds use the existing ABI-versioned library (`libqvf.so.1` on
Linux). Toolkit examples, unit tests, fuzzing and Fortran are disabled.

## Install and check locally

Obtain vcpkg outside the source tree or in the ignored `build/` directory.
From the repository root:

```sh
git clone https://github.com/microsoft/vcpkg.git build/tools/vcpkg
build/tools/vcpkg/bootstrap-vcpkg.sh -disableMetrics
export VCPKG_ROOT="$PWD/build/tools/vcpkg"

# Current checkout, static and shared
QVF_VCPKG_SOURCE_DIR="$PWD" "$VCPKG_ROOT/vcpkg" install qvf \
  --overlay-ports=packaging/vcpkg/ports --binarysource=clear \
  --x-install-root="$PWD/build/vcpkg-installed" --triplet x64-linux
QVF_VCPKG_SOURCE_DIR="$PWD" "$VCPKG_ROOT/vcpkg" install qvf \
  --overlay-ports=packaging/vcpkg/ports --binarysource=clear \
  --x-install-root="$PWD/build/vcpkg-installed" --triplet x64-linux-dynamic

# Released tag archive, in a separate install root
"$VCPKG_ROOT/vcpkg" install qvf --overlay-ports=packaging/vcpkg/ports \
  --binarysource=clear --x-install-root="$PWD/build/vcpkg-installed-tag" \
  --triplet x64-linux
```

Pick the triplet of your platform (for example `arm64-osx` or
`x64-windows`). Run native commands under your environment's native-build
lock. Local and tag builds have the same port version, so keep them in
separate install roots.

Then build the C99 and C++ consumers from the Conan `test_package` sources
out of tree. Copy `consume.c`, `consume.cpp` and `CMakeLists.txt` from
[`../conan/test_package`](../conan/test_package) to a scratch directory
outside the checkout and run:

```sh
cmake -S "$CONSUMER" -B "$CONSUMER/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=OFF -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DVCPKG_INSTALLED_DIR="<repository>/build/vcpkg-installed"
cmake --build "$CONSUMER/build"
"$CONSUMER/build/consume_c" && "$CONSUMER/build/consume_cpp"
```

Both programs write and read a small archive and exit with status 0. The
consumer CMake project requests C++11, so a successful build also shows that
`qvf::qvf` raises its users to C++17. For the shared triplet, check that the
programs resolve `libqvf.so.1` from the vcpkg install root.

To refresh the tag archive hash after a new release tag, download
`https://github.com/vibyverse/qvf/archive/v<VERSION>.tar.gz`, compute
`sha512sum` and update `SHA512` together with the `version` field. The
release owner does this after the tag exists. Never point the port at a
moving branch.

## Consume the package

With the vcpkg toolchain file (manifest or classic mode) use:

```cmake
find_package(qvf CONFIG REQUIRED)
target_link_libraries(your_app PRIVATE qvf::qvf)
```

The installed configuration in `share/qvf` preserves the C++17 requirement
and, for shared triplets, the `QVF_SHARED=1` import definition. A C-only
target still needs the C++ runtime. Set its `LINKER_LANGUAGE` to `CXX`, as
the consumer does. The port also installs `lib/pkgconfig/qvf.pc` (fixed up by
vcpkg), the single header `include/qvf/qvf_single.hpp` and
`share/qvf/copyright` with LICENSE and NOTICE. With that `pkgconfig`
directory on `PKG_CONFIG_PATH`, compile C++ clients with `-std=c++17` and
`pkg-config --cflags --libs qvf`. Compile C translation units as C and link
with the C++ driver.

The port does not package compiler-specific Fortran modules. Fortran users
should follow the [separate binding guide](../../docs/fortran.md). A
submission to the curated vcpkg registry or a custom Git registry, and
platform acceptance beyond the checked Linux triplets, remain WP-49 work.
