# Overlay port for the QVF C/C++ library. It is not a vcpkg registry entry.
#
# Source selection:
#   QVF_VCPKG_SOURCE_DIR unset: the GitHub archive of tag v${VERSION},
#                               verified by the SHA512 below.
#   QVF_VCPKG_SOURCE_DIR set:   a local QVF source tree (Git checkout or
#                               extracted source archive). Its VERSION must
#                               match this port. vcpkg does not hash local
#                               sources into the package ABI, so build such
#                               trees with --binarysource=clear.
if(DEFINED ENV{QVF_VCPKG_SOURCE_DIR} AND NOT "$ENV{QVF_VCPKG_SOURCE_DIR}" STREQUAL "")
    file(REAL_PATH "$ENV{QVF_VCPKG_SOURCE_DIR}" QVF_LOCAL_ROOT)
    if(NOT EXISTS "${QVF_LOCAL_ROOT}/cpp/CMakeLists.txt" OR NOT EXISTS "${QVF_LOCAL_ROOT}/VERSION")
        message(FATAL_ERROR "QVF_VCPKG_SOURCE_DIR must name a QVF source tree: ${QVF_LOCAL_ROOT}")
    endif()
    file(READ "${QVF_LOCAL_ROOT}/VERSION" QVF_LOCAL_VERSION)
    string(STRIP "${QVF_LOCAL_VERSION}" QVF_LOCAL_VERSION)
    if(NOT QVF_LOCAL_VERSION STREQUAL VERSION)
        message(FATAL_ERROR "QVF source VERSION ${QVF_LOCAL_VERSION} does not match port version ${VERSION}")
    endif()
    # Copy only the library inputs into the buildtree. The checkout itself is
    # never configured or written, and local build directories stay out.
    set(SOURCE_PATH "${CURRENT_BUILDTREES_DIR}/src/local")
    file(REMOVE_RECURSE "${SOURCE_PATH}")
    file(MAKE_DIRECTORY "${SOURCE_PATH}/cpp")
    file(COPY "${QVF_LOCAL_ROOT}/LICENSE" "${QVF_LOCAL_ROOT}/NOTICE" "${QVF_LOCAL_ROOT}/VERSION"
         DESTINATION "${SOURCE_PATH}")
    file(COPY "${QVF_LOCAL_ROOT}/cpp/CMakeLists.txt" "${QVF_LOCAL_ROOT}/cpp/qvf_single.hpp"
              "${QVF_LOCAL_ROOT}/cpp/cmake" "${QVF_LOCAL_ROOT}/cpp/include" "${QVF_LOCAL_ROOT}/cpp/src"
         DESTINATION "${SOURCE_PATH}/cpp")
    message(STATUS "Using local QVF source ${QVF_LOCAL_ROOT}")
else()
    vcpkg_from_github(
        OUT_SOURCE_PATH SOURCE_PATH
        REPO vibyverse/qvf
        REF "v${VERSION}"
        SHA512 0259b483080c249a2a7fd1db3a44e629a7a9d5aafc4820ecece36f014bae2685f4a22dfae1fcdb8660a22ede4b8e7cbe826869fe70ab52c1a1cbd02005063370
        HEAD_REF main
    )
endif()

string(COMPARE EQUAL "${VCPKG_LIBRARY_LINKAGE}" "dynamic" QVF_BUILD_SHARED)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}/cpp"
    OPTIONS
        -DQVF_BUILD_SHARED=${QVF_BUILD_SHARED}
        -DQVF_BUILD_TESTS=OFF
        -DQVF_BUILD_EXAMPLES=OFF
        -DQVF_FORTRAN=OFF
        -DQVF_FUZZ=OFF
)
vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME qvf CONFIG_PATH lib/cmake/qvf)
vcpkg_fixup_pkgconfig()
vcpkg_copy_pdbs()

file(REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
    "${CURRENT_PACKAGES_DIR}/debug/share"
)

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE" "${SOURCE_PATH}/NOTICE")
