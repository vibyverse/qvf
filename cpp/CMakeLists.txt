cmake_minimum_required(VERSION 3.16)
project(qvf_writer LANGUAGES CXX C VERSION 0.6.0)

# Standalone development includes checks and examples. Embedded consumers only
# build the library unless they explicitly opt into those targets.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(QVF_STANDALONE ON)
else()
  set(QVF_STANDALONE OFF)
endif()
option(QVF_BUILD_TESTS "Build QVF tests and register them with CTest" ${QVF_STANDALONE})
option(QVF_BUILD_EXAMPLES "Build QVF example programs" ${QVF_STANDALONE})

# Zero external dependencies — only the C++17 standard library.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  # The sources are UTF-8 (comments and one string literal contain non-ASCII
  # text). Without /utf-8 MSVC reads them in the active ANSI code page, which
  # raises warning C4819 on some locales (an error under /WX).
  add_compile_options(/utf-8)
endif()

# --- fuzzing (off by default) -------------------------------------------------
# -DQVF_FUZZ=ON (clang only) builds the libFuzzer targets in cpp/fuzz:
# zip_reader_fuzzer over raw archive bytes (ZIP reader and validator) and
# manifest_fuzzer over manifest.json and member payloads in a consistent ZIP.
# Everything in such a build is compiled with ASan + UBSan and coverage
# instrumentation, so the fuzzers reach into the library. Where clang ships no libFuzzer runtime (Apple
# clang), the targets link fuzz/standalone_fuzz_main.cpp instead.
option(QVF_FUZZ "Build the fuzz targets in cpp/fuzz (clang)" OFF)
if(QVF_FUZZ)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "QVF_FUZZ=ON needs clang (-fsanitize=fuzzer)")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-sanitize-recover=undefined
                      -fno-omit-frame-pointer -g -fsanitize=fuzzer-no-link)
  add_link_options(-fsanitize=address,undefined)
endif()

# --- library ---------------------------------------------------------------
# src/qvf.cpp is the C++ writer; src/qvf_c.cpp is the extern "C" ABI wrapper;
# src/qvf_reader.cpp is the C++ reader (manifest, members, validation);
# src/qvf_container.cpp is the ZIP layer both use (qvf/container.hpp).
option(QVF_BUILD_SHARED "Build the shared C/C++ library (ABI 1)" OFF)
if(QVF_BUILD_SHARED)
  set(QVF_LIBRARY_TYPE SHARED)
else()
  set(QVF_LIBRARY_TYPE STATIC)
endif()
add_library(qvf ${QVF_LIBRARY_TYPE} src/qvf.cpp src/qvf_c.cpp src/qvf_c_facade.cpp src/qvf_reader.cpp
                       src/qvf_container.cpp src/qvf_reader_facade.cpp
                       src/qvf_validation.cpp src/qvf_tensor_view.cpp
                       src/qvf_writer.cpp src/qvf_editor.cpp)
include(GNUInstallDirs)
target_compile_features(qvf PUBLIC cxx_std_17)
target_compile_definitions(qvf PRIVATE QVF_TOOLKIT_VERSION="${PROJECT_VERSION}")
if(QVF_BUILD_SHARED)
  target_compile_definitions(qvf PUBLIC QVF_SHARED=1 PRIVATE QVF_BUILDING_LIBRARY=1)
  set_target_properties(qvf PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1
                        WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
target_include_directories(qvf
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)   # internal sha256/zip headers
if(MSVC)
  target_compile_options(qvf PRIVATE /W4)
else()
  target_compile_options(qvf PRIVATE -Wall -Wextra)
endif()

# Alias so downstream `target_link_libraries(foo qvf::qvf)` works when this is
# added via add_subdirectory().
add_library(qvf::qvf ALIAS qvf)

# Optional ISO_C_BINDING facade; the C++ library retains its usual dependencies.
option(QVF_FORTRAN "Build the Fortran modules and optional examples/tests" OFF)
if(QVF_FORTRAN)
  enable_language(Fortran)
  set(QVF_FORTRAN_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
  set(QVF_FORTRAN_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/qvf/fortran/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
  file(MAKE_DIRECTORY "${QVF_FORTRAN_MODULE_DIR}")
  add_library(qvf_fortran STATIC ../fortran/qvf_c.f90 ../fortran/qvf.f90)
  add_library(qvf::fortran ALIAS qvf_fortran)
  set_target_properties(qvf_fortran PROPERTIES EXPORT_NAME fortran
    Fortran_MODULE_DIRECTORY "${QVF_FORTRAN_MODULE_DIR}")
  target_link_libraries(qvf_fortran PUBLIC qvf)
  target_include_directories(qvf_fortran PUBLIC
    $<BUILD_INTERFACE:${QVF_FORTRAN_MODULE_DIR}>
    $<INSTALL_INTERFACE:${QVF_FORTRAN_INSTALL_DIR}>)
  install(TARGETS qvf_fortran EXPORT qvfTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY "${QVF_FORTRAN_MODULE_DIR}/" DESTINATION "${QVF_FORTRAN_INSTALL_DIR}"
          FILES_MATCHING PATTERN "*.mod")
  install(FILES ../fortran/qvf_c.f90 ../fortran/qvf.f90 DESTINATION ${CMAKE_INSTALL_DATADIR}/qvf/fortran)
  if(QVF_BUILD_EXAMPLES OR QVF_BUILD_TESTS)
    add_executable(write_sections_fortran ../fortran/examples/write_sections.f90)
    target_link_libraries(write_sections_fortran PRIVATE qvf_fortran)
  endif()
  if(QVF_BUILD_TESTS)
    enable_testing()
    add_executable(qvf_fortran_tests ../fortran/tests/test_binding.f90 ../fortran/tests/abi_sizes.c)
    target_link_libraries(qvf_fortran_tests PRIVATE qvf_fortran)
    add_test(NAME qvf_fortran_binding COMMAND qvf_fortran_tests "${CMAKE_CURRENT_BINARY_DIR}/fortran-refusals.qvf")
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME qvf_fortran_layout COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../fortran/tests/check_layout.py"
      $<TARGET_FILE:write_sections_fortran> "${CMAKE_CURRENT_BINARY_DIR}/fortran-sections.qvf")
  endif()
endif()

# --- example ---------------------------------------------------------------
if(QVF_BUILD_EXAMPLES)
add_executable(write_h2o examples/write_h2o.cpp)
target_link_libraries(write_h2o PRIVATE qvf)

add_executable(write_all_kinds examples/write_all_kinds.cpp)
target_link_libraries(write_all_kinds PRIVATE qvf)

# The typed reader CLI; keep read_qvf as the compatibility executable.
add_executable(qvf-inspect examples/read_qvf.cpp)
target_link_libraries(qvf-inspect PRIVATE qvf)
add_executable(read_qvf examples/read_qvf.cpp)
target_link_libraries(read_qvf PRIVATE qvf)

# A pure-C example driving the extern "C" ABI. Linked with the C++ driver so the
# C++ runtime is pulled in.
add_executable(write_h2o_c examples/write_h2o_c.c)
target_link_libraries(write_h2o_c PRIVATE qvf)
set_target_properties(write_h2o_c PROPERTIES LINKER_LANGUAGE CXX)

add_executable(write_all_kinds_c examples/write_all_kinds_c.c)
target_link_libraries(write_all_kinds_c PRIVATE qvf)
set_target_properties(write_all_kinds_c PROPERTIES LINKER_LANGUAGE CXX)
add_executable(read_qvf_c examples/read_qvf_c.c)
target_link_libraries(read_qvf_c PRIVATE qvf)
set_target_properties(read_qvf_c PROPERTIES LINKER_LANGUAGE CXX)

# Single-header amalgamation demo. Compiled against ONLY qvf_single.hpp (the
# amalgamated header lives at the cpp/ root) — it does NOT see include/ or src/
# and does NOT link the qvf library, so building it proves the single header is
# fully self-contained (STB-style, with QVF_IMPLEMENTATION in this one TU).
add_executable(write_single_header examples/write_single_header.cpp)
target_include_directories(write_single_header PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# CTest needs the runner even when example programs are disabled.
if(QVF_BUILD_EXAMPLES OR QVF_BUILD_TESTS)
  add_executable(run_conformance examples/run_conformance.cpp)
  target_link_libraries(run_conformance PRIVATE qvf)
endif()

# --- tests -----------------------------------------------------------------
if(QVF_BUILD_TESTS)
enable_testing()
add_executable(qvf_c_api_tests tests/test_c_api.c)
target_link_libraries(qvf_c_api_tests PRIVATE qvf)
set_target_properties(qvf_c_api_tests PROPERTIES LINKER_LANGUAGE CXX)
add_test(NAME qvf_c_api_tests COMMAND qvf_c_api_tests
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_executable(qvf_tests tests/test_writer.cpp)
target_link_libraries(qvf_tests PRIVATE qvf)
target_include_directories(qvf_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME qvf_tests COMMAND qvf_tests)

add_executable(qvf_writer_facade_tests tests/test_writer_facade.cpp)
target_link_libraries(qvf_writer_facade_tests PRIVATE qvf)
add_test(NAME qvf_writer_facade_tests
         COMMAND qvf_writer_facade_tests ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/corpus
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# The container layer (ZIP64, data descriptors, streaming, caps) and the
# streaming / updating writer entry points. Scratch files go to the build tree;
# the corpus' zip64.qvf and data_descriptor.qvf are updated too when present.
add_executable(qvf_container_tests tests/test_container.cpp)
target_link_libraries(qvf_container_tests PRIVATE qvf)
target_include_directories(qvf_container_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME qvf_container_tests
         COMMAND qvf_container_tests ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/corpus
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Validator hardening: malformed archives built in memory, plus the negative
# conformance corpus (every archive must yield exactly its expected finding
# codes). A standalone copy of cpp/ without ../conformance skips the latter.
add_executable(qvf_reader_tests tests/test_reader.cpp)
target_link_libraries(qvf_reader_tests PRIVATE qvf)
target_include_directories(qvf_reader_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME qvf_reader_tests COMMAND qvf_reader_tests)
# Lazy typed access, bounded tensor streaming and cumulative validation stages.
add_executable(qvf_reader_facade_tests tests/test_reader_facade.cpp)
target_link_libraries(qvf_reader_facade_tests PRIVATE qvf)
target_include_directories(qvf_reader_facade_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME qvf_reader_facade_tests
         COMMAND qvf_reader_facade_tests ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/negative)
# Drives the generated JSON Schema checker for the parity test against the
# Python reference (python/tests/test_schema_checker_parity.py).
add_executable(qvf_schema_check tests/schema_check_tool.cpp)
target_link_libraries(qvf_schema_check PRIVATE qvf)
target_include_directories(qvf_schema_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME qvf_negative_corpus
         COMMAND qvf_reader_tests --negative
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/negative)
set_tests_properties(qvf_negative_corpus PROPERTIES SKIP_RETURN_CODE 77)
# The typed manifest model (qvf/model.hpp) round-trips every conformance
# manifest losslessly.
add_test(NAME qvf_model_roundtrip
         COMMAND qvf_reader_tests --model-roundtrip
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/corpus
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/negative
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/legacy)
set_tests_properties(qvf_model_roundtrip PROPERTIES SKIP_RETURN_CODE 77)

# The whole conformance suite through the C++ reader: every positive archive
# must decode to its expected.json values, every must-fail archive must be
# rejected with its codes, and the legacy archives (qvf_version 2) judged as
# their expected.json says: valid, with the warning W_VERSION_WITHDRAWN. A
# standalone copy of cpp/ without ../conformance skips it.
add_test(NAME qvf_conformance
         COMMAND run_conformance
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/corpus
                 --negative ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/negative
                 --legacy ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/legacy)
set_tests_properties(qvf_conformance PROPERTIES SKIP_RETURN_CODE 77)
# The same without the generated JSON Schema check: the behaviour of a
# validator without a schema engine (W_SCHEMA_SKIPPED) stays testable.
add_test(NAME qvf_conformance_no_schema
         COMMAND run_conformance
                 ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/corpus
                 --negative ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/negative
                 --legacy ${CMAKE_CURRENT_SOURCE_DIR}/../conformance/legacy
                 --no-schema-check)
set_tests_properties(qvf_conformance_no_schema PROPERTIES SKIP_RETURN_CODE 77)
endif()

# --- fuzz targets (QVF_FUZZ=ON) ---------------------------------------------
if(QVF_FUZZ)
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=fuzzer")
  set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=fuzzer")
  check_cxx_source_compiles(
    "#include <cstddef>\n#include <cstdint>\nextern \"C\" int LLVMFuzzerTestOneInput(const uint8_t*, size_t) { return 0; }"
    QVF_HAVE_LIBFUZZER)
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LINK_OPTIONS)
  foreach(fuzzer zip_reader_fuzzer manifest_fuzzer)
    add_executable(${fuzzer} fuzz/${fuzzer}.cpp)
    target_link_libraries(${fuzzer} PRIVATE qvf)
    if(QVF_HAVE_LIBFUZZER)
      target_compile_options(${fuzzer} PRIVATE -fsanitize=fuzzer)
      target_link_options(${fuzzer} PRIVATE -fsanitize=fuzzer)
    else()
      target_sources(${fuzzer} PRIVATE fuzz/standalone_fuzz_main.cpp)
    endif()
  endforeach()
  if(NOT QVF_HAVE_LIBFUZZER)
    message(WARNING "QVF_FUZZ: this clang has no libFuzzer runtime; the fuzz targets "
                    "link fuzz/standalone_fuzz_main.cpp (replay and random mutation, "
                    "no coverage feedback)")
  endif()
endif()

# --- install / export ------------------------------------------------------
# Lets adopters consume the library with `find_package(qvf CONFIG)` or
# FetchContent + `target_link_libraries(app PRIVATE qvf::qvf)`.
include(CMakePackageConfigHelpers)

install(TARGETS qvf EXPORT qvfTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/qvf DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# The single-header amalgamation ships alongside for header-only consumers.
install(FILES qvf_single.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qvf)
install(EXPORT qvfTargets
        NAMESPACE qvf::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qvf
        FILE qvfTargets.cmake)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/qvfConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qvfConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/qvfConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qvf)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/qvfConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/qvfConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qvf)

# Keep the usual relative GNUInstallDirs layout relocatable. Explicit absolute
# library/include destinations retain their requested locations.
file(RELATIVE_PATH QVF_PC_PREFIX_RELATIVE
     "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig" "${CMAKE_INSTALL_PREFIX}")
foreach(dir LIBDIR INCLUDEDIR)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_${dir}}")
    set(QVF_PC_${dir} "${CMAKE_INSTALL_${dir}}")
  else()
    set(QVF_PC_${dir} "\${prefix}/${CMAKE_INSTALL_${dir}}")
  endif()
endforeach()
set(QVF_PC_DEFINES "")
if(QVF_BUILD_SHARED)
  set(QVF_PC_DEFINES " -DQVF_SHARED=1")
endif()
configure_file(cmake/qvf.pc.in qvf.pc @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qvf.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Also usable straight from a build tree via FetchContent (the qvf::qvf alias
# above resolves without an install step).
export(EXPORT qvfTargets NAMESPACE qvf::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/qvfTargets.cmake")
