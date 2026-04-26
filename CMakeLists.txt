cmake_minimum_required(VERSION 3.16)
project(WonderVol VERSION 1.0.0 LANGUAGES C)

# Options
option(WV_BUILD_TESTS      "Build tests"                       ON)
option(WV_BUILD_BENCHMARKS "Build benchmarks"                  ON)
option(WV_BUILD_EXAMPLES   "Build examples"                    ON)
option(WV_ENABLE_OPENMP    "Enable OpenMP parallelization"     ON)
option(WV_STRICT_TRAP      "Enable strict trapping on violations" OFF)

# Language standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Header-only library
add_library(wondervol INTERFACE)
add_library(WonderVol::wondervol ALIAS wondervol)

target_include_directories(wondervol INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Feature test macros for aligned_alloc
target_compile_definitions(wondervol INTERFACE _GNU_SOURCE=1)

# Platform detection and SIMD flags
include(CheckCCompilerFlag)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    target_compile_options(wondervol INTERFACE -mcpu=native)
    target_compile_definitions(wondervol INTERFACE WV_TARGET_ARM64=1)
    message(STATUS "WonderVol: ARM64 detected, using -mcpu=native (NEON enabled)")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    check_c_compiler_flag("-mavx512f" WV_HAS_AVX512F)
    if(WV_HAS_AVX512F)
        target_compile_options(wondervol INTERFACE -mavx512f -mavx512dq)
        message(STATUS "WonderVol: AVX512F enabled")
    else()
        check_c_compiler_flag("-mavx2" WV_HAS_AVX2)
        if(WV_HAS_AVX2)
            target_compile_options(wondervol INTERFACE -mavx2 -mfma)
            message(STATUS "WonderVol: AVX2+FMA enabled")
        else()
            message(STATUS "WonderVol: No SIMD detected, using scalar fallback")
        endif()
    endif()
else()
    message(STATUS "WonderVol: Unknown arch ${CMAKE_SYSTEM_PROCESSOR}, using scalar fallback")
endif()

# OpenMP
if(WV_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        target_link_libraries(wondervol INTERFACE OpenMP::OpenMP_C)
    endif()
endif()

# Compile options
if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(wondervol INTERFACE
        -O3
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# Strict trap
if(WV_STRICT_TRAP)
    target_compile_definitions(wondervol INTERFACE WV_STRICT_TRAP=1)
endif()

# Link math library
target_link_libraries(wondervol INTERFACE m)

# Testing
if(WV_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Benchmarks
if(WV_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Examples
if(WV_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Convenience target: run all tests
if(WV_BUILD_TESTS)
    add_custom_target(check
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS test_bsm test_iv test_dispersion test_platform
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS wondervol
    EXPORT WonderVolTargets
)

install(EXPORT WonderVolTargets
    FILE WonderVolTargets.cmake
    NAMESPACE WonderVol::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/WonderVol
)

install(FILES include/wondervol.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package config
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/WonderVolConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_file(cmake/WonderVolConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/WonderVolConfig.cmake"
    @ONLY
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/WonderVolConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/WonderVolConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/WonderVol
)
