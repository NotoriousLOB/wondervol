cmake_minimum_required(VERSION 3.16)
project(DispersionScanner VERSION 1.0.0 LANGUAGES C)

# Options
option(DISP_BUILD_TESTS "Build tests" ON)
option(DISP_BUILD_BENCHMARKS "Build benchmarks" ON)
option(DISP_BUILD_EXAMPLES "Build examples" ON)
option(DISP_ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(DISP_STRICT_TRAP "Enable strict trapping on violations" OFF)

# Language standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Header-only library
add_library(dispersion_scanner INTERFACE)
add_library(Dispersion::Scanner ALIAS dispersion_scanner)

target_include_directories(dispersion_scanner INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Platform detection and flags
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    set(DISP_IS_ARM64 TRUE)
    target_compile_definitions(dispersion_scanner INTERFACE DISP_TARGET_ARM64=1)
else()
    set(DISP_IS_ARM64 FALSE)
endif()

# Feature test macros for aligned_alloc
target_compile_definitions(dispersion_scanner INTERFACE _GNU_SOURCE=1)

# OpenMP
if(DISP_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        target_link_libraries(dispersion_scanner INTERFACE OpenMP::OpenMP_C)
    endif()
endif()

# Compile options
if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(dispersion_scanner INTERFACE
        -O3
        -Wall
        -Wextra
        -Wpedantic
    )
    
    if(DISP_IS_ARM64)
        target_compile_options(dispersion_scanner INTERFACE -mcpu=native)
    endif()
endif()

# Strict trap
if(DISP_STRICT_TRAP)
    target_compile_definitions(dispersion_scanner INTERFACE DISP_STRICT_TRAP=1)
endif()

# Link math library
target_link_libraries(dispersion_scanner INTERFACE m)

# Testing
if(DISP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Benchmarks
if(DISP_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Examples
if(DISP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Convenience target: run all tests
if(DISP_BUILD_TESTS)
    add_custom_target(check
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS test_bsm test_iv test_dispersion
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS dispersion_scanner
    EXPORT DispersionScannerTargets
)

install(EXPORT DispersionScannerTargets
    FILE DispersionScannerTargets.cmake
    NAMESPACE Dispersion::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DispersionScanner
)

install(FILES include/dispersion_scanner.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package config
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DispersionScannerConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_file(cmake/DispersionScannerConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/DispersionScannerConfig.cmake"
    @ONLY
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DispersionScannerConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DispersionScannerConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DispersionScanner
)
