# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for HIP requires CMake 3.21.
cmake_minimum_required(VERSION 3.21) # HIP langauge support requires minimum 3.21

find_package(HIPToolkit)

# Enable HIP as a language.
enable_language(HIP)

# Set the HIP build flags.
include(detray-compiler-options-hip)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    detray_test_hip
    STATIC
    "material_validation.hpp"
    "material_validation.hip"
    "navigation_validation.hpp"
    "navigation_validation.hip"
)

add_library(detray::test_hip ALIAS detray_test_hip)

target_link_libraries(
    detray_test_hip
    PUBLIC
        HIP::hiprt
        vecmem::hip
        covfie::hip
        detray::core_array
        detray::test_device
        detray::test_cpu
        detray::validation_utils
        HIP::hiprt
)
